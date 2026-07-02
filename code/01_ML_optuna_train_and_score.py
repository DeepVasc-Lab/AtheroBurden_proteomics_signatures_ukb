#!/usr/bin/env python3
"""
ML workflow for AtheroBurden signatures.

What it does:
1. Read one matched discovery dataset:
   Participant ID | protein columns | group
2. Split into train/test: 80/20, stratified, random_state=42.
3. On the training set, tune each ML classifier with Optuna using
   5-fold x 10-repeat stratified CV.
4. Evaluate all tuned classifiers on the held-out test set.
5. Refit CatBoost on the training set and export CatBoost raw/z scores.

Final downstream score:
  CatBoost RawFormulaVal, not probability.

How to run:
  cd atheroburden-code-release

  # Arterial-enriched panel, 248 proteins
  python code/01_ML_optuna_train_and_score.py \
    --input-csv 0720_olink248.csv \
    --panel-name arterial_248 \
    --output-dir results/arterial_248_optuna \
    --n-features 248

  # MR-derived/genetic panel, 402 proteins
  python code/01_ML_optuna_train_and_score.py \
    --input-csv data/ml_genetic_402.csv \
    --panel-name genetic_402 \
    --output-dir results/genetic_402_optuna \
    --n-features 402

  # Atherosclerosis-related/mechanistic panel, 680 proteins
  python code/01_ML_optuna_train_and_score.py \
    --input-csv 0720_olink680.csv \
    --panel-name mechanistic_680 \
    --output-dir results/mechanistic_680_optuna \
    --n-features 680

  # Whole proteome panel, 2920 proteins
  python code/01_ML_optuna_train_and_score.py \
    --input-csv 0720_olink2920.csv \
    --panel-name whole_2920 \
    --output-dir results/whole_2920_optuna \
    --n-features 2920
  
  # If you want to score a broader UKB matrix after training:
  add: --score-csv /path/to/OlinkNDimpute.csv

Main outputs:
  cv_summary_by_model.csv          classifier comparison from repeated CV
  cv_results_each_fold.csv         fold-level CV metrics
  test_results_all_models.csv      held-out test metrics for all classifiers
  catboost_scores_all_samples.csv  final CatBoost raw/z scores only
  CatBoost_model.cbm               final CatBoost model
  manifest.json                    compact run record
"""

import argparse
import json
import warnings
from pathlib import Path

import joblib
import numpy as np
import optuna
import pandas as pd
from catboost import CatBoostClassifier
from lightgbm import LGBMClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from xgboost import XGBClassifier

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=ConvergenceWarning)


MODELS = [
    "SGDClassifier",
    "Random Forest",
    "Logistic Regression",
    "ElasticNET",
    "SVM",
    "MLP",
    "LGBM",
    "CatBoost",
    "XGBoost",
]

METRICS = ["Accuracy", "Precision", "Recall", "F1", "ROC_AUC"]


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input-csv", required=True)
    p.add_argument("--panel-name", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--score-csv", default=None, help="Optional broader UKB matrix to score.")
    p.add_argument("--id-col", default="Participant ID")
    p.add_argument("--label-col", default=None, help="Default: last column.")
    p.add_argument("--n-features", type=int, default=None, help="Default: all columns except ID/label.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n-trials", type=int, default=120)
    p.add_argument("--n-jobs", type=int, default=4)
    p.add_argument(
        "--primary-metric",
        default="ROC_AUC",
        choices=METRICS,
        help="Use Accuracy if you want to reproduce the earliest scripts exactly.",
    )
    return p.parse_args()


def load_data(path, id_col, label_col, n_features):
    df = pd.read_csv(path)
    df[id_col] = df[id_col].astype(str)

    if label_col is None:
        label_col = [c for c in df.columns if c != id_col][-1]

    feature_cols = [c for c in df.columns if c not in [id_col, label_col]]
    if n_features is not None:
        feature_cols = feature_cols[:n_features]

    x = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    y = pd.to_numeric(df[label_col], errors="raise").astype(int).to_numpy()
    if sorted(np.unique(y).tolist()) != [0, 1]:
        raise ValueError("The label must be binary 0/1.")

    return df, x, y, feature_cols, label_col


def params_for_trial(trial, model_name, seed, n_jobs):
    if model_name == "SGDClassifier":
        return {
            "alpha": trial.suggest_float("alpha", 1e-5, 1e-1, log=True),
            "penalty": trial.suggest_categorical("penalty", ["l1", "l2", "elasticnet"]),
        }

    if model_name == "Random Forest":
        return {
            "n_estimators": trial.suggest_int("n_estimators", 50, 200),
            "max_depth": trial.suggest_int("max_depth", 3, 10),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 1, 10),
            "random_state": seed,
            "n_jobs": n_jobs,
        }

    if model_name == "Logistic Regression":
        return {
            "C": trial.suggest_float("C", 1e-5, 10.0, log=True),
            "penalty": "l1",
            "solver": "liblinear",
            "max_iter": 5000,
            "random_state": seed,
        }

    if model_name == "ElasticNET":
        return {
            "C": trial.suggest_float("C", 1e-5, 10.0, log=True),
            "l1_ratio": trial.suggest_float("l1_ratio", 0.05, 0.95),
            "penalty": "elasticnet",
            "solver": "saga",
            "max_iter": 5000,
            "random_state": seed,
        }

    if model_name == "SVM":
        return {
            "C": trial.suggest_float("C", 1e-5, 10.0, log=True),
            "kernel": "linear",
            "probability": True,
            "random_state": seed,
        }

    if model_name == "MLP":
        hidden = trial.suggest_categorical("hidden_layer_sizes", ["64", "128", "64,32"])
        return {
            "hidden_layer_sizes": tuple(int(v) for v in hidden.split(",")),
            "activation": trial.suggest_categorical("activation", ["relu", "logistic", "tanh"]),
            "alpha": trial.suggest_float("alpha", 1e-5, 1e-1, log=True),
            "learning_rate_init": trial.suggest_float("learning_rate_init", 1e-5, 1e-2, log=True),
            "max_iter": 1000,
            "early_stopping": True,
            "random_state": seed,
        }

    if model_name == "LGBM":
        return {
            "objective": "binary",
            "metric": "auc",
            "n_estimators": trial.suggest_int("n_estimators", 50, 200),
            "learning_rate": trial.suggest_float("learning_rate", 1e-5, 1e-1, log=True),
            "max_depth": trial.suggest_int("max_depth", 3, 10),
            "num_leaves": trial.suggest_int("num_leaves", 7, 63),
            "verbosity": -1,
            "random_state": seed,
            "n_jobs": n_jobs,
        }

    if model_name == "CatBoost":
        return {
            "loss_function": "Logloss",
            "eval_metric": "AUC",
            "iterations": trial.suggest_int("iterations", 50, 200),
            "learning_rate": trial.suggest_float("learning_rate", 1e-5, 1e-1, log=True),
            "depth": trial.suggest_int("depth", 3, 10),
            "l2_leaf_reg": trial.suggest_float("l2_leaf_reg", 1e-3, 20.0, log=True),
            "bootstrap_type": trial.suggest_categorical("bootstrap_type", ["Bayesian", "Bernoulli", "No"]),
            "random_seed": seed,
            "thread_count": n_jobs,
            "verbose": False,
            "allow_writing_files": False,
        }

    if model_name == "XGBoost":
        return {
            "n_estimators": trial.suggest_int("n_estimators", 50, 200),
            "learning_rate": trial.suggest_float("learning_rate", 1e-5, 1e-1, log=True),
            "max_depth": trial.suggest_int("max_depth", 3, 10),
            "subsample": trial.suggest_float("subsample", 0.60, 1.00),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.60, 1.00),
            "eval_metric": "logloss",
            "random_state": seed,
            "n_jobs": n_jobs,
        }

    raise ValueError(model_name)


def build_model(model_name, params):
    if model_name == "SGDClassifier":
        return make_pipeline(
            SimpleImputer(strategy="median"),
            StandardScaler(),
            SGDClassifier(loss="log_loss", max_iter=10000, learning_rate="optimal", **params),
        )
    if model_name == "Random Forest":
        return make_pipeline(SimpleImputer(strategy="median"), RandomForestClassifier(**params))
    if model_name == "Logistic Regression":
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(**params))
    if model_name == "ElasticNET":
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(**params))
    if model_name == "SVM":
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), SVC(**params))
    if model_name == "MLP":
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), MLPClassifier(**params))
    if model_name == "LGBM":
        return LGBMClassifier(**params)
    if model_name == "CatBoost":
        return CatBoostClassifier(**params)
    if model_name == "XGBoost":
        return XGBClassifier(**params)
    raise ValueError(model_name)


def prediction_score(model, x):
    if hasattr(model, "predict_proba"):
        return model.predict_proba(x)[:, 1]
    return model.decision_function(x)


def calc_metrics(y_true, y_pred, score):
    return {
        "Accuracy": accuracy_score(y_true, y_pred),
        "Precision": precision_score(y_true, y_pred, zero_division=0),
        "Recall": recall_score(y_true, y_pred, zero_division=0),
        "F1": f1_score(y_true, y_pred, zero_division=0),
        "ROC_AUC": roc_auc_score(y_true, score),
    }


def cv_score(model_name, params, x_train, y_train, folds, primary_metric):
    rows = []
    for fold_id, (i_train, i_valid) in enumerate(folds, start=1):
        model = build_model(model_name, params)
        model.fit(x_train.iloc[i_train], y_train[i_train])
        pred = model.predict(x_train.iloc[i_valid])
        score = prediction_score(model, x_train.iloc[i_valid])
        rows.append({"Fold": fold_id, **calc_metrics(y_train[i_valid], pred, score)})
    cv = pd.DataFrame(rows)
    return float(cv[primary_metric].mean()), cv


def tune_model(model_name, x_train, y_train, folds, args):
    def objective(trial):
        params = params_for_trial(trial, model_name, args.seed, args.n_jobs)
        mean_score, _ = cv_score(model_name, params, x_train, y_train, folds, args.primary_metric)
        return mean_score

    study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=args.seed))
    study.optimize(objective, n_trials=args.n_trials, show_progress_bar=True)
    best_params = params_for_trial_from_best(study.best_trial, model_name, args.seed, args.n_jobs)
    _, cv = cv_score(model_name, best_params, x_train, y_train, folds, args.primary_metric)
    return best_params, cv, study.trials_dataframe()


def params_for_trial_from_best(best_trial, model_name, seed, n_jobs):
    params = dict(best_trial.params)

    if model_name == "SGDClassifier":
        return params
    if model_name == "Random Forest":
        params.update({"random_state": seed, "n_jobs": n_jobs})
        return params
    if model_name == "Logistic Regression":
        params.update({"penalty": "l1", "solver": "liblinear", "max_iter": 5000, "random_state": seed})
        return params
    if model_name == "ElasticNET":
        params.update({"penalty": "elasticnet", "solver": "saga", "max_iter": 5000, "random_state": seed})
        return params
    if model_name == "SVM":
        params.update({"kernel": "linear", "probability": True, "random_state": seed})
        return params
    if model_name == "MLP":
        params["hidden_layer_sizes"] = tuple(int(v) for v in str(params["hidden_layer_sizes"]).split(","))
        params.update({"max_iter": 1000, "early_stopping": True, "random_state": seed})
        return params
    if model_name == "LGBM":
        params.update({"objective": "binary", "metric": "auc", "verbosity": -1, "random_state": seed, "n_jobs": n_jobs})
        return params
    if model_name == "CatBoost":
        params.update({
            "loss_function": "Logloss",
            "eval_metric": "AUC",
            "random_seed": seed,
            "thread_count": n_jobs,
            "verbose": False,
            "allow_writing_files": False,
        })
        return params
    if model_name == "XGBoost":
        params.update({"eval_metric": "logloss", "random_state": seed, "n_jobs": n_jobs})
        return params

    raise ValueError(model_name)


def catboost_raw(model, x):
    return np.asarray(model.predict(x, prediction_type="RawFormulaVal")).reshape(-1)


def main():
    args = get_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    df, x, y, features, label_col = load_data(args.input_csv, args.id_col, args.label_col, args.n_features)

    train_idx, test_idx = train_test_split(
        np.arange(len(y)),
        test_size=0.20,
        random_state=args.seed,
        stratify=y,
    )
    x_train, x_test = x.iloc[train_idx].reset_index(drop=True), x.iloc[test_idx].reset_index(drop=True)
    y_train, y_test = y[train_idx], y[test_idx]

    folds = list(RepeatedStratifiedKFold(
        n_splits=5,
        n_repeats=10,
        random_state=args.seed,
    ).split(x_train, y_train))

    cv_all = []
    trial_all = []
    best_params = {}
    cv_summary = []

    for model_name in MODELS:
        print(f"\nTuning {model_name}")
        params, cv, trials = tune_model(model_name, x_train, y_train, folds, args)
        best_params[model_name] = params

        cv.insert(0, "Model", model_name)
        cv_all.append(cv)

        trials.insert(0, "Model", model_name)
        trial_all.append(trials)

        row = {"Model": model_name}
        for metric in METRICS:
            row[f"CV_mean_{metric}"] = cv[metric].mean()
            row[f"CV_sd_{metric}"] = cv[metric].std(ddof=1)
        row["Best_params_json"] = json.dumps(params, sort_keys=True)
        cv_summary.append(row)

    cv_summary = pd.DataFrame(cv_summary).sort_values(f"CV_mean_{args.primary_metric}", ascending=False)
    cv_summary.to_csv(out / "cv_summary_by_model.csv", index=False)
    pd.concat(cv_all, ignore_index=True).to_csv(out / "cv_results_each_fold.csv", index=False)
    pd.concat(trial_all, ignore_index=True).to_csv(out / "optuna_trials.csv", index=False)

    test_rows = []
    catboost_model = None
    for model_name in MODELS:
        model = build_model(model_name, best_params[model_name])
        model.fit(x_train, y_train)
        pred = model.predict(x_test)
        score = prediction_score(model, x_test)
        test_rows.append({"Model": model_name, **calc_metrics(y_test, pred, score)})
        if model_name == "CatBoost":
            catboost_model = model

    pd.DataFrame(test_rows).sort_values(args.primary_metric, ascending=False).to_csv(
        out / "test_results_all_models.csv",
        index=False,
    )

    train_raw = catboost_raw(catboost_model, x_train)
    train_mean = train_raw.mean()
    train_sd = train_raw.std(ddof=1)

    score_df = pd.read_csv(args.score_csv) if args.score_csv else df
    score_df[args.id_col] = score_df[args.id_col].astype(str)
    score_x = score_df[features].apply(pd.to_numeric, errors="coerce")
    raw = catboost_raw(catboost_model, score_x)

    scores = pd.DataFrame({
        args.id_col: score_df[args.id_col],
        f"{args.panel_name}_catboost_raw": raw,
        f"{args.panel_name}_catboost_z": (raw - train_mean) / train_sd,
    })
    scores.to_csv(out / "catboost_scores_all_samples.csv", index=False)

    catboost_model.save_model(str(out / "CatBoost_model.cbm"))
    joblib.dump(catboost_model, out / "CatBoost_model.joblib")

    pd.DataFrame({"feature": features, "order": np.arange(1, len(features) + 1)}).to_csv(
        out / "feature_order.csv",
        index=False,
    )

    manifest = {
        "input_csv": args.input_csv,
        "score_csv": args.score_csv,
        "panel_name": args.panel_name,
        "label_col": label_col,
        "seed": args.seed,
        "train_test_split": "80/20 stratified",
        "cv": "5-fold stratified CV repeated 10 times on the training set",
        "n_trials": args.n_trials,
        "primary_metric": args.primary_metric,
        "final_score": "CatBoost RawFormulaVal; z-score uses train-set raw mean/SD",
        "n_samples": int(len(y)),
        "n_train": int(len(y_train)),
        "n_test": int(len(y_test)),
        "n_features": int(len(features)),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("\nDone.")
    print("CV summary:", out / "cv_summary_by_model.csv")
    print("Test results:", out / "test_results_all_models.csv")
    print("CatBoost scores:", out / "catboost_scores_all_samples.csv")


if __name__ == "__main__":
    main()
