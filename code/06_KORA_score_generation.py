"""
Generate KORA AtheroBurden scores for external validation.

Input note:
This script starts from an anonymous KORA protein matrix that has already been
renamed to UKB/Olink-style protein feature names. It does not include raw KORA
data cleaning, identifier handling, or protein-name mapping.

Required local input files:
- data/kora_validation_proteins_renamed.csv
  Anonymous KORA protein matrix. Required columns: Participant ID plus protein
  features using the same names as the UKB training panels.
- data/panel_features_arterial_248.csv
- data/panel_features_mechanistic_680.csv
- data/panel_features_genetic_402.csv
- data/panel_features_whole_proteome_2920.csv
  One-column or multi-column files whose column names define the feature order
  used by each CatBoost model.
- models/Artery_enriched_248.cbm
- models/Artherosclerosis_680.cbm
- models/MR_derived_402.cbm
- models/Whole_proteome_2920.cbm
  CatBoost model files generated during model training. These binary model
  files are not intended to be committed to GitHub.

Output:
- results/kora_atheroburden_scores.csv
  Standardized RawFormulaVal scores for the four AtheroBurden panels.
"""

from pathlib import Path

import pandas as pd
from catboost import CatBoostClassifier
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
MODEL_DIR = ROOT / "models"
RESULTS_DIR = ROOT / "results"
RESULTS_DIR.mkdir(exist_ok=True)


PANELS = {
    "arterial_248": {
        "features": DATA_DIR / "panel_features_arterial_248.csv",
        "model": MODEL_DIR / "Artery_enriched_248.cbm",
        "score": "predictions_248",
    },
    "mechanistic_680": {
        "features": DATA_DIR / "panel_features_mechanistic_680.csv",
        "model": MODEL_DIR / "Artherosclerosis_680.cbm",
        "score": "predictions_680",
    },
    "genetic_402": {
        "features": DATA_DIR / "panel_features_genetic_402.csv",
        "model": MODEL_DIR / "MR_derived_402.cbm",
        "score": "predictions_402",
    },
    "whole_proteome_2920": {
        "features": DATA_DIR / "panel_features_whole_proteome_2920.csv",
        "model": MODEL_DIR / "Whole_proteome_2920.cbm",
        "score": "predictions_2920",
    },
}


def read_feature_order(path: Path) -> list[str]:
    feature_table = pd.read_csv(path)
    for column in ("feature", "protein", "name"):
        if column in feature_table.columns:
            features = feature_table[column].dropna().astype(str).tolist()
            return [feature for feature in features if feature != "Participant ID"]
    return [feature for feature in feature_table.columns if feature != "Participant ID"]


def align_features(data: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    aligned = data.copy()
    missing = [feature for feature in features if feature not in aligned.columns]
    for feature in missing:
        aligned[feature] = float("nan")
    return aligned[features]


def score_panel(protein_data: pd.DataFrame, panel_info: dict[str, Path | str]) -> pd.Series:
    features = read_feature_order(panel_info["features"])
    model_input = align_features(protein_data, features)

    model = CatBoostClassifier()
    model.load_model(str(panel_info["model"]))

    raw_score = model.predict(model_input, prediction_type="RawFormulaVal")
    return pd.Series(raw_score, index=protein_data.index, name=panel_info["score"])


def main() -> None:
    protein_data = pd.read_csv(DATA_DIR / "kora_validation_proteins_renamed.csv")
    protein_data = protein_data.set_index("Participant ID")

    scores = pd.concat(
        [score_panel(protein_data, panel_info) for panel_info in PANELS.values()],
        axis=1,
    )

    score_columns = ["predictions_248", "predictions_680", "predictions_402", "predictions_2920"]
    scores[score_columns] = StandardScaler().fit_transform(scores[score_columns])
    scores.to_csv(RESULTS_DIR / "kora_atheroburden_scores.csv", index=True)


if __name__ == "__main__":
    main()
