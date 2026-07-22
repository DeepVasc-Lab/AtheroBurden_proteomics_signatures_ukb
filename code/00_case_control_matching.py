#!/usr/bin/env python3
"""Define prevalent atherosclerosis cases and build a healthy 1:1 control set."""

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


INPUT_FILES = {
    "participants": "participant_baseline.csv",
    "proteomics": "olink_npx.csv",
    "episodes": "hospital_episodes.csv",
    "diagnoses": "hospital_diagnoses.csv",
    "operations": "hospital_operations.csv",
}


# Periods are removed and codes are matched by prefix.
ATHEROSCLEROSIS_CODES = {
    "coronary": {
        "icd9": ("4109", "4119", "4129", "4139", "4140", "4148", "4149"),
        "icd10": (
            "I20", "I21", "I22", "I23", "I24",
            "I250", "I251", "I252", "I255", "I256", "I258", "I259",
        ),
        "opcs4": (
            "K40", "K41", "K42", "K43", "K44", "K45", "K46", "K471",
            "K49", "K501", "K502", "K504", "K75",
        ),
    },
    "cerebrovascular": {
        "icd9": ("4331", "4339", "4349", "4359", "4370", "4371"),
        "icd10": (
            "I630", "I631", "I632", "I633", "I634", "I635", "I638", "I639",
            "I65", "I66", "I672",
        ),
        "opcs4": ("L29", "L301", "L302", "L303", "L311", "L313", "L314", "L343"),
    },
    "aortic": {
        "icd9": ("4400", "441", "4440"),
        "icd10": ("I700", "I713", "I714"),
        "opcs4": (
            "L251", "L252", "L253", "L254", "L261", "L262", "L263", "L265",
            "L266", "L267", "L27", "L28", "L45", "L461", "L471",
        ),
    },
    "peripheral": {
        "icd9": ("4402", "4442"),
        "icd10": ("I702",),
        "opcs4": (),
    },
    "other_artery": {
        "icd9": ("4401", "4408", "4409", "4448"),
        "icd10": ("I701", "I708", "I709"),
        "opcs4": (
            "L383", "L391", "L392", "L393", "L41", "L421", "L431", "L432",
            "L435", "L37", "L38", "L395", "L48", "L49", "L50", "L51",
            "L52", "L53", "L541", "L542", "L544", "L60", "L621", "L622",
            "L631", "L632", "L635", "L661", "L662", "L665", "L667", "L681",
            "L682", "L701",
        ),
    },
}


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Folder containing the five descriptively named input CSV files.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Default: <data-dir>/atherosclerosis_psm_1to1",
    )
    return parser.parse_args()


def clean_icd10(values):
    return (
        values.fillna("")
        .astype(str)
        .str.strip()
        .str.split()
        .str[0]
        .fillna("")
        .str.upper()
        .str.replace(".", "", regex=False)
    )


def build_eligible_participants(data_dir):
    olink = pd.read_csv(data_dir / INPUT_FILES["proteomics"], dtype={"eid": str})
    protein_columns = [column for column in olink.columns if column != "eid"]
    protein_missing = olink[protein_columns].isna().mean()
    kept_proteins = protein_missing[protein_missing <= 0.30].index
    participant_missing = olink[kept_proteins].isna().mean(axis=1)
    eligible_ids = set(olink.loc[participant_missing <= 0.30, "eid"].astype(str))

    columns = [
        "Participant ID",
        "Sex",
        "Age at recruitment",
        "Age when attended assessment centre | Instance 0",
        "Date of attending assessment centre | Instance 0",
    ]
    people = pd.read_csv(data_dir / INPUT_FILES["participants"], usecols=columns, dtype=str)
    people = people.rename(
        columns={
            "Participant ID": "eid",
            "Sex": "sex",
            "Age at recruitment": "age",
            "Age when attended assessment centre | Instance 0": "assessment_age",
            "Date of attending assessment centre | Instance 0": "baseline_date",
        }
    )
    people["age"] = pd.to_numeric(people["age"], errors="coerce")
    assessment_age = pd.to_numeric(people["assessment_age"], errors="coerce")
    people["age"] = people["age"].fillna(assessment_age)
    people["sex"] = people["sex"].fillna("").astype(str).str.strip()
    people["baseline_date"] = pd.to_datetime(people["baseline_date"], errors="coerce")
    people = people[people["eid"].isin(eligible_ids)]
    people = people.dropna(subset=["age", "baseline_date"])
    people = people[people["sex"].ne("")]

    print(f"Olink proteins retained after 30% QC: {len(kept_proteins):,}")
    print(f"Olink participants retained after 30% QC: {len(people):,}")
    return people[["eid", "age", "sex", "baseline_date"]].copy()


def classify_vascular_bed(codes, system):
    bed = pd.Series("", index=codes.index, dtype=object)
    for name, code_sets in ATHEROSCLEROSIS_CODES.items():
        prefixes = code_sets[system]
        if prefixes:
            bed.loc[codes.str.startswith(prefixes)] = name
    return bed


def build_cases_and_icd10_burden(people, data_dir):
    hesin = pd.read_csv(
        data_dir / INPUT_FILES["episodes"],
        usecols=["dnx_hesin_id", "epistart"],
        dtype=str,
    )
    hesin["epistart"] = pd.to_datetime(hesin["epistart"], errors="coerce")
    episode_date = (
        hesin.drop_duplicates("dnx_hesin_id")
        .set_index("dnx_hesin_id")["epistart"]
    )

    eligible = set(people["eid"])
    baseline = people.set_index("eid")["baseline_date"]
    entry_count = defaultdict(int)
    case_events = []

    columns = ["dnx_hesin_id", "eid", "level", "diag_icd9", "diag_icd10"]
    for chunk in pd.read_csv(
        data_dir / INPUT_FILES["diagnoses"],
        usecols=columns,
        dtype=str,
        chunksize=500_000,
    ):
        chunk["eid"] = chunk["eid"].astype(str).str.strip()
        chunk = chunk[chunk["eid"].isin(eligible)]
        chunk = chunk[pd.to_numeric(chunk["level"], errors="coerce").eq(1)]
        if chunk.empty:
            continue

        icd10 = clean_icd10(chunk["diag_icd10"])
        has_icd10 = icd10.ne("")
        for eid, count in chunk.loc[has_icd10].groupby("eid").size().items():
            entry_count[eid] += int(count)

        bed10 = classify_vascular_bed(icd10, "icd10")
        mask10 = bed10.ne("")
        if mask10.any():
            events10 = chunk.loc[mask10, ["eid", "dnx_hesin_id"]].copy()
            events10["code"] = icd10.loc[mask10]
            events10["vascular_bed"] = bed10.loc[mask10]
            events10["event_date"] = events10["dnx_hesin_id"].map(episode_date)
            events10["baseline_date"] = events10["eid"].map(baseline)
            events10 = events10[
                events10["event_date"].notna()
                & events10["event_date"].le(events10["baseline_date"])
            ]
            events10["source"] = "ICD10"
            case_events.append(events10[["eid", "event_date", "source", "code", "vascular_bed"]])

        icd9 = clean_icd10(chunk["diag_icd9"])
        bed9 = classify_vascular_bed(icd9, "icd9")
        mask9 = bed9.ne("")
        if mask9.any():
            # ICD-9 evidence is treated as pre-baseline when its date is unavailable.
            events9 = chunk.loc[mask9, ["eid"]].copy()
            events9["event_date"] = pd.NaT
            events9["source"] = "ICD9"
            events9["code"] = icd9.loc[mask9]
            events9["vascular_bed"] = bed9.loc[mask9]
            case_events.append(events9)

    for chunk in pd.read_csv(
        data_dir / INPUT_FILES["operations"],
        usecols=["eid", "opdate", "oper4"],
        dtype=str,
        chunksize=500_000,
    ):
        chunk["eid"] = chunk["eid"].astype(str).str.strip()
        chunk = chunk[chunk["eid"].isin(eligible)]
        opcs4 = clean_icd10(chunk["oper4"])
        bed = classify_vascular_bed(opcs4, "opcs4")
        mask = bed.ne("")
        if not mask.any():
            continue
        events = chunk.loc[mask, ["eid", "opdate"]].copy()
        events["event_date"] = pd.to_datetime(events["opdate"], errors="coerce")
        events["baseline_date"] = events["eid"].map(baseline)
        events = events[
            events["event_date"].notna()
            & events["event_date"].le(events["baseline_date"])
        ]
        events["source"] = "OPCS4"
        events["code"] = opcs4.loc[events.index]
        events["vascular_bed"] = bed.loc[events.index]
        case_events.append(events[["eid", "event_date", "source", "code", "vascular_bed"]])

    people["icd10_main_alltime_count"] = people["eid"].map(entry_count).fillna(0).astype(int)
    if not case_events:
        raise ValueError("No prevalent atherosclerosis cases were found.")
    events = pd.concat(case_events, ignore_index=True)
    events["source_code"] = events["source"] + ":" + events["code"]

    defining_icd10 = (
        events[events["source"].eq("ICD10")]
        .groupby("eid", as_index=False)
        .agg(
            first_icd10_date=("event_date", "min"),
            defining_icd10_codes=("code", lambda x: ";".join(sorted(set(x)))),
        )
    )
    case_ids = set(defining_icd10["eid"])
    auxiliary_evidence = events[events["eid"].isin(case_ids)]
    cases = (
        auxiliary_evidence.groupby("eid", as_index=False)
        .agg(
            vascular_beds=("vascular_bed", lambda x: ";".join(sorted(set(x)))),
            sources=("source", lambda x: ";".join(sorted(set(x)))),
            codes_used=("source_code", lambda x: ";".join(sorted(set(x)))),
        )
        .merge(defining_icd10, on="eid", how="inner")
        .merge(people[["eid", "age", "sex", "baseline_date"]], on="eid", how="left")
    )
    people["case"] = people["eid"].isin(cases["eid"]).astype(int)
    return people, cases


def add_propensity_score(people):
    age = people[["age"]].astype(float).reset_index(drop=True)
    sex = pd.get_dummies(people["sex"], prefix="sex", dtype=float).reset_index(drop=True)
    model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    people = people.reset_index(drop=True)
    model.fit(pd.concat([age, sex], axis=1), people["case"])
    people["propensity_score"] = model.predict_proba(pd.concat([age, sex], axis=1))[:, 1]
    return people


def match_controls(people):
    candidate_rows = []
    pair_number = 0

    for sex_value in sorted(people["sex"].unique()):
        cases = people[(people["sex"] == sex_value) & (people["case"] == 1)].copy()
        controls = people[(people["sex"] == sex_value) & (people["case"] == 0)].copy()
        cases = cases.sort_values(["propensity_score", "eid"], ascending=[False, True])
        controls = controls.sort_values(["propensity_score", "eid"]).reset_index(drop=True)
        available = np.ones(len(controls), dtype=bool)

        for case in cases.itertuples(index=False):
            available_index = np.flatnonzero(available)
            if len(available_index) < 6:
                continue

            pool = controls.loc[available_index].copy()
            pool["ps_diff"] = (pool["propensity_score"] - case.propensity_score).abs()
            pool["age_diff"] = (pool["age"] - case.age).abs()
            nearest = pool.sort_values(["ps_diff", "age_diff", "eid"]).head(6)
            selected_eid = str(
                nearest.sort_values(
                    [
                        "icd10_main_alltime_count",
                        "ps_diff",
                        "age_diff",
                        "eid",
                    ]
                ).iloc[0]["eid"]
            )

            pair_number += 1
            pair_id = f"pair_{pair_number:06d}"
            for rank, control in enumerate(nearest.itertuples(index=True), start=1):
                candidate_rows.append(
                    {
                        "pair_id": pair_id,
                        "case_eid": case.eid,
                        "control_eid": control.eid,
                        "candidate_rank": rank,
                        "selected": int(control.eid == selected_eid),
                        "case_age": case.age,
                        "control_age": control.age,
                        "sex": case.sex,
                        "ps_diff": control.ps_diff,
                        "control_icd10_main_alltime_count": control.icd10_main_alltime_count,
                    }
                )
            available[nearest.index] = False

    candidates = pd.DataFrame(candidate_rows)
    pairs = candidates[candidates["selected"] == 1].copy()
    return candidates, pairs


def main():
    args = get_args()
    data_dir = Path(args.data_dir)
    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else data_dir / "atherosclerosis_psm_1to1"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    people = build_eligible_participants(data_dir)
    people, cases = build_cases_and_icd10_burden(people, data_dir)
    people = add_propensity_score(people)
    candidates, pairs = match_controls(people)
    pairs = pairs.merge(
        cases[
            [
                "eid",
                "first_icd10_date",
                "defining_icd10_codes",
                "vascular_beds",
                "sources",
                "codes_used",
            ]
        ],
        left_on="case_eid",
        right_on="eid",
        how="left",
    ).drop(columns="eid")

    cases.to_csv(output_dir / "atherosclerosis_cases.csv", index=False)
    candidates.to_csv(output_dir / "psm_1to6_candidates.csv", index=False)
    pairs.to_csv(output_dir / "matched_1to1_pairs.csv", index=False)

    print(f"Prevalent atherosclerosis cases: {len(cases):,}")
    print(f"Matched 1:1 pairs: {len(pairs):,}")
    print(f"Candidate controls: {len(candidates):,}")


if __name__ == "__main__":
    main()
