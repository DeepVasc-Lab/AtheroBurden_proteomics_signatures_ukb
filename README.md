# Proteomic Signatures as Biomarkers of Atherosclerosis Burden

This repository contains the analysis code for our study on proteomic
signatures as biomarkers of atherosclerosis burden. The code supports model
training, score generation, association analyses, prospective validation,
longitudinal analyses, external KORA validation, and manuscript tables/figures.

The repository starts from anonymous analysis-ready files. Raw UK Biobank,
KORA, imaging, registry, and proteomics source data are not included.

## Analysis Workflow

### 1. Protein panel selection

Select the artery-enriched protein panel from GTEx tissue-specific expression
results and select the MR-derived protein panel from CAD-to-protein MR
screening results (`00_panel_selection.R`). The whole-proteome panel used all
2,920 measured proteins without prior feature selection. The
atherosclerosis-related panel was generated from Enrichr atherosclerosis gene
sets mapped to the UKB Olink proteome.

### 2. AtheroBurden model training and score generation

Train CatBoost models for the AtheroBurden panels and generate model scores
(`01_ML_optuna_train_and_score.py`).

### 3. Model performance and Figure 2

Evaluate model performance and generate Figure 2 analyses, including ROC,
confusion matrix, SHAP summary, score density, and vascular-bed burden plots
(`02_Figure2.R`).

### 4. Plaque association analysis

Test associations between AtheroBurden scores and carotid plaque outcomes using
anonymous plaque analysis files (`03_Figure3_plaque_regression.R`).

### 5. UK Biobank prospective analyses

Run Cox models, Kaplan-Meier plots, time-dependent ROC analyses, SCORE2
calculation, C-index, NRI/IDI, and calibration analyses for prospective UKB
outcomes (`04_Figure4_UKB_prospective.R`; `08_SCORE2_calculation.R`;
`09_C_index_NRI_IDI_Calibration.R`).

### 6. Longitudinal AtheroBurden trajectories

Analyze longitudinal changes in AtheroBurden scores across SCORE2 risk groups
and MACE status using mixed-effects models (`05_Figure5_longitudinal.R`).

### 7. KORA external validation

Generate restricted KORA AtheroBurden scores and perform external validation
analyses in KORA S4 and Age1 (`06_KORA_score_generation.py`;
`06_Figure6_KORA_validation.R`).

### 8. Manuscript tables

Generate the main descriptive tables from anonymous analysis-ready files
(`07_Table1_summary.R`).

## Requirements

Analyses were conducted in R and Python. Required packages include common
statistical and machine-learning libraries used in the scripts, including
`data.table`, `dplyr`, `ggplot2`, `survival`, `survminer`, `catboost`,
`scikit-learn`, `pandas`, and related packages.

## Data

Individual-level cohort data and generated output files are not distributed in
this repository. Scripts are intended to be run with locally available
anonymous analysis-ready files.

## Contact

For questions about the analysis or data, please contact:

* Lanyue Zhang - Lanyue.Zhang@med.uni-muenchen.de / zhanglanyue1996@gmail.com
