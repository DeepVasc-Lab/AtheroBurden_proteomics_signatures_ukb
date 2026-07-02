# Protein panel selection for AtheroBurden signatures.
#
# This script contains the panel-selection steps used before model training.
# It starts from public/summary-level or anonymous analysis-ready files and
# does not include individual-level raw data processing.
#
# Notes for the two panels not recalculated here:
# - Whole Proteome Panel: no prior feature selection was applied. All 2,920
#   Olink plasma proteins were used directly for model development.
# - Atherosclerosis-Related Protein Panel: atherosclerosis-related gene sets
#   were collected from Enrichr, combined across relevant terms/pathways, and
#   mapped to the UKB Olink proteome to define the 680-protein panel.
#
# How to run selected sections:
#   Rscript code/00_panel_selection.R artery
#   Rscript code/00_panel_selection.R mr
#
# Artery-enriched workflow:
# 1. Compare GTEx Blood Vessel samples (artery aorta, coronary artery, and
#    tibial artery) with all other GTEx tissues.
# 2. Run differential expression with easyTCGA for count and TPM matrices.
# 3. Intersect GTEx limma results with Olink 2920 protein names.
# 4. Keep proteins with Blood Vessel enrichment greater than three-fold.
# 5. Use the count-based limma hits as the artery-enriched panel.
#
# MR-derived workflow:
# 1. Read CAD GWAS instruments from GCST90132315 and keep p < 5e-8.
# 2. LD-clump the instruments using a European PLINK reference.
# 3. Match clumped SNPs to CARDIoGRAMplusC4D and use those beta/se values as
#    the CAD exposure for MR.
# 4. Run IVW MR for each UKB-PPP protein GWAS outcome file.
# 5. Select proteins using the raw MR p value.
# 6. Export all protein MR results plus raw-p-significant hits.
#
# Section 1: artery-enriched protein selection
# Required input files:
# - data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
#   GTEx v8 gene TPM matrix. Required columns: Name, Description, and samples.
# - data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
#   GTEx v8 gene read-count matrix with the same gene/sample structure.
# - data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
#   GTEx sample annotation table. Required columns: SAMPID and SMTS.
# - data/Olink2920rm_missing.csv
#   Olink 2920 protein matrix. The first column is the participant identifier;
#   the next protein columns provide the Olink feature names to intersect with
#   GTEx differential-expression results.
#
# Section 2: MR-derived protein selection
# Required input files:
# - data/mr_cad_aragam_gwas.tsv.gz
#   CAD GWAS summary statistics from GCST90132315/Aragam et al. Required
#   columns: rsid, beta, standard_error, effect_allele, other_allele,
#   effect_allele_frequency, p_value.
# - data/mr_cardiogramplusc4d_gwas.txt
#   CARDIoGRAMplusC4D CAD GWAS summary statistics for the clumped instruments.
#   Required columns: markername, beta, se_dgc, effect_allele,
#   noneffect_allele, effect_allele_freq, p_dgc.
# - data/UKB_PPP_GWAS_CAD_snp/*_CVD.csv
#   One protein GWAS outcome file per protein. Required columns: SNP, beta, se,
#   effect_allele, other_allele, eaf, and pval.
# - data/plink/EUR.{bed,bim,fam}
#   European LD reference files for local PLINK clumping.

# Local input/output folders. Data files are not committed to GitHub.
DATA_DIR <- "data"
RESULTS_DIR <- "results"

# Selection thresholds used in the manuscript analysis.
ARTERY_FOLD_CHANGE_CUTOFF <- 3
ARTERY_LOG2FC_CUTOFF <- log2(ARTERY_FOLD_CHANGE_CUTOFF)
MR_EXPOSURE_P_CUTOFF <- 5e-8
MR_RAW_P_CUTOFF <- 0.05

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)


require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Please install the R package: ", package, call. = FALSE)
  }
}


check_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing required input files:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
}


read_gtex_gct <- function(path) {
  gct <- tryCatch(
    data.table::fread(path, skip = 2, data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(gct) || !all(c("Name", "Description") %in% names(gct))) {
    gct <- data.table::fread(path, data.table = FALSE)
  }
  if (!all(c("Name", "Description") %in% names(gct))) {
    stop("GTEx matrix must contain Name and Description columns: ", path, call. = FALSE)
  }
  gct
}


prepare_expression_matrix <- function(gct) {
  gct <- dplyr::distinct(gct, .data$Description, .keep_all = TRUE)
  sample_cols <- setdiff(names(gct), c("Name", "Description"))
  expr <- as.data.frame(gct[, sample_cols, drop = FALSE])
  rownames(expr) <- gct$Description
  as.matrix(expr)
}


align_expression_with_group <- function(expr, annotation) {
  sample_table <- data.frame(SAMPID = colnames(expr), stringsAsFactors = FALSE)
  sample_table <- dplyr::inner_join(sample_table, annotation, by = "SAMPID")

  if (nrow(sample_table) == 0) {
    stop("No expression samples matched the GTEx sample annotation table.", call. = FALSE)
  }

  expr <- expr[, sample_table$SAMPID, drop = FALSE]
  group <- factor(sample_table$SMTS, levels = c("Other", "Blood Vessel"))
  list(expr = expr, group = group)
}


normalize_limma_result <- function(result) {
  result <- as.data.frame(result)
  if (!("genesymbol" %in% names(result))) {
    result$genesymbol <- rownames(result)
  }
  result
}


write_feature_list <- function(features, path) {
  data.table::fwrite(
    data.frame(feature = unique(features), stringsAsFactors = FALSE),
    path
  )
}


read_olink_features <- function(path) {
  olink <- data.table::fread(path, nrows = 1, data.table = FALSE)
  columns <- names(olink)
  if ("Participant ID" %in% columns) {
    features <- columns[columns != "Participant ID"]
  } else {
    features <- columns[-1]
  }
  features <- features[seq_len(min(2920, length(features)))]
  features <- setdiff(features, c("group", "label", "burden"))
  data.frame(genesymbol = features, stringsAsFactors = FALSE)
}


select_olink_enriched_genes <- function(limma_result, olink_features) {
  limma_result <- normalize_limma_result(limma_result)
  if (!("logFC" %in% names(limma_result))) {
    stop("The limma result must contain a logFC column.", call. = FALSE)
  }

  dplyr::inner_join(olink_features, limma_result, by = "genesymbol") |>
    dplyr::filter(.data$logFC > ARTERY_LOG2FC_CUTOFF) |>
    dplyr::arrange(dplyr::desc(.data$logFC))
}


run_artery_enriched_selection <- function() {
  require_package("data.table")
  require_package("dplyr")
  require_package("easyTCGA")

  tpm_file <- file.path(DATA_DIR, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
  counts_file <- file.path(DATA_DIR, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
  annotation_file <- file.path(DATA_DIR, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
  olink_file <- file.path(DATA_DIR, "Olink2920rm_missing.csv")

  check_files(c(tpm_file, counts_file, annotation_file, olink_file))

  tpm_matrix <- prepare_expression_matrix(read_gtex_gct(tpm_file))
  counts_matrix <- prepare_expression_matrix(read_gtex_gct(counts_file))

  annotation <- data.table::fread(annotation_file, data.table = FALSE)
  annotation <- annotation[, c("SAMPID", "SMTS")]
  annotation$SMTS <- ifelse(annotation$SMTS == "Blood Vessel", "Blood Vessel", "Other")

  counts_input <- align_expression_with_group(counts_matrix, annotation)
  tpm_input <- align_expression_with_group(tpm_matrix, annotation)

  diff_res_count <- easyTCGA::diff_analysis(
    exprset = counts_input$expr,
    group = counts_input$group,
    is_count = TRUE
  )
  diff_res_tpm <- easyTCGA::diff_analysis(
    exprset = tpm_input$expr,
    group = tpm_input$group,
    is_count = FALSE
  )

  olink_features <- read_olink_features(olink_file)
  count_hits <- select_olink_enriched_genes(diff_res_count[["deg_limma"]], olink_features)
  tpm_hits <- select_olink_enriched_genes(diff_res_tpm[["deg_limma"]], olink_features)

  artery_panel <- count_hits |>
    dplyr::distinct(.data$genesymbol, .keep_all = TRUE)

  data.table::fwrite(
    count_hits,
    file.path(RESULTS_DIR, "panel_selection_artery_enriched_count_limma_fold_change_gt3.csv")
  )
  data.table::fwrite(
    tpm_hits,
    file.path(RESULTS_DIR, "panel_selection_artery_enriched_tpm_limma_fold_change_gt3.csv")
  )
  write_feature_list(
    artery_panel$genesymbol,
    file.path(RESULTS_DIR, "panel_features_arterial_248.csv")
  )
}


run_mr_derived_selection <- function() {
  require_package("data.table")
  require_package("dplyr")
  require_package("tidyr")
  require_package("TwoSampleMR")
  require_package("ieugwasr")
  require_package("genetics.binaRies")

  cad_instrument_gwas_file <- file.path(DATA_DIR, "mr_cad_aragam_gwas.tsv.gz")
  cardiogram_c4d_gwas_file <- file.path(DATA_DIR, "mr_cardiogramplusc4d_gwas.txt")
  protein_outcome_dir <- file.path(DATA_DIR, "UKB_PPP_GWAS_CAD_snp")
  plink_bfile <- file.path(DATA_DIR, "plink", "EUR")

  check_files(c(cad_instrument_gwas_file, cardiogram_c4d_gwas_file))
  if (!dir.exists(protein_outcome_dir)) {
    stop("Missing protein outcome directory: ", protein_outcome_dir, call. = FALSE)
  }
  check_files(paste0(plink_bfile, c(".bed", ".bim", ".fam")))

  cad_instrument_gwas <- TwoSampleMR::read_exposure_data(
    filename = cad_instrument_gwas_file,
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )

  genome_wide_hits <- dplyr::filter(cad_instrument_gwas, .data$pval.exposure < MR_EXPOSURE_P_CUTOFF)
  if (nrow(genome_wide_hits) == 0) {
    stop("No genome-wide significant exposure variants were found.", call. = FALSE)
  }

  clumped <- ieugwasr::ld_clump(
    dplyr::tibble(rsid = genome_wide_hits$SNP, pval = genome_wide_hits$pval.exposure),
    clump_r2 = 0.001,
    clump_p = 1,
    clump_kb = 10000,
    bfile = plink_bfile,
    plink_bin = genetics.binaRies::get_plink_binary()
  )

  cardiogram_c4d_gwas <- TwoSampleMR::read_exposure_data(
    filename = cardiogram_c4d_gwas_file,
    sep = "\t",
    snp_col = "markername",
    beta_col = "beta",
    se_col = "se_dgc",
    effect_allele_col = "effect_allele",
    other_allele_col = "noneffect_allele",
    eaf_col = "effect_allele_freq",
    pval_col = "p_dgc"
  )

  exposure_data <- dplyr::filter(cardiogram_c4d_gwas, .data$SNP %in% clumped$rsid)
  exposure_data$exposure <- "CAD"

  data.table::fwrite(exposure_data, file.path(RESULTS_DIR, "mr_cardiogramplusc4d_matched_exposure.csv"))

  protein_files <- list.files(
    path = protein_outcome_dir,
    pattern = "_CVD\\.csv$",
    full.names = TRUE
  )
  if (length(protein_files) == 0) {
    stop("No *_CVD.csv protein outcome files were found in: ", protein_outcome_dir, call. = FALSE)
  }

  process_protein_file <- function(file) {
    tryCatch({
      outcome <- data.table::fread(file, data.table = FALSE)
      outcome_data <- TwoSampleMR::format_data(outcome, type = "outcome")
      harmonised <- TwoSampleMR::harmonise_data(exposure_data, outcome_data, action = 1)
      mr_result <- TwoSampleMR::mr(harmonised, method_list = c("mr_ivw"))

      if (nrow(mr_result) == 0) {
        return(NULL)
      }

      mr_result$feature_id <- sub("_CVD\\.csv$", "", basename(file))
      mr_result
    }, error = function(e) {
      message("Skipping ", basename(file), ": ", e$message)
      NULL
    })
  }

  mr_results <- lapply(protein_files, process_protein_file)
  mr_results <- dplyr::bind_rows(mr_results)

  if (nrow(mr_results) == 0) {
    stop("No valid MR results were generated.", call. = FALSE)
  }

  result_df <- mr_results |>
    dplyr::rename(P_Value = pval) |>
    tidyr::separate(
      feature_id,
      into = c("Protein", "OID"),
      sep = "_",
      remove = FALSE,
      extra = "merge",
      fill = "right"
    ) |>
    dplyr::mutate(padjust = p.adjust(.data$P_Value, method = "fdr")) |>
    dplyr::arrange(.data$P_Value)

  raw_p_hits <- dplyr::filter(result_df, .data$P_Value < MR_RAW_P_CUTOFF)
  fdr_hits <- dplyr::filter(result_df, .data$padjust < 0.05)

  data.table::fwrite(result_df, file.path(RESULTS_DIR, "MRderived_all_results.csv"))
  data.table::fwrite(raw_p_hits, file.path(RESULTS_DIR, "MRderived_list_rawPsignificant.csv"))
  write_feature_list(raw_p_hits$Protein, file.path(RESULTS_DIR, "panel_features_genetic_402.csv"))
  data.table::fwrite(fdr_hits, file.path(RESULTS_DIR, "MRderived_list_FDRsignificant_reference.csv"))
}


args <- commandArgs(trailingOnly = TRUE)
allowed_args <- c("artery", "mr")
if (length(args) != 1 || !(args %in% allowed_args)) {
  stop(
    "Run one section at a time: Rscript code/00_panel_selection.R artery OR Rscript code/00_panel_selection.R mr",
    call. = FALSE
  )
}

if (args == "artery") {
  run_artery_enriched_selection()
}
if (args == "mr") {
  run_mr_derived_selection()
}
