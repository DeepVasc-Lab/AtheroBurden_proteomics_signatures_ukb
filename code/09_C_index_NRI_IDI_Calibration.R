# Input note:
# This script starts from an anonymous UKB prospective analysis-ready file:
# data/ukb_prospective_analysis.csv
#
# Required columns:
# - MACE1_days, MACE1_event: follow-up time and event indicator for MACE.
# - SCORE2_linear_predictor: baseline SCORE2 model linear predictor.
# - Artery_enriched_248, Atherosclerosis_680, MR_derived_402,
#   Whole_proteome_2920: AtheroBurden score columns from the four CatBoost panels.
#
# Analyses included:
# - C-index comparison: SCORE2 versus SCORE2 plus each AtheroBurden score.
# - NRI and IDI at 10 years.
# - 10-year calibration plots for SCORE2 and SCORE2-plus-AtheroBurden models.

library(data.table)
library(dplyr)
library(tibble)
library(survival)
library(survcomp)
library(nricens)
library(survIDINRI)
library(riskRegression)
library(survival.calib)
library(tidyr)
library(ggplot2)
library(patchwork)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

df <- fread("data/ukb_prospective_analysis.csv", data.table = FALSE)

horizon_days <- 3652.5
risk_cuts <- c(0.075, 0.10)

score_info <- tibble::tribble(
  ~score_var, ~model_name, ~label,
  "Artery_enriched_248", "arterial_248", "SCORE2 + AtheroBurden Arterial Signature",
  "Atherosclerosis_680", "mechanistic_680", "SCORE2 + AtheroBurden Mechanistic Signature",
  "MR_derived_402", "genetic_402", "SCORE2 + AtheroBurden Genetic Signature",
  "Whole_proteome_2920", "whole_proteome_2920", "SCORE2 + AtheroBurden WholeProteome Signature"
)

make_complete_data <- function(score_var = NULL) {
  required_vars <- c("MACE1_days", "MACE1_event", "SCORE2_linear_predictor", score_var)
  df %>%
    mutate(
      MACE1_days = suppressWarnings(as.numeric(MACE1_days)),
      MACE1_event = as.integer(MACE1_event)
    ) %>%
    filter(!is.na(MACE1_days) & MACE1_days > 0, !is.na(MACE1_event)) %>%
    filter(if_all(all_of(required_vars), ~ !is.na(.x)))
}

fit_score2_model <- function(data) {
  coxph(
    Surv(MACE1_days, MACE1_event) ~ SCORE2_linear_predictor,
    data = data,
    x = TRUE,
    y = TRUE
  )
}

fit_score2_plus_model <- function(data, score_var) {
  fml <- as.formula(
    paste0("Surv(MACE1_days, MACE1_event) ~ SCORE2_linear_predictor + ", score_var)
  )
  coxph(fml, data = data, x = TRUE, y = TRUE)
}

safe_value <- function(x, name) {
  if (!is.null(x[[name]])) {
    x[[name]]
  } else {
    NA_real_
  }
}

# C-index comparison.
run_cindex <- function(score_var, model_name) {
  d <- make_complete_data(score_var)
  model_score2 <- fit_score2_model(d)
  model_plus <- fit_score2_plus_model(d, score_var)

  cindex_score2 <- concordance.index(
    predict(model_score2, type = "lp"),
    surv.time = d$MACE1_days,
    surv.event = d$MACE1_event,
    method = "noether"
  )

  cindex_plus <- concordance.index(
    predict(model_plus, type = "lp"),
    surv.time = d$MACE1_days,
    surv.event = d$MACE1_event,
    method = "noether"
  )

  cindex_test <- cindex.comp(cindex_plus, cindex_score2)

  tibble(
    model_name = model_name,
    score_var = score_var,
    n = nrow(d),
    events = sum(d$MACE1_event == 1, na.rm = TRUE),
    cindex_score2 = safe_value(cindex_score2, "c.index"),
    cindex_plus = safe_value(cindex_plus, "c.index"),
    delta_cindex = safe_value(cindex_plus, "c.index") - safe_value(cindex_score2, "c.index"),
    p_value = safe_value(cindex_test, "p.value")
  )
}

cindex_results <- purrr::map2_dfr(
  score_info$score_var,
  score_info$model_name,
  run_cindex
)

write.csv(cindex_results, "results/C_index_results.csv", row.names = FALSE)

# NRI and IDI at 10 years.
run_reclassification <- function(score_var, model_name) {
  d <- make_complete_data(score_var)

  nri_result <- nricens(
    time = d$MACE1_days,
    event = d$MACE1_event,
    z.std = d[, "SCORE2_linear_predictor", drop = FALSE],
    z.new = d[, c("SCORE2_linear_predictor", score_var), drop = FALSE],
    t0 = horizon_days,
    cut = risk_cuts,
    niter = 1000,
    updown = "category"
  )

  idi_result <- IDI.INF(
    indata = as.matrix(d[, c("MACE1_days", "MACE1_event")]),
    covs0 = as.matrix(d[, "SCORE2_linear_predictor", drop = FALSE]),
    covs1 = as.matrix(d[, c("SCORE2_linear_predictor", score_var), drop = FALSE]),
    t0 = horizon_days,
    npert = 300,
    alpha = 0.05
  )

  nri_text <- capture.output(print(nri_result))
  idi_text <- capture.output(IDI.INF.OUT(idi_result))

  writeLines(
    c(
      paste0("Model: ", model_name),
      paste0("Score: ", score_var),
      "",
      "NRI result:",
      nri_text,
      "",
      "IDI result:",
      idi_text
    ),
    con = file.path("results", paste0("NRI_IDI_", model_name, ".txt"))
  )

  tibble(
    model_name = model_name,
    score_var = score_var,
    n = nrow(d),
    events_10y = sum(d$MACE1_days <= horizon_days & d$MACE1_event == 1, na.rm = TRUE),
    nri_output_file = file.path("results", paste0("NRI_IDI_", model_name, ".txt"))
  )
}

reclassification_manifest <- purrr::map2_dfr(
  score_info$score_var,
  score_info$model_name,
  run_reclassification
)

write.csv(reclassification_manifest, "results/NRI_IDI_manifest.csv", row.names = FALSE)

# 10-year calibration plots.
plot_calibration <- function(risk_pr, t_var, out_var, subtitle) {
  scalib_object <- scalib(
    pred_risk = risk_pr,
    pred_horizon = 3650,
    event_time = t_var,
    event_status = out_var
  )

  scalib_slope <- scalib_hare(scalib_object = scalib_object)
  slope_data <- scalib_slope %>%
    getElement("data_outputs") %>%
    select(._id_., hare_data_plot) %>%
    unnest(hare_data_plot)

  calibration_curve <- ggplot(slope_data, aes(x = predicted, y = observed)) +
    geom_line(color = "blue") +
    geom_abline(color = "grey", linetype = 2, intercept = 0, slope = 1) +
    scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 1, by = 0.1)) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 1, by = 0.1))

  grouped_calibration <- scalib_gnd(
    scalib_object = scalib_object,
    group_count_init = 10
  ) %>%
    getElement("data_outputs") %>%
    unnest(gnd_data)

  calibration_curve +
    geom_pointrange(
      data = grouped_calibration,
      mapping = aes(
        x = percent_expected,
        y = percent_observed,
        ymin = percent_observed - 1.96 * sqrt(variance),
        ymax = percent_observed + 1.96 * sqrt(variance)
      ),
      shape = 21,
      fill = "grey",
      color = "black"
    ) +
    labs(subtitle = subtitle, x = "Predicted 10-year risk", y = "Observed 10-year risk") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      plot.subtitle = element_text(size = 14, face = "bold")
    )
}

calibration_data <- make_complete_data()
calibration_models <- list(
  score2 = list(
    label = "SCORE2",
    model = fit_score2_model(calibration_data)
  )
)

for (i in seq_len(nrow(score_info))) {
  score_var <- score_info$score_var[i]
  model_name <- score_info$model_name[i]
  d <- make_complete_data(score_var)
  calibration_models[[model_name]] <- list(
    label = score_info$label[i],
    model = fit_score2_plus_model(d, score_var),
    data = d
  )
}

calibration_plots <- lapply(names(calibration_models), function(model_name) {
  model_item <- calibration_models[[model_name]]
  newdata <- if (!is.null(model_item$data)) model_item$data else calibration_data
  predicted_risk <- predictRisk(model_item$model, newdata = newdata, times = 3650)
  plot_calibration(
    risk_pr = predicted_risk,
    t_var = newdata$MACE1_days,
    out_var = newdata$MACE1_event,
    subtitle = model_item$label
  )
})

combined_calibration_plot <- wrap_plots(calibration_plots, nrow = 1)

ggsave(
  "figures/Calibration_SCORE2_plus_AtheroBurden.png",
  combined_calibration_plot,
  width = 24,
  height = 5,
  dpi = 300
)

print(combined_calibration_plot)
