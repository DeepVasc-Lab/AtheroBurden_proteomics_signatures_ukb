# Input note:
# This script starts from an anonymous UKB prospective analysis-ready file:
# data/ukb_prospective_analysis.csv
#
# Main required columns:
# - Artery_enriched_248, Atherosclerosis_680, MR_derived_402,
#   Whole_proteome_2920: AtheroBurden score columns from the four CatBoost panels.
# - MACE_event_new, MACE_days_new, Stroke_event_new, Stroke_days_new,
#   AMI_event_new, AMI_days_new, CVdeath_event_new, CVdeath_days_new:
#   prospective endpoint indicators and follow-up time in days.
# - Age, Sex, ethnicity, Cholesterol, HDL, SBP, BMI, LDL, Triglycerides,
#   eGFRcr_cys, HbA1c, diabetes_status, hypert_status:
#   clinical covariates used in adjusted Cox models.
# - Current_smoke, Unknown_smoke, Current_drink, Unknown_drink, fh_cvd:
#   smoking, drinking, and family-history covariates.
#
# The 10-year ROC panel uses data/ukb_timeROC_curve_10y.csv, a curve-point
# table generated from the SCORE2 and SCORE2-plus-AtheroBurden timeROC models.
# Required columns for that file: time, model, FPR, TPR.

library(data.table)
library(dplyr)
library(purrr)
library(survival)
library(survminer)
library(ggplot2)
library(ggtext)
library(patchwork)
library(tibble)
library(export)

dat <- fread("data/ukb_prospective_analysis.csv", data.table = FALSE)
data_for_timeROC <- dat

# Calculate res2: Cox models for each endpoint, adjustment model, and signature.
score_vars <- c("Artery_enriched_248", "MR_derived_402", "Atherosclerosis_680", "Whole_proteome_2920")

score_name_map <- c(
  Artery_enriched_248 = "AtheroBurden Arterial Signature",
  MR_derived_402 = "AtheroBurden Genetic Signature",
  Atherosclerosis_680 = "AtheroBurden Mechanistic Signature",
  Whole_proteome_2920 = "AtheroBurden WholeProteome Signature"
)

score_order <- unname(score_name_map[score_vars])

endpoints <- tribble(
  ~Endpoint,  ~event_col,           ~time_col,
  "MACE",     "MACE_event_new",     "MACE_days_new",
  "Stroke",   "Stroke_event_new",   "Stroke_days_new",
  "AMI",      "AMI_event_new",      "AMI_days_new",
  "CV death", "CVdeath_event_new",  "CVdeath_days_new"
)

dat <- dat %>%
  mutate(
    Sex = factor(Sex),
    ethnicity6 = as.character(ethnicity),
    ethnicity6 = ifelse(is.na(ethnicity6) | ethnicity6 == "", "Unknown", ethnicity6),
    ethnicity6 = factor(ethnicity6, levels = c("White", "Mixed", "Asian", "Black", "Other", "Unknown")),
    fh_cvd = as.character(fh_cvd),
    fh_cvd = ifelse(is.na(fh_cvd) | fh_cvd == "", "Unknown", fh_cvd),
    fh_cvd = factor(fh_cvd, levels = c("No", "Yes", "Unknown")),
    Current_smoke = ifelse(is.na(Current_smoke), 0L, as.integer(Current_smoke)),
    Unknown_smoke = ifelse(is.na(Unknown_smoke), 0L, as.integer(Unknown_smoke)),
    Current_drink = ifelse(is.na(Current_drink), 0L, as.integer(Current_drink)),
    Unknown_drink = ifelse(is.na(Unknown_drink), 0L, as.integer(Unknown_drink)),
    diabetes_status = factor(diabetes_status),
    hypert_status = factor(hypert_status)
  )

model_covars <- list(
  model1 = c("Age", "Sex", "ethnicity6"),
  model2 = c("Age", "Sex", "Cholesterol", "HDL", "SBP", "Current_smoke", "Unknown_smoke"),
  model3 = c(
    "Age", "Sex", "SBP", "BMI", "Current_smoke", 
    "LDL", "Triglycerides", "eGFRcr_cys", "HbA1c",
    "diabetes_status", "hypert_status",
    "ethnicity6", "Current_drink",  "fh_cvd"
  )
)

fmt_hr_ci <- function(beta, se, digits = 3) {
  hr <- exp(beta)
  lo <- exp(beta - 1.96 * se)
  hi <- exp(beta + 1.96 * se)
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), hr, lo, hi)
}

fit_one <- function(dat, time_col, event_col, score_var, covars, model_name, endpoint_name) {
  d <- dat %>%
    mutate(
      .time = suppressWarnings(as.numeric(.data[[time_col]])),
      .event = ifelse(is.na(.data[[event_col]]), 0L, as.integer(.data[[event_col]]))
    ) %>%
    filter(!is.na(.time) & .time > 0) %>%
    filter(!is.na(.event)) %>%
    filter(if_all(all_of(c(score_var, covars)), ~ !is.na(.x)))

  rhs <- paste(c(score_var, covars), collapse = " + ")
  fml <- as.formula(paste0("Surv(.time, .event) ~ ", rhs))

  tryCatch({
    fit <- coxph(fml, data = d, ties = "efron")
    s <- summary(fit)

    beta <- s$coefficients[score_var, "coef"]
    se <- s$coefficients[score_var, "se(coef)"]
    p <- s$coefficients[score_var, "Pr(>|z|)"]

    tibble(
      `AtheroBurden Signatures` = unname(score_name_map[[score_var]]),
      estimate = beta,
      std.error = se,
      `p value` = p,
      model = model_name,
      Endpoint = endpoint_name,
      `HR (95% CI)` = fmt_hr_ci(beta, se)
    )
  }, error = function(e) {
    tibble(
      `AtheroBurden Signatures` = unname(score_name_map[[score_var]]),
      estimate = NA_real_,
      std.error = NA_real_,
      `p value` = NA_real_,
      model = model_name,
      Endpoint = endpoint_name,
      `HR (95% CI)` = NA_character_
    )
  })
}

cox_grid <- expand.grid(
  Endpoint = endpoints$Endpoint,
  model = names(model_covars),
  score = score_vars,
  stringsAsFactors = FALSE
)

res <- pmap_dfr(cox_grid, function(Endpoint, model, score) {
  ep <- endpoints %>% filter(Endpoint == !!Endpoint)
  fit_one(
    dat = dat,
    time_col = ep$time_col[[1]],
    event_col = ep$event_col[[1]],
    score_var = score,
    covars = model_covars[[model]],
    model_name = model,
    endpoint_name = Endpoint
  )
})

res2 <- res %>%
  group_by(Endpoint, model) %>%
  mutate(
    `p-adjust` = p.adjust(`p value`, method = "fdr"),
    Significance = case_when(
      `p-adjust` < 0.001 ~ "***",
      `p-adjust` < 0.01 ~ "**",
      `p-adjust` < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup() %>%
  mutate(
    `AtheroBurden Signatures` = factor(`AtheroBurden Signatures`, levels = score_order),
    model = factor(model, levels = c("model1", "model2", "model3")),
    Endpoint = factor(Endpoint, levels = c("MACE", "Stroke", "AMI", "CV death"))
  ) %>%
  arrange(Endpoint, model, `AtheroBurden Signatures`) %>%
  mutate(
    `AtheroBurden Signatures` = as.character(`AtheroBurden Signatures`),
    model = as.character(model),
    Endpoint = as.character(Endpoint)
  ) %>%
  select(
    `AtheroBurden Signatures`, estimate, std.error, `p value`,
    model, `p-adjust`, `HR (95% CI)`, Significance, Endpoint
  )

dir.create("results", showWarnings = FALSE)
write.csv(res2, "results/Figure4_cox_results.csv", row.names = FALSE)

#========================

#========================
base_theme <- theme_classic(base_size = 22) +
  theme(
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text  = element_text(size = 22, face = "bold"),
    legend.text  = element_text(size = 22),
    legend.title = element_text(size = 22, face = "bold")
  )

#========================

#========================
cox_plot_df <- res2 %>%
  mutate(
    Endpoint = case_when(
      Endpoint %in% c("CV death", "CV Death") ~ "CV Death",
      TRUE ~ Endpoint
    ),
    HR = exp(estimate),
    CI_lower = exp(estimate - 1.96 * std.error),
    CI_upper = exp(estimate + 1.96 * std.error)
  ) %>%
  filter(!is.na(HR) & !is.na(CI_lower) & !is.na(CI_upper))

#========================

#========================
x_max_fixed <- 2.408
x_breaks_fixed <- c(1.0, 1.5, 2.0, 2.4)


# max(cox_plot_df$CI_upper)


x_min_fixed <- floor(min(cox_plot_df$CI_lower, na.rm = TRUE) * 10) / 10
x_min_fixed <- max(0.6, x_min_fixed)

#========================

#========================
create_forest_plot <- function(data, title) {
  

  data$Endpoint <- factor(data$Endpoint, levels = c("CV Death", "AMI", "Stroke", "MACE"))
  

  data$y_base <- as.numeric(data$Endpoint)
  

  data$y_position <- data$y_base +
    ifelse(data$model == "model1",  0.30,
           ifelse(data$model == "model2", 0.00,
                  ifelse(data$model == "model3", -0.30, 0.00)))
  

  background_data <- data.frame(
    Endpoint = factor(c("CV Death", "AMI", "Stroke", "MACE"),
                      levels = c("CV Death", "AMI", "Stroke", "MACE")),
    ymin = c(0.5, 1.5, 2.5, 3.5),
    ymax = c(1.5, 2.5, 3.5, 4.5)
  )
  

  data$star_x <- pmin(data$CI_upper + 0.03, x_max_fixed - 0.01)
  
  p <- ggplot() +
    geom_rect(
      data = background_data,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = Endpoint),
      alpha = 0.10
    ) +

    geom_segment(
      data = data,
      aes(x = CI_lower, xend = CI_upper, y = y_position, yend = y_position),
      linewidth = 1, color = "black"
    ) +
    geom_point(
      data = data,
      aes(x = HR, y = y_position, color = model),
      size = 6
    ) +

    geom_text(
      data = data,
      aes(x = star_x, y = y_position, label = Significance, color = model),
      size = 6, hjust = 0
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    scale_fill_manual(values = c("MACE" = "#bdbdbd", "Stroke" = "#e6eaf2",
                                 "AMI" = "#bdbdbd", "CV Death" = "#e6eaf2"),
                      guide = "none") +
    scale_color_manual(
      values = c("model1" = "#e66101", "model2" = "#fdb863", "model3" = "#5e3c99"),
      labels = c(
        "model1" = "Adjusted for Age, Sex and Ethnicity",
        "model2" = "Adjusted for SCORE2 Variables",
        "model3" = "Adjusted for VRFs"
      )
    ) +

    scale_y_continuous(
      breaks = 1:4,
      labels = c("CV Death", "AMI", "Ischemic\nStroke", "MACE"),
      limits = c(0.5, 4.5),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      breaks = x_breaks_fixed,
      labels = function(x) format(x, nsmall = 1)
    ) +
    coord_cartesian(xlim = c(x_min_fixed, x_max_fixed)) +
    labs(title = title, x = "Hazard Ratio (95% CI)", y = "") +
    base_theme +
    theme(
      strip.text = element_text(face = "bold", size = 22, hjust = 0.5),

      legend.position = "none",
      plot.title = element_text(size = 22, face = "bold", hjust = 0)
    )
  
  p
}

#========================

#========================
cox_results_combined <- cox_plot_df

data_Artery <- cox_results_combined %>% filter(`AtheroBurden Signatures` == "AtheroBurden Arterial Signature")
data_Atherosclerosis <- cox_results_combined %>% filter(`AtheroBurden Signatures` == "AtheroBurden Mechanistic Signature")
data_MR <- cox_results_combined %>% filter(`AtheroBurden Signatures` == "AtheroBurden Genetic Signature")
data_Whole <- cox_results_combined %>% filter(`AtheroBurden Signatures` == "AtheroBurden WholeProteome Signature")

plot1 <- create_forest_plot(data_Artery, "AtheroBurden Arterial Signature")
plot2 <- create_forest_plot(data_MR, "AtheroBurden Genetic Signature")
plot3 <- create_forest_plot(data_Atherosclerosis, "AtheroBurden Mechanistic Signature")
plot4 <- create_forest_plot(data_Whole, "AtheroBurden WholeProteome Signature")

plot1
plot2
plot3
plot4


get_legend <- function(a_plot) {
  tmp <- ggplotGrob(a_plot + theme(legend.position = "bottom"))
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}
legend_plot <- get_legend(plot1)
#======================

#======================
df <- data_for_timeROC

#======================

#======================
base_theme <- theme_classic(base_size = 22) +
  theme(
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22, face = "bold")
  )

#======================

#======================
df <- df %>%
  mutate(
    MACE_days_new_y = MACE_days_new / 365.25,
    MACE_event_new  = as.integer(MACE_event_new)
  )

#======================


#======================
df <- df %>%
  mutate(
    Artery_enriched_248_quartile = ifelse(is.na(Artery_enriched_248), NA_integer_,
                                          dplyr::ntile(Artery_enriched_248, 4)),
    Atherosclerosis_680_quartile = ifelse(is.na(Atherosclerosis_680), NA_integer_,
                                          dplyr::ntile(Atherosclerosis_680, 4)),
    MR_derived_402_quartile      = ifelse(is.na(MR_derived_402), NA_integer_,
                                          dplyr::ntile(MR_derived_402, 4)),
    Whole_proteome_2920_quartile = ifelse(is.na(Whole_proteome_2920), NA_integer_,
                                          dplyr::ntile(Whole_proteome_2920, 4))
  ) %>%
  mutate(
    Artery_enriched_248_quartile = factor(Artery_enriched_248_quartile, levels = 1:4, labels = c("Q1","Q2","Q3","Q4")),
    Atherosclerosis_680_quartile = factor(Atherosclerosis_680_quartile, levels = 1:4, labels = c("Q1","Q2","Q3","Q4")),
    MR_derived_402_quartile      = factor(MR_derived_402_quartile, levels = 1:4, labels = c("Q1","Q2","Q3","Q4")),
    Whole_proteome_2920_quartile = factor(Whole_proteome_2920_quartile, levels = 1:4, labels = c("Q1","Q2","Q3","Q4"))
  )


# table(df$Artery_enriched_248_quartile, useNA="ifany")

# =====================================================================

# =====================================================================

# ----- p1: Artery_enriched_248_quartile -----
log_rank_test <- survdiff(Surv(MACE_days_new_y, MACE_event_new) ~ Artery_enriched_248_quartile,
                          data = df, na.action = na.exclude)
p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

fit1 <- survfit(Surv(MACE_days_new_y, MACE_event_new) ~ Artery_enriched_248_quartile,
                data = df, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain", na.action = na.exclude)

cox_model1 <- coxph(Surv(MACE_days_new_y, MACE_event_new) ~ Artery_enriched_248_quartile +
                      Age + Sex + Current_smoke + HDL + Cholesterol + SBP,
                    data = df)
cox_summary1 <- summary(cox_model1)
coef_indices_1 <- grep("Artery_enriched_248_quartile", rownames(cox_summary1$coefficients))
hr_vals_1 <- exp(cox_summary1$coefficients[coef_indices_1, "coef"])
ci_vals_1 <- cox_summary1$conf.int[coef_indices_1, c("lower .95", "upper .95")]

legend.labs1 <- c(
  "Q1, Reference",
  paste0("Q2, HR = ", round(hr_vals_1[1], 2), " (95% CI, ", round(ci_vals_1[1,1], 2), "-", round(ci_vals_1[1,2], 2), ")"),
  paste0("Q3, HR = ", round(hr_vals_1[2], 2), " (95% CI, ", round(ci_vals_1[2,1], 2), "-", round(ci_vals_1[2,2], 2), ")"),
  paste0("Q4, HR = ", round(hr_vals_1[3], 2), " (95% CI, ", round(ci_vals_1[3,1], 2), "-", round(ci_vals_1[3,2], 2), ")")
)

p1_surv <- ggsurvplot(
  fit1,
  fun = "event",
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  ncensor.plot = FALSE,
  censor.shape = NA,
  linetype = "solid",
  ggtheme = base_theme,
  palette = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99"),
  xlab = "Time (years)",
  ylab = "Cumulative incidence rate of MACE (%)",
  legend.labs = legend.labs1,
  risk.table.font = 7.5,
  legend.title = "ArtheroBurden Arterial Signature",
  xlim = c(0, 15),
  break.time.by = 3
)

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.0001, " < 0.0001",
                       paste0(" = ", round(p.val, 3))))

p1_surv$plot <- p1_surv$plot +
  scale_y_continuous(labels = function(x) x * 100) +
  annotate("text", x = 0, y = 0.1, label = p.lab, hjust = 0, fontface = "bold", size = 8) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.direction = "vertical",
    plot.title = element_text(size = 24, face = "bold")
  )

p1_surv$table <- p1_surv$table +
  scale_y_discrete(labels = function(x) {
    groups <- c("Q1", "Q2", "Q3", "Q4")
    colors <- c("#e66101", "#fdb863", "#b2abd2", "#5e3c99")
    rev(paste0("<span style='color:", colors, "'>", groups, "</span>"))
  }) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown(size = 24),
    plot.margin = margin(0, 0, 0, 0)
  )

# ----- p2: Atherosclerosis_680_quartile -----
log_rank_test <- survdiff(Surv(MACE_days_new_y, MACE_event_new) ~ Atherosclerosis_680_quartile,
                          data = df, na.action = na.exclude)
p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

fit2 <- survfit(Surv(MACE_days_new_y, MACE_event_new) ~ Atherosclerosis_680_quartile,
                data = df, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain", na.action = na.exclude)

cox_model2 <- coxph(Surv(MACE_days_new_y, MACE_event_new) ~ Atherosclerosis_680_quartile +
                      Age + Sex + Current_smoke + HDL + Cholesterol + SBP,
                    data = df)
cox_summary2 <- summary(cox_model2)
coef_indices_2 <- grep("Atherosclerosis_680_quartile", rownames(cox_summary2$coefficients))
hr_vals_2 <- exp(cox_summary2$coefficients[coef_indices_2, "coef"])
ci_vals_2 <- cox_summary2$conf.int[coef_indices_2, c("lower .95", "upper .95")]

legend.labs2 <- c(
  "Q1, Reference",
  paste0("Q2, HR = ", round(hr_vals_2[1], 2), " (95% CI, ", round(ci_vals_2[1,1], 2), "-", round(ci_vals_2[1,2], 2), ")"),
  paste0("Q3, HR = ", round(hr_vals_2[2], 2), " (95% CI, ", round(ci_vals_2[2,1], 2), "-", round(ci_vals_2[2,2], 2), ")"),
  paste0("Q4, HR = ", round(hr_vals_2[3], 2), " (95% CI, ", round(ci_vals_2[3,1], 2), "-", round(ci_vals_2[3,2], 2), ")")
)

p2_surv <- ggsurvplot(
  fit2,
  fun = "event",
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  ncensor.plot = FALSE,
  censor.shape = NA,
  linetype = "solid",
  ggtheme = base_theme,
  palette = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99"),
  xlab = "Time (years)",
  ylab = "Cumulative incidence rate of MACE (%)",
  legend.labs = legend.labs2,
  risk.table.font = 7.5,
  legend.title = "ArtheroBurden Mechanistic Signature",
  xlim = c(0, 15),
  break.time.by = 3
)

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.0001, " < 0.0001",
                       paste0(" = ", round(p.val, 3))))

p2_surv$plot <- p2_surv$plot +
  scale_y_continuous(labels = function(x) x * 100) +
  annotate("text", x = 0, y = 0.1, label = p.lab, hjust = 0, fontface = "bold", size = 8) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.direction = "vertical",
    plot.title = element_text(size = 24, face = "bold")
  )

p2_surv$table <- p2_surv$table +
  scale_y_discrete(labels = function(x) {
    groups <- c("Q1", "Q2", "Q3", "Q4")
    colors <- c("#e66101", "#fdb863", "#b2abd2", "#5e3c99")
    rev(paste0("<span style='color:", colors, "'>", groups, "</span>"))
  }) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown(size = 24),
    plot.margin = margin(0, 0, 0, 0)
  )

# ----- p3: MR_derived_402_quartile -----
log_rank_test <- survdiff(Surv(MACE_days_new_y, MACE_event_new) ~ MR_derived_402_quartile,
                          data = df, na.action = na.exclude)
p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

fit3 <- survfit(Surv(MACE_days_new_y, MACE_event_new) ~ MR_derived_402_quartile,
                data = df, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain", na.action = na.exclude)

cox_model3 <- coxph(Surv(MACE_days_new_y, MACE_event_new) ~ MR_derived_402_quartile +
                      Age + Sex + Current_smoke + HDL + Cholesterol + SBP,
                    data = df)
cox_summary3 <- summary(cox_model3)
coef_indices_3 <- grep("MR_derived_402_quartile", rownames(cox_summary3$coefficients))
hr_vals_3 <- exp(cox_summary3$coefficients[coef_indices_3, "coef"])
ci_vals_3 <- cox_summary3$conf.int[coef_indices_3, c("lower .95", "upper .95")]

legend.labs3 <- c(
  "Q1, Reference",
  paste0("Q2, HR = ", round(hr_vals_3[1], 2), " (95% CI, ", round(ci_vals_3[1,1], 2), "-", round(ci_vals_3[1,2], 2), ")"),
  paste0("Q3, HR = ", round(hr_vals_3[2], 2), " (95% CI, ", round(ci_vals_3[2,1], 2), "-", round(ci_vals_3[2,2], 2), ")"),
  paste0("Q4, HR = ", round(hr_vals_3[3], 2), " (95% CI, ", round(ci_vals_3[3,1], 2), "-", round(ci_vals_3[3,2], 2), ")")
)

p3_surv <- ggsurvplot(
  fit3,
  fun = "event",
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  ncensor.plot = FALSE,
  censor.shape = NA,
  linetype = "solid",
  ggtheme = base_theme,
  palette = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99"),
  xlab = "Time (years)",
  ylab = "Cumulative incidence rate of MACE (%)",
  legend.labs = legend.labs3,
  risk.table.font = 7.5,
  legend.title = "ArtheroBurden Genetic Signature",
  xlim = c(0, 15),
  break.time.by = 3
)

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.0001, " < 0.0001",
                       paste0(" = ", round(p.val, 3))))

p3_surv$plot <- p3_surv$plot +
  scale_y_continuous(labels = function(x) x * 100) +
  annotate("text", x = 0, y = 0.1, label = p.lab, hjust = 0, fontface = "bold", size = 8) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.direction = "vertical",
    plot.title = element_text(size = 24, face = "bold")
  )

p3_surv$table <- p3_surv$table +
  scale_y_discrete(labels = function(x) {
    groups <- c("Q1", "Q2", "Q3", "Q4")
    colors <- c("#e66101", "#fdb863", "#b2abd2", "#5e3c99")
    rev(paste0("<span style='color:", colors, "'>", groups, "</span>"))
  }) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown(size = 24),
    plot.margin = margin(0, 0, 0, 0)
  )

# ----- p4: Whole_proteome_2920_quartile -----
log_rank_test <- survdiff(Surv(MACE_days_new_y, MACE_event_new) ~ Whole_proteome_2920_quartile,
                          data = df, na.action = na.exclude)
p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

fit4 <- survfit(Surv(MACE_days_new_y, MACE_event_new) ~ Whole_proteome_2920_quartile,
                data = df, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain", na.action = na.exclude)

cox_model4 <- coxph(Surv(MACE_days_new_y, MACE_event_new) ~ Whole_proteome_2920_quartile +
                      Age + Sex + Current_smoke + HDL + Cholesterol + SBP,
                    data = df)
cox_summary4 <- summary(cox_model4)
coef_indices_4 <- grep("Whole_proteome_2920_quartile", rownames(cox_summary4$coefficients))
hr_vals_4 <- exp(cox_summary4$coefficients[coef_indices_4, "coef"])
ci_vals_4 <- cox_summary4$conf.int[coef_indices_4, c("lower .95", "upper .95")]

legend.labs4 <- c(
  "Q1, Reference",
  paste0("Q2, HR = ", round(hr_vals_4[1], 2), " (95% CI, ", round(ci_vals_4[1,1], 2), "-", round(ci_vals_4[1,2], 2), ")"),
  paste0("Q3, HR = ", round(hr_vals_4[2], 2), " (95% CI, ", round(ci_vals_4[2,1], 2), "-", round(ci_vals_4[2,2], 2), ")"),
  paste0("Q4, HR = ", round(hr_vals_4[3], 2), " (95% CI, ", round(ci_vals_4[3,1], 2), "-", round(ci_vals_4[3,2], 2), ")")
)

p4_surv <- ggsurvplot(
  fit4,
  fun = "event",
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  ncensor.plot = FALSE,
  censor.shape = NA,
  linetype = "solid",
  ggtheme = base_theme,
  palette = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99"),
  xlab = "Time (years)",
  ylab = "Cumulative incidence rate of MACE (%)",
  legend.labs = legend.labs4,
  risk.table.font = 7.5,
  legend.title = "ArtheroBurden WholeProteome Signature",
  xlim = c(0, 15),
  break.time.by = 3
)

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.0001, " < 0.0001",
                       paste0(" = ", round(p.val, 3))))

p4_surv$plot <- p4_surv$plot +
  scale_y_continuous(labels = function(x) x * 100) +
  annotate("text", x = 0, y = 0.1, label = p.lab, hjust = 0, fontface = "bold", size = 8) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.direction = "vertical",
    plot.title = element_text(size = 24, face = "bold")
  )

p4_surv$table <- p4_surv$table +
  scale_y_discrete(labels = function(x) {
    groups <- c("Q1", "Q2", "Q3", "Q4")
    colors <- c("#e66101", "#fdb863", "#b2abd2", "#5e3c99")
    rev(paste0("<span style='color:", colors, "'>", groups, "</span>"))
  }) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown(size = 24),
    plot.margin = margin(0, 0, 0, 0)
  )



col1 <- (plot1 / plot_spacer() / p1_surv$plot / p1_surv$table) +
  plot_layout(heights = c(0.8, 0.3, 1, 0.3))
col2 <- (plot2 / plot_spacer() / p2_surv$plot / p2_surv$table) +
  plot_layout(heights = c(0.8, 0.3, 1, 0.3))
col3 <- (plot3 / plot_spacer() / p3_surv$plot / p3_surv$table) +
  plot_layout(heights = c(0.8, 0.3, 1, 0.3))
col4 <- (plot4 / plot_spacer() / p4_surv$plot / p4_surv$table) +
  plot_layout(heights = c(0.8, 0.3, 1, 0.3))

final_surv_plot <- ((col4 | col3 | col2 | col1) / legend_plot) +
  plot_layout(heights = c(1, 0.1))

final_surv_plot


get_km_cuminc_percent <- function(fit, t = 16) {
  s <- summary(fit, times = t, extend = TRUE)
  out <- data.frame(
    strata  = s$strata,
    time    = s$time,
    n_risk  = s$n.risk,
    cuminc  = 1 - s$surv
  )
  out$cuminc_percent <- 100 * out$cuminc
  out
}

km1_16y <- get_km_cuminc_percent(fit1, t = 16)
km1_16y

extract_Q1_Q4 <- function(km_df) {
  q1 <- km_df[km_df$strata %in% grep("Q1", km_df$strata, value = TRUE), ]
  q4 <- km_df[km_df$strata %in% grep("Q4", km_df$strata, value = TRUE), ]
  data.frame(
    Q1_percent = q1$cuminc_percent,
    Q4_percent = q4$cuminc_percent,
    Q1_n_risk  = q1$n_risk,
    Q4_n_risk  = q4$n_risk
  )
}

extract_Q1_Q4(km1_16y)

fits <- list(
  Arterial       = fit1,
  Mechanistic    = fit2,
  Genetic        = fit3,
  WholeProteome  = fit4
)

km_Q1Q4_16y <- do.call(rbind, lapply(names(fits), function(nm) {
  km <- get_km_cuminc_percent(fits[[nm]], t = 16)
  q  <- extract_Q1_Q4(km)
  cbind(Signature = nm, q)
}))

km_Q1Q4_16y
range(km_Q1Q4_16y$Q4_percent, na.rm = TRUE)
range(km_Q1Q4_16y$Q1_percent, na.rm = TRUE)


timeROCdata <- fread("data/ukb_timeROC_curve_10y.csv", data.table = FALSE)
df_filtered_3652 <- subset(timeROCdata, time == 3652.500)
df_filtered_3652$model <- factor(
  df_filtered_3652$model,
  levels = c("Score2", 
             "B1", 
             "B2", 
             "B3", 
             "B4")
)

df_filtered_248 <- df_filtered_3652 %>% filter(model %in% c("Score2", "B1"))
df_filtered_680 <- df_filtered_3652 %>% filter(model %in% c("Score2", "B2"))
df_filtered_402 <- df_filtered_3652 %>% filter(model %in% c("Score2", "B3"))
df_filtered_2920 <- df_filtered_3652 %>% filter(model %in% c("Score2", "B4"))

df_filtered_248$model  <- factor(df_filtered_248$model, levels = c("Score2", "B1"))
df_filtered_680$model  <- factor(df_filtered_680$model, levels = c("Score2", "B2"))
df_filtered_402$model  <- factor(df_filtered_402$model, levels = c("Score2", "B3"))
df_filtered_2920$model  <- factor(df_filtered_2920$model, levels = c("Score2", "B4"))

p_248 <- ggplot(df_filtered_248, aes(x = FPR, y = TPR, color = model)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, color = "grey10", linetype = 2) +
  scale_color_manual(name = NULL,values = c("#e66101", "#5e3c99"),
                     labels = c("SCORE2: 0.705 (0.694 - 0.716)", 
                                "plus Arterial Signature: 0.744 (0.733 - 0.755)")) +
  base_theme +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = "ArtheroBurden Arterial Signature") +
  theme(
    legend.position = c(0.55, 0.1),
    legend.box.background = element_rect(color = "grey80", fill = "white", size = 0.5),
    legend.key = element_rect(fill = "white", color = NA),
    legend.spacing.x = unit(1, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 28, face = "bold")
  )
print(p_248)

p_680 <- ggplot(df_filtered_680, aes(x = FPR, y = TPR, color = model)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, color = "grey10", linetype = 2) +
  scale_color_manual(name = NULL,values = c("#e66101", "#5e3c99"),
                     labels = c("SCORE2: 0.705 (0.694 - 0.716)", 
                                "plus Mechanistic Signature: 0.742 (0.731 - 0.753)")) +
  base_theme +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = "ArtheroBurden Mechanistic Signature") +
  theme(
    legend.position = c(0.55, 0.1),
    legend.box.background = element_rect(color = "grey80", fill = "white", size = 0.5),
    legend.key = element_rect(fill = "white", color = NA),
    legend.spacing.x = unit(1, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 28, face = "bold")
  )
print(p_680)

p_402 <- ggplot(df_filtered_402, aes(x = FPR, y = TPR, color = model)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, color = "grey10", linetype = 2) +
  scale_color_manual(name = NULL,values = c("#e66101", "#5e3c99"),
                     labels = c("SCORE2: 0.705 (0.694 - 0.716)", 
                                "plus Genetic Signature: 0.748 (0.738 - 0.759)")) +
  base_theme +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = "ArtheroBurden Genetic Signature") +
  theme(
    legend.position = c(0.55, 0.1),
    legend.box.background = element_rect(color = "grey80", fill = "white", size = 0.5),
    legend.key = element_rect(fill = "white", color = NA),
    legend.spacing.x = unit(1, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 28, face = "bold")
  )
print(p_402)

p_2920 <- ggplot(df_filtered_2920, aes(x = FPR, y = TPR, color = model)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, color = "grey10", linetype = 2) +
  scale_color_manual(name = NULL,values = c("#e66101", "#5e3c99"),
                     labels = c("SCORE2: 0.705 (0.694 - 0.716)", 
                                "plus WholeProteome Signature: 0.750 (0.740 - 0.761)")) +
  base_theme +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", title = "ArtheroBurden WholeProteome Signature") +
  theme(
    legend.position = c(0.55, 0.1),
    legend.box.background = element_rect(color = "grey80", fill = "white", size = 0.5),
    legend.key = element_rect(fill = "white", color = NA),
    legend.spacing.x = unit(1, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 28, face = "bold")
  )
print(p_2920)


final_ROC_plot <- (p_2920 | p_402 | p_680 | p_248)
final_ROC_plot

graph2ppt(final_surv_plot,file = "figures/Figure4ab.pptx",height = 100, width = 180)
graph2ppt(final_ROC_plot,file = "figures/Figure4c.pptx",height = 50, width = 180)

final_surv_plot
final_ROC_plot

# 42*28
# 42*10.5
