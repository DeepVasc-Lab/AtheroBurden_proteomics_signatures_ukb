
# Input note:
# This script starts from anonymous KORA external-validation analysis files.
#
# Score generation is shown in 06_KORA_score_generation.py. The generated
# KORA AtheroBurden scores should be merged with anonymous KORA clinical
# follow-up files before running this Figure 6 script.
#
# Main required input files:
# - data/kora_restricted_score_correlation.csv:
#   paired full-panel and restricted-panel KORA scores for the correlation row.
#   Required columns: Artery_enriched_248, Atherosclerosis_680,
#   MR_derived_402, Whole_proteome_2920, pred_artery, pred_arthero,
#   pred_mr, pred_whole.
# - data/kora_s4_validation_analysis.csv:
#   KORA S4 analysis-ready file with restricted AtheroBurden scores, MI/stroke
#   follow-up, and S4 covariates.
# - data/kora_age1_validation_analysis.csv:
#   KORA Age1 analysis-ready file with restricted AtheroBurden scores,
#   MI/stroke follow-up, and Age1 covariates.
#
# Required outcome columns for S4 and Age1 files: mi_apo_time and inz_mi_apo.
# S4 score columns: lo_predict_248_scale, lo_predict_680_scale,
# lo_predict_402_scale, lo_predict_2920_scale.
# Age1 score columns: wo_predict_248_scale, wo_predict_680_scale,
# wo_predict_402_scale, wo_predict_2920_scale.

library(data.table)
library(dplyr)
library(purrr)
library(tibble)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(forcats)
library(grid)
library(gridExtra)
library(ggtext)
library(patchwork)
library(RColorBrewer)
library(export)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

font_size <- 16
border_width <- 0.8

analysis_df <- fread("data/kora_restricted_score_correlation.csv", data.table = FALSE)
S4_cleaned2 <- fread("data/kora_s4_validation_analysis.csv", data.table = FALSE)
Age1_cleaned2 <- fread("data/kora_age1_validation_analysis.csv", data.table = FALSE)

S4_cleaned2$current_smoking <- ifelse(S4_cleaned2$LTCIGREG_SF %in% c(1, 2), 1, 0)
Age1_cleaned2$current_smoking <- ifelse(Age1_cleaned2$K3TCIGREG_SF %in% c(1, 2), 1, 0)

score_terms <- c(
  arterial = "Restricted Artery-Enriched ProbRS",
  mechanistic = "Restricted Atherosclerosis-Related ProbRS",
  genetic = "Restricted MR-Derived ProbRS",
  whole = "Restricted Whole Proteome ProbRS"
)

fit_kora_cox <- function(dat, score_vars, covar_sets, cohort_name) {
  imap_dfr(covar_sets, function(covars, model_name) {
    imap_dfr(score_vars, function(score_var, score_key) {
      needed <- c("mi_apo_time", "inz_mi_apo", score_var, covars)
      d <- dat %>%
        mutate(
          .time = suppressWarnings(as.numeric(mi_apo_time)),
          .event = as.integer(inz_mi_apo)
        ) %>%
        filter(!is.na(.time) & .time > 0, !is.na(.event)) %>%
        filter(if_all(all_of(needed), ~ !is.na(.x)))

      rhs <- paste(c(score_var, covars), collapse = " + ")
      fml <- as.formula(paste0("Surv(.time, .event) ~ ", rhs))

      tryCatch({
        fit <- coxph(fml, data = d, ties = "efron")
        sm <- summary(fit)
        beta <- sm$coefficients[score_var, "coef"]
        se <- sm$coefficients[score_var, "se(coef)"]
        p_value <- sm$coefficients[score_var, "Pr(>|z|)"]

        tibble(
          term = unname(score_terms[[score_key]]),
          estimate = beta,
          std.error = se,
          p.value = p_value,
          model = model_name,
          cohort = cohort_name,
          HR = exp(beta),
          CI_lower = exp(beta - 1.96 * se),
          CI_upper = exp(beta + 1.96 * se)
        )
      }, error = function(e) {
        tibble(
          term = unname(score_terms[[score_key]]),
          estimate = NA_real_,
          std.error = NA_real_,
          p.value = NA_real_,
          model = model_name,
          cohort = cohort_name,
          HR = NA_real_,
          CI_lower = NA_real_,
          CI_upper = NA_real_
        )
      })
    })
  })
}

s4_score_vars <- c(
  arterial = "lo_predict_248_scale",
  mechanistic = "lo_predict_680_scale",
  genetic = "lo_predict_402_scale",
  whole = "lo_predict_2920_scale"
)

age1_score_vars <- c(
  arterial = "wo_predict_248_scale",
  mechanistic = "wo_predict_680_scale",
  genetic = "wo_predict_402_scale",
  whole = "wo_predict_2920_scale"
)

s4_covars <- list(
  model1 = c("LTALTER", "LCSEX"),
  model2 = c("LTALTER", "LCSEX", "current_smoking", "LL_HDLN", "LL_CHOLN", "LTSYSMM"),
  model3 = c(
    "LTALTER", "LCSEX", "LTSYSMM", "LTBMI", "current_smoking",
    "LL_LDLN", "LL_TRIN", "LTGFR_CKD_CR", "LL_HBAVA", "LTDIABET", "LTWHRATC"
  )
)

age1_covars <- list(
  model1 = c("WTALTER", "WCSEX"),
  model2 = c("WTALTER", "WCSEX", "current_smoking", "wl_hdln", "wl_choln", "WTSYSMM"),
  model3 = c(
    "WTALTER", "WCSEX", "WTSYSMM", "WTBMI", "current_smoking",
    "wl_ldln", "wl_trin", "WTGFR_CKD_CR", "wl_hbava", "WTDIABET", "WTHYACT"
  )
)

combined_results <- bind_rows(
  fit_kora_cox(S4_cleaned2, s4_score_vars, s4_covars, "KORA S4"),
  fit_kora_cox(Age1_cleaned2, age1_score_vars, age1_covars, "KORA Age1")
) %>%
  group_by(cohort) %>%
  mutate(
    pajust = p.adjust(p.value, method = "fdr"),
    Significance = case_when(
      pajust < 0.001 ~ "***",
      pajust < 0.01 ~ "**",
      pajust < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

write.csv(combined_results, "results/Figure6_KORA_cox_results.csv", row.names = FALSE)

data_Artery <- combined_results %>% filter(term == score_terms[["arterial"]])
data_Atherosclerosis <- combined_results %>% filter(term == score_terms[["mechanistic"]])
data_MR <- combined_results %>% filter(term == score_terms[["genetic"]])
data_Whole <- combined_results %>% filter(term == score_terms[["whole"]])


common_theme <- theme_classic(base_size = font_size) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = border_width),
    legend.position = "right",
    text = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    axis.title = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    legend.title = element_text(size = font_size),
    plot.subtitle = element_text(size = font_size, face = "bold"),

    aspect.ratio = 0.8
  )


my_core_colors <- c("#5e3c99", "#b2abd2", "#fdb863", "#e66101")
my_hex_colors  <- colorRampPalette(my_core_colors)(10)

max_hex_count <- function(data, xvar, yvar, bins = 40) {
  df_sub <- data %>% select(all_of(c(xvar, yvar))) %>% na.omit()
  p_tmp <- ggplot(df_sub, aes_string(x = xvar, y = yvar)) +
    stat_bin2d(bins = bins)
  layer_data <- ggplot_build(p_tmp)$data[[1]]
  maxcount <- max(layer_data$count, na.rm = TRUE)
  return(maxcount)
}

max_arterial  <- max_hex_count(analysis_df, "Artery_enriched_248", "pred_artery", bins = 40)
max_mech      <- max_hex_count(analysis_df, "Atherosclerosis_680", "pred_arthero", bins = 40)
max_genetic   <- max_hex_count(analysis_df, "MR_derived_402", "pred_mr", bins = 40)
max_fullprot  <- max_hex_count(analysis_df, "Whole_proteome_2920", "pred_whole", bins = 40)
globalMaxCount <- max(max_arterial, max_mech, max_genetic, max_fullprot)

plot_hex_lm_r2 <- function(data, xvar, yvar,
                           xlab = "", ylab = "",
                           bins = 40,
                           maxcount = 1000,
                           xlim = c(-3,5), 
                           ylim = c(-3,5)) {
  ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_hex(bins = bins, color = "white", size = 0.2) +
    scale_fill_gradientn(
      colours = my_hex_colors,
      name = "Count",
      limits = c(0, maxcount)
    ) +
    geom_smooth(
      method = "lm",
      color = "#e66101",
      fill  = "#b2abd2",
      alpha = 0.3,
      size = 0.8
    ) +
    stat_cor(
      method = "pearson",
      aes(label = ifelse(
        ..p.. < 2.2e-16,
        paste0("italic(R) == ", round(..r..,2),
               "~','~~italic(p)~'<'~2.2e-16"),
        paste0("italic(R) == ", round(..r..,2),
               "~','~~italic(p) == ", format(..p.., digits = 2, scientific = TRUE))
      )),
      parse = TRUE,
      label.x.npc = "left",
      label.y.npc = "top",
      size = 5
    ) +
    labs(x = xlab, y = ylab) +
    common_theme +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = xlim) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = ylim)
}


global_xmin <- min(c(analysis_df$Artery_enriched_248, 
                     analysis_df$Atherosclerosis_680, 
                     analysis_df$MR_derived_402, 
                     analysis_df$Whole_proteome_2920), na.rm = TRUE)
global_xmax <- max(c(analysis_df$Artery_enriched_248, 
                     analysis_df$Atherosclerosis_680, 
                     analysis_df$MR_derived_402, 
                     analysis_df$Whole_proteome_2920), na.rm = TRUE)


global_ymin <- min(c(analysis_df$pred_artery, 
                     analysis_df$pred_arthero, 
                     analysis_df$pred_mr, 
                     analysis_df$pred_whole), na.rm = TRUE)
global_ymax <- max(c(analysis_df$pred_artery, 
                     analysis_df$pred_arthero, 
                     analysis_df$pred_mr, 
                     analysis_df$pred_whole), na.rm = TRUE)


p_arterial <- plot_hex_lm_r2(
  data = analysis_df,
  xvar = "Artery_enriched_248",
  yvar = "pred_artery",
  xlab = "Full AtheroBurden Arterial Signature",
  ylab = "Restricted AtheroBurden Arterial Signature",
  bins = 40,
  maxcount = globalMaxCount,
  xlim = c(global_xmin, global_xmax),
  ylim = c(global_ymin, global_ymax)
)

p_mech <- plot_hex_lm_r2(
  data = analysis_df,
  xvar = "Atherosclerosis_680",
  yvar = "pred_arthero",
  xlab = "Full AtheroBurden Mechanistic Signature",
  ylab = "Restricted AtheroBurden Mechanistic Signature",
  bins = 40,
  maxcount = globalMaxCount,
  xlim = c(global_xmin, global_xmax),
  ylim = c(global_ymin, global_ymax)
)

p_genetic <- plot_hex_lm_r2(
  data = analysis_df,
  xvar = "MR_derived_402",
  yvar = "pred_mr",
  xlab = "Full AtheroBurden Genetic Signature",
  ylab = "Restricted AtheroBurden Genetic Signature",
  bins = 40,
  maxcount = globalMaxCount,
  xlim = c(global_xmin, global_xmax),
  ylim = c(global_ymin, global_ymax)
)

p_fullprot <- plot_hex_lm_r2(
  data = analysis_df,
  xvar = "Whole_proteome_2920",
  yvar = "pred_whole",
  xlab = "Full AtheroBurden WholeProteome Signature",
  ylab = "Restricted AtheroBurden WholeProteome Signature",
  bins = 40,
  maxcount = globalMaxCount,
  xlim = c(global_xmin, global_xmax),
  ylim = c(global_ymin, global_ymax)
)


hex_plots_combined <- (p_arterial | p_mech | p_genetic | p_fullprot) +
  plot_layout(guides = "collect")




library(ggplot2)
library(dplyr)

create_forest_plot <- function(data, title) {
  # Ensure cohort is a factor and properly ordered
  data$cohort <- factor(data$cohort, levels = c("KORA Age1", "KORA S4"))
  
  # Add y_position offset
  data$y_position <- as.numeric(data$cohort) + 
    ifelse(data$model == "model1", 0.3,
           ifelse(data$model == "model2", 0,
                  ifelse(data$model == "model3", -0.3, 0)))
  
  # Define background shadow data
  background_data <- data.frame(
    cohort = c("KORA Age1", "KORA S4"),
    ymin = c(0.5, 1.5),
    ymax = c(1.5, 2.5)
  )
  


  
  # Build the plot
  p <- ggplot() +
    geom_rect(
      data = background_data,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = cohort),
      alpha = 0.1
    ) +  # Add background shadow
    geom_errorbarh(
      data = data,
      aes(xmin = CI_lower, xmax = CI_upper, y = y_position),
      height = 0.1, color = "black",
      size = 1
    ) +  # Add error bars
    geom_point(
      data = data,
      aes(x = HR, y = y_position, color = model),
      size = 8
    ) +  # Add points
    scale_y_continuous(
      breaks = c(1, 2),
      labels = c("MI/Stroke\n(KORA S4)", "MI/Stroke\n(KORA Age1)"),
      limits = c(0.4, 2.6)
    ) +  # Y-axis settings
    scale_fill_manual(
      values = c("KORA S4" = "#bdbdbd", "KORA Age1" = "#e6eaf2"),
      guide = "none"
    ) +  # Background fill
    scale_color_manual(
      values = c("model1" = "#e66101", "model2" = "#fdb863", "model3" = "#5e3c99"),
      labels = c("Adjusted for Age and Sex", "Adjusted for SCORE2 Variables", "Adjusted for VRFs")
    ) +  # Point colors
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey", size = 0.8) +
    scale_x_continuous(breaks = seq(1.0, 2.44, by = 0.5)) +
    theme_classic(base_size = 16) +
    labs(title = title, x = "Hazard Ratio (95% CI)", y = "") +
    theme(
      axis.text = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0)
    )
  

  p <- p + geom_text(
    data = data, 
    aes(x = CI_upper + 0.05, y = y_position, label = Significance, color = model), 
    size = 5,
    hjust = 0
  )
  
  return(p)
}



# global_max_ci <- max(
#   max(data_Artery$CI_upper, na.rm = TRUE),
#   max(data_Atherosclerosis$CI_upper, na.rm = TRUE),
#   max(data_MR$CI_upper, na.rm = TRUE),
#   max(data_Whole$CI_upper, na.rm = TRUE)
# ) + 0.3


plot1 <- create_forest_plot(data_Artery, "Restricted AtheroBurden Arterial")
plot2 <- create_forest_plot(data_Atherosclerosis, "Restricted AtheroBurden Mechanistic")
plot3 <- create_forest_plot(data_MR, "Restricted AtheroBurden Genetic")
plot4 <- create_forest_plot(data_Whole, "Restricted AtheroBurden WholeProteome")



library(survival)
library(survminer)
library(ggplot2)
library(forcats)
library(gridExtra)
library(ggtext)  
library(patchwork)
library(data.table)
library(dplyr)
library(RColorBrewer)

km = S4_cleaned2


base_theme <- theme_classic(base_size = font_size) +
  theme(
    axis.line = element_line(linewidth = border_width, color = "black"),
    axis.ticks = element_line(linewidth = border_width, color = "black"),
    axis.title = element_text(size = font_size, face = "bold"),
    axis.text = element_text(size = font_size, face = "bold"),
    legend.text = element_text(size = font_size),
    legend.title = element_text(size = font_size, face = "bold"),

    aspect.ratio = 1.2
  )


km$lo_predict_248_scale_q <- cut(km$lo_predict_248_scale, 
                                 breaks = quantile(km$lo_predict_248_scale, probs = seq(0, 1, by = 0.25), na.rm = TRUE), 
                                 include.lowest = TRUE, 
                                 labels = c("Q1", "Q2", "Q3", "Q4"))

km$lo_predict_680_scale_q <- cut(km$lo_predict_680_scale, 
                                 breaks = quantile(km$lo_predict_680_scale, probs = seq(0, 1, by = 0.25), na.rm = TRUE), 
                                 include.lowest = TRUE, 
                                 labels = c("Q1", "Q2", "Q3", "Q4"))

km$lo_predict_402_scale_q <- cut(km$lo_predict_402_scale, 
                                 breaks = quantile(km$lo_predict_402_scale, probs = seq(0, 1, by = 0.25), na.rm = TRUE), 
                                 include.lowest = TRUE, 
                                 labels = c("Q1", "Q2", "Q3", "Q4"))

km$lo_predict_2920_scale_q <- cut(km$lo_predict_2920_scale, 
                                  breaks = quantile(km$lo_predict_2920_scale, probs = seq(0, 1, by = 0.25), na.rm = TRUE), 
                                  include.lowest = TRUE, 
                                  labels = c("Q1", "Q2", "Q3", "Q4"))
km$mi_apo_time_y <- km$mi_apo_time/365

# ----- p1: lo_predict_248_scale_q -----
log_rank_test <- survdiff(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_248_scale_q, data = km, na.action = na.exclude)
p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

fit1 <- survfit(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_248_scale_q,
                data = km, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain", na.action = na.exclude)

cox_model1 <- coxph(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_248_scale_q + LTALTER + LCSEX + 
                      current_smoking + LL_HDLN + LL_CHOLN + LTSYSMM, data = km)
cox_summary1 <- summary(cox_model1)
coef_indices_1 <- grep("lo_predict_248_scale_q", rownames(cox_summary1$coefficients))
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
  fun = "cumhaz",
  pval = FALSE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ncensor.plot = FALSE,
  censor.shape = NA,
  linetype = "solid",
  ggtheme = base_theme,
  palette = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99"),
  xlab = "Time (years)",
  ylab = "Cumulative incidence rate of MI/Stroke (%)",
  legend.labs = legend.labs1,
  risk.table.font = font_size * 0.5,
  legend.title = "Restricted ArtheroBurden Arterial",
  xlim = c(0, 16),
  break.time.by = 4,
  tables.theme = theme_survminer(base_size = font_size)
)

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.0001, " < 0.0001", 
                       paste0(" = ", round(p.val, 3))))

p1_surv$plot <- p1_surv$plot +
  scale_y_continuous(labels = function(x) x * 100) +
  annotate("text", x = 0, y = 0.2, label = p.lab, hjust = 0, fontface = "bold", size = font_size/3) +
  theme(
    legend.position = c(0.02, 0.98), 
    legend.justification = c("left", "top"), 
    legend.direction = "vertical",
    plot.title = element_text(size = font_size, face = "bold")
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
    axis.text.y = element_markdown(size = font_size),

    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

# ----- p2: Atherosclerosis_680_quartile -----
log_rank_test <- survdiff(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_680_scale_q, data = km, na.action = na.exclude)
p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

fit1 <- survfit(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_680_scale_q,
                data = km, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain", na.action = na.exclude)

cox_model1 <- coxph(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_680_scale_q + LTALTER + LCSEX + 
                      current_smoking + LL_HDLN + LL_CHOLN + LTSYSMM, data = km)
cox_summary1 <- summary(cox_model1)
coef_indices_1 <- grep("lo_predict_680_scale_q", rownames(cox_summary1$coefficients))
hr_vals_1 <- exp(cox_summary1$coefficients[coef_indices_1, "coef"])
ci_vals_1 <- cox_summary1$conf.int[coef_indices_1, c("lower .95", "upper .95")]

legend.labs1 <- c(
  "Q1, Reference",
  paste0("Q2, HR = ", round(hr_vals_1[1], 2), " (95% CI, ", round(ci_vals_1[1,1], 2), "-", round(ci_vals_1[1,2], 2), ")"),
  paste0("Q3, HR = ", round(hr_vals_1[2], 2), " (95% CI, ", round(ci_vals_1[2,1], 2), "-", round(ci_vals_1[2,2], 2), ")"),
  paste0("Q4, HR = ", round(hr_vals_1[3], 2), " (95% CI, ", round(ci_vals_1[3,1], 2), "-", round(ci_vals_1[3,2], 2), ")")
)

p2_surv <- ggsurvplot(
  fit1,
  fun = "cumhaz",
  pval = FALSE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ncensor.plot = FALSE,
  censor.shape = NA,
  linetype = "solid",
  ggtheme = base_theme,
  palette = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99"),
  xlab = "Time (years)",
  ylab = "Cumulative incidence rate of MI/Stroke (%)",
  legend.labs = legend.labs1,
  risk.table.font = font_size * 0.5,
  legend.title = "Restricted ArtheroBurden Mechanistic",
  xlim = c(0, 16),
  break.time.by = 4,
  tables.theme = theme_survminer(base_size = font_size)
)

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.0001, " < 0.0001", 
                       paste0(" = ", round(p.val, 3))))

p2_surv$plot <- p2_surv$plot +
  scale_y_continuous(labels = function(x) x * 100) +
  annotate("text", x = 0, y = 0.2, label = p.lab, hjust = 0, fontface = "bold", size = font_size/3) +
  theme(
    legend.position = c(0.02, 0.98), 
    legend.justification = c("left", "top"), 
    legend.direction = "vertical",
    plot.title = element_text(size = font_size, face = "bold")
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
    axis.text.y = element_markdown(size = font_size),

    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

# ----- p3: MR_derived_402_quartile -----
log_rank_test <- survdiff(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_402_scale_q, data = km, na.action = na.exclude)
p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

fit1 <- survfit(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_402_scale_q,
                data = km, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain", na.action = na.exclude)

cox_model1 <- coxph(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_402_scale_q + LTALTER + LCSEX + 
                      current_smoking + LL_HDLN + LL_CHOLN + LTSYSMM, data = km)
cox_summary1 <- summary(cox_model1)
coef_indices_1 <- grep("lo_predict_402_scale_q", rownames(cox_summary1$coefficients))
hr_vals_1 <- exp(cox_summary1$coefficients[coef_indices_1, "coef"])
ci_vals_1 <- cox_summary1$conf.int[coef_indices_1, c("lower .95", "upper .95")]

legend.labs1 <- c(
  "Q1, Reference",
  paste0("Q2, HR = ", round(hr_vals_1[1], 2), " (95% CI, ", round(ci_vals_1[1,1], 2), "-", round(ci_vals_1[1,2], 2), ")"),
  paste0("Q3, HR = ", round(hr_vals_1[2], 2), " (95% CI, ", round(ci_vals_1[2,1], 2), "-", round(ci_vals_1[2,2], 2), ")"),
  paste0("Q4, HR = ", round(hr_vals_1[3], 2), " (95% CI, ", round(ci_vals_1[3,1], 2), "-", round(ci_vals_1[3,2], 2), ")")
)

p3_surv <- ggsurvplot(
  fit1,
  fun = "cumhaz",
  pval = FALSE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ncensor.plot = FALSE,
  censor.shape = NA,
  linetype = "solid",
  ggtheme = base_theme,
  palette = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99"),
  xlab = "Time (years)",
  ylab = "Cumulative incidence rate of MI/Stroke (%)",
  legend.labs = legend.labs1,
  risk.table.font = font_size * 0.5,
  legend.title = "Restricted ArtheroBurden Genetic",
  xlim = c(0, 16),
  break.time.by = 4,
  tables.theme = theme_survminer(base_size = font_size)
)

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.0001, " < 0.0001", 
                       paste0(" = ", round(p.val, 3))))

p3_surv$plot <- p3_surv$plot +
  scale_y_continuous(labels = function(x) x * 100) +
  annotate("text", x = 0, y = 0.2, label = p.lab, hjust = 0, fontface = "bold", size = font_size/3) +
  theme(
    legend.position = c(0.02, 0.98), 
    legend.justification = c("left", "top"), 
    legend.direction = "vertical",
    plot.title = element_text(size = font_size, face = "bold")
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
    axis.text.y = element_markdown(size = font_size),

    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

# ----- p4: Whole_proteome_2920_quartile -----
log_rank_test <- survdiff(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_2920_scale_q, data = km, na.action = na.exclude)
p.val <- 1 - pchisq(log_rank_test$chisq, length(log_rank_test$n) - 1)

fit1 <- survfit(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_2920_scale_q,
                data = km, type = "kaplan-meier", error = "greenwood",
                conf.type = "plain", na.action = na.exclude)

cox_model1 <- coxph(Surv(mi_apo_time_y, inz_mi_apo) ~ lo_predict_2920_scale_q + LTALTER + LCSEX + 
                      current_smoking + LL_HDLN + LL_CHOLN + LTSYSMM, data = km)
cox_summary1 <- summary(cox_model1)
coef_indices_1 <- grep("lo_predict_2920_scale_q", rownames(cox_summary1$coefficients))
hr_vals_1 <- exp(cox_summary1$coefficients[coef_indices_1, "coef"])
ci_vals_1 <- cox_summary1$conf.int[coef_indices_1, c("lower .95", "upper .95")]

legend.labs1 <- c(
  "Q1, Reference",
  paste0("Q2, HR = ", round(hr_vals_1[1], 2), " (95% CI, ", round(ci_vals_1[1,1], 2), "-", round(ci_vals_1[1,2], 2), ")"),
  paste0("Q3, HR = ", round(hr_vals_1[2], 2), " (95% CI, ", round(ci_vals_1[2,1], 2), "-", round(ci_vals_1[2,2], 2), ")"),
  paste0("Q4, HR = ", round(hr_vals_1[3], 2), " (95% CI, ", round(ci_vals_1[3,1], 2), "-", round(ci_vals_1[3,2], 2), ")")
)

p4_surv <- ggsurvplot(
  fit1,
  fun = "cumhaz",
  pval = FALSE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  ncensor.plot = FALSE,
  censor.shape = NA,
  linetype = "solid",
  ggtheme = base_theme,
  palette = c("#e66101", "#fdb863", "#b2abd2", "#5e3c99"),
  xlab = "Time (years)",
  ylab = "Cumulative incidence rate of MI/Stroke (%)",
  legend.labs = legend.labs1,
  risk.table.font = font_size * 0.5,
  legend.title = "Restricted ArtheroBurden WholeProteome",
  xlim = c(0, 16),
  break.time.by = 4,
  tables.theme = theme_survminer(base_size = font_size)
)

p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.0001, " < 0.0001", 
                       paste0(" = ", round(p.val, 3))))

p4_surv$plot <- p4_surv$plot +
  scale_y_continuous(labels = function(x) x * 100) +
  annotate("text", x = 0, y = 0.2, label = p.lab, hjust = 0, fontface = "bold", size = font_size/3) +
  theme(
    legend.position = c(0.02, 0.98), 
    legend.justification = c("left", "top"), 
    legend.direction = "vertical",
    plot.title = element_text(size = font_size, face = "bold")
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
    axis.text.y = element_markdown(size = font_size),

    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )


library(patchwork)
library(gridExtra)


font_size <- 16
border_width <- 0.8



p_arterial_with_legend <- p_arterial + 
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    legend.key.width = unit(2, "cm")
  )


png("figures/Figure6_hex_legend.png", width = 800, height = 100, res = 100)
tmp <- ggplot_gtable(ggplot_build(p_arterial_with_legend))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
if(length(leg) > 0) {
  legend <- tmp$grobs[[leg]]
  grid::grid.draw(legend)
}
dev.off()



get_legend <- function(a_plot) {
  tmp <- ggplotGrob(a_plot + 
                      scale_color_manual(
                        values = c("model1" = "#e66101", "model2" = "#fdb863", "model3" = "#5e3c99"),
                        labels = c("Adjusted for Age and Sex", "Adjusted for SCORE2 Variables", "Adjusted for VRFs"),
                        name = NULL
                      ) + 
                      theme(legend.position = "bottom"))
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

forest_legend <- get_legend(plot1)


hex_legend_img <- png::readPNG("figures/Figure6_hex_legend.png")
hex_legend <- rasterGrob(hex_legend_img)


p_arterial <- p_arterial + theme(legend.position = "none")
p_mech <- p_mech + theme(legend.position = "none")
p_genetic <- p_genetic + theme(legend.position = "none")
p_fullprot <- p_fullprot + theme(legend.position = "none")

plot1 <- plot1 + theme(legend.position = "none")
plot2 <- plot2 + theme(legend.position = "none")
plot3 <- plot3 + theme(legend.position = "none")
plot4 <- plot4 + theme(legend.position = "none")



col1 <- (p_arterial / plot_spacer() / plot1 / plot_spacer() / p1_surv$plot / p1_surv$table) +
  plot_layout(heights = c(1, 0.2, 1, 0.2, 1.6, 0.3))


col2 <- (p_mech / plot_spacer() / plot2 / plot_spacer() / p2_surv$plot / p2_surv$table) + 
  plot_layout(heights = c(1, 0.3, 1, 0.2, 1.6, 0.4))


col3 <- (p_genetic / plot_spacer() / plot3 / plot_spacer() / p3_surv$plot / p3_surv$table) + 
  plot_layout(heights = c(1, 0.3, 1, 0.2, 1.6, 0.4))


col4 <- (p_fullprot / plot_spacer() / plot4 / plot_spacer() / p4_surv$plot / p4_surv$table) + 
  plot_layout(heights = c(1, 0.3, 1, 0.2, 1.6, 0.4))



final_plot <- ((col1 | col2 | col3 | col4) / hex_legend / forest_legend) + 
  plot_layout(heights = c(6, 0.3, 0.3), widths = c(1, 1, 1, 1))


print(final_plot)
library(export)
graph2ppt(final_plot, file = "figures/Figure6.pptx", width = 32, height = 30)


ggsave("figures/final_aligned_plot.pdf", final_plot, width = 32, height = 28, limitsize = FALSE)
