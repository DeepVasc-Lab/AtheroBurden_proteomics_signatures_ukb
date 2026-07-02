# Complete Figure 2 plotting code from anonymous analysis-ready files.

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(pROC)
library(ggsignif)

data_dir <- "data"
figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# note:
# Figure 2A-C use output files generated from CatBoost model training.
# Figure 2D-E use a score file containing AtheroBurden signature from CatBoost model.

panel_info <- data.frame(
  panel = c("248", "680", "402", "2920"),
  score_col = c("predictions_248", "predictions_680", "predictions_402", "predictions_2920"),
  label = c(
    "Artery-enriched panel (n=248)",
    "Atherosclerosis-related panel (n=680)",
    "MR-derived panel (n=402)",
    "Whole proteome panel (n=2920)"
  ),
  signature_label = c(
    "AtheroBurden Arterial Signature",
    "AtheroBurden Mechanistic Signature",
    "AtheroBurden Genetic Signature",
    "AtheroBurden Complete Signature"
  ),
  stringsAsFactors = FALSE
)

read_input <- function(file) {
  fread(file.path(data_dir, file), data.table = FALSE)
}

save_figure <- function(plot, file, width, height) {
  ggsave(file.path(figures_dir, file), plot = plot, width = width, height = height, dpi = 300)
}

plot_confusion_matrix <- function(conf_matrix, title) {
  conf_matrix_percent <- round(conf_matrix / sum(conf_matrix) * 100, 2)
  colnames(conf_matrix) <- c("Positive", "Negative")
  rownames(conf_matrix) <- c("Positive", "Negative")
  conf_matrix_df <- as.data.frame(as.table(conf_matrix))
  conf_matrix_percent_df <- as.data.frame(as.table(conf_matrix_percent))
  conf_matrix_df$percent <- conf_matrix_percent_df$Freq
  conf_matrix_df$Var1 <- factor(conf_matrix_df$Var1, levels = rev(levels(conf_matrix_df$Var1)))

  ggplot(conf_matrix_df, aes(Var2, Var1, fill = Freq)) +
    geom_tile(color = "black", size = 0.6) +
    geom_text(aes(label = paste0(Freq, "\n(", percent, "%)")), color = "black", size = 6) +
    scale_fill_gradient(low = "white", high = "#cd4a33") +
    labs(x = "Actual Label", y = "Predicted Label") +
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold"),
      axis.text.y = element_text(angle = 90, hjust = 0.5, size = 18, face = "bold"),
      panel.grid = element_blank()
    )
}

confusion_matrices <- list(
  "248" = matrix(c(268, 89, 62, 248), nrow = 2, byrow = TRUE),
  "680" = matrix(c(281, 65, 49, 272), nrow = 2, byrow = TRUE),
  "402" = matrix(c(283, 67, 47, 270), nrow = 2, byrow = TRUE),
  "2920" = matrix(c(281, 63, 49, 274), nrow = 2, byrow = TRUE)
)

confusion_plots <- lapply(panel_info$panel, function(panel) {
  plot_confusion_matrix(
    confusion_matrices[[panel]],
    panel_info$label[panel_info$panel == panel]
  ) + theme(legend.position = "none")
})
figure2a <- grid.arrange(grobs = confusion_plots, nrow = 1)
save_figure(figure2a, "Figure2A_confusion_matrix.pdf", width = 22, height = 6)

create_roc_plot <- function(data, model_name, color, score_roc_df, score_auc_value, score_auc_ci) {
  roc_curve <- roc(data$True_Label, data$Predicted_Probability)
  roc_df <- data.frame(TPR = roc_curve$sensitivities, FPR = 1 - roc_curve$specificities)
  auc_value <- auc(roc_curve)
  auc_ci <- ci.auc(roc_curve, conf.level = 0.95)

  ggplot() +
    geom_line(data = roc_df, aes(x = FPR, y = TPR), color = color, size = 1.5) +
    geom_line(data = score_roc_df, aes(x = FPR, y = TPR), color = "#0072B2", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "darkgray") +
    labs(title = model_name, x = "False Positive Rate", y = "True Positive Rate") +
    theme_classic() +
    theme(
      axis.line = element_line(size = 0.6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0, size = 18, face = "bold"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 18)
    ) +
    annotate("text", x = 0.15, y = 0.12, label = "AUC (95% CI):",
             size = 4.5, color = "black", fontface = "bold", hjust = 0) +
    geom_segment(aes(x = 0.15, y = 0.06, xend = 0.25, yend = 0.06),
                 color = "#0072B2", size = 5) +
    annotate("text", x = 0.26, y = 0.06,
             label = paste("SCORE2 variables:", round(score_auc_value, 3),
                           "(", round(score_auc_ci[1], 3), "-", round(score_auc_ci[3], 3), ")"),
             size = 4.5, color = "black", fontface = "bold", hjust = 0) +
    geom_segment(aes(x = 0.15, y = 0, xend = 0.25, yend = 0),
                 color = color, size = 5) +
    annotate("text", x = 0.26, y = 0,
             label = paste("Protein Panel:", round(auc_value, 3),
                           "(", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), ")"),
             size = 4.5, color = "black", fontface = "bold", hjust = 0)
}

roc_files <- paste0("roc_data_", panel_info$panel, ".csv")
# ROC files are CatBoost model-training outputs.
if (all(file.exists(file.path(data_dir, roc_files))) &&
    file.exists(file.path(data_dir, "score2_roc_df_test.csv"))) {
  score2_roc_df <- read_input("score2_roc_df_test.csv")
  score2_roc_curve <- roc(score2_roc_df$True_Label, score2_roc_df$Predicted_Probability)
  score_roc_df <- data.frame(
    TPR = score2_roc_curve$sensitivities,
    FPR = 1 - score2_roc_curve$specificities
  )
  score_auc_value <- auc(score2_roc_curve)
  score_auc_ci <- ci.auc(score2_roc_curve, conf.level = 0.95)

  roc_plots <- lapply(seq_len(nrow(panel_info)), function(i) {
    dat <- read_input(roc_files[i])
    create_roc_plot(dat, panel_info$label[i], "#cd4a33", score_roc_df, score_auc_value, score_auc_ci)
  })
  figure2b <- grid.arrange(grobs = roc_plots, nrow = 1)
  save_figure(figure2b, "Figure2B_ROC.pdf", width = 22, height = 6)
}

create_shap_bar_plot <- function(data, title) {
  ggplot(data, aes(x = reorder(Feature, Mean_SHAP_Value), y = Mean_SHAP_Value, fill = Mean_SHAP_Value)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_gradient(low = "#0072B2", high = "#cd4a33") +
    labs(title = title, x = "Top 10 Feature", y = "Mean Absolute SHAP Value") +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(size = 0.6),
      plot.title = element_text(size = 18, face = "bold", hjust = 0),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 18, face = "bold"),
      legend.position = "none"
    )
}

shap_files <- paste0("mean_shap_values_", panel_info$panel, ".csv")
# SHAP files are CatBoost model-training outputs.
if (all(file.exists(file.path(data_dir, shap_files)))) {
  shap_plots <- lapply(seq_len(nrow(panel_info)), function(i) {
    shap_data <- read_input(shap_files[i])
    top10 <- shap_data[order(-shap_data$Mean_SHAP_Value), ][1:10, ]
    create_shap_bar_plot(top10, panel_info$label[i])
  })
  figure2c <- grid.arrange(grobs = shap_plots, nrow = 1)
  save_figure(figure2c, "Figure2C_SHAP_top10.pdf", width = 24, height = 6)
}

get_legend <- function(a_plot) {
  tmp <- ggplotGrob(a_plot + theme(legend.position = "bottom"))
  gtable::gtable_filter(tmp, "guide-box")
}

figure2_data_file <- file.path(data_dir, "figure2_plot_data.csv")
# This file should contain AtheroBurden scores for density and burden plots.
if (file.exists(figure2_data_file)) {
  figure2_data <- fread(figure2_data_file, data.table = FALSE)

  if (all(c("Group", panel_info$score_col) %in% colnames(figure2_data))) {
    create_density_plot <- function(score_column, title) {
      ggplot(figure2_data, aes(x = .data[[score_column]], fill = Group)) +
        geom_density(alpha = 0.5, color = "#000000") +
        scale_fill_manual(
          values = c("UKB" = "#756bb1", "Case" = "#cd4a33", "Control" = "#0072B2"),
          labels = c(
            "UKB" = "UKB cohort",
            "Case" = "Atherosclerosis case",
            "Control" = "Healthy control"
          )
        ) +
        labs(title = title, x = "AtheroBurden Score (Z-score)", y = "Density") +
        theme_classic(base_size = 18) +
        theme(
          axis.line = element_line(size = 0.6),
          plot.title = element_text(size = 18, face = "bold", hjust = 0),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 18),
          legend.title = element_blank(),
          legend.text = element_text(size = 18),
          legend.position = "none"
        )
    }
    density_plots <- lapply(seq_len(nrow(panel_info)), function(i) {
      create_density_plot(panel_info$score_col[i], panel_info$signature_label[i])
    })
    density_legend <- get_legend(density_plots[[1]])
    figure2d <- grid.arrange(arrangeGrob(grobs = density_plots, ncol = 4), density_legend,
                             nrow = 2, heights = c(10, 4))
    save_figure(figure2d, "Figure2D_density.pdf", width = 22, height = 6.4)
  }

  if (all(c("burden", panel_info$score_col) %in% colnames(figure2_data))) {
    score_long <- figure2_data %>%
      select(any_of(c("Participant ID", "participant_id")), burden, all_of(panel_info$score_col)) %>%
      mutate(burden = ifelse(burden == 3, 2, burden)) %>%
      pivot_longer(cols = all_of(panel_info$score_col), names_to = "Score_Type", values_to = "Score_Value")

    create_box_violin_plot <- function(score_type, title) {
      df_filtered <- score_long %>% filter(Score_Type == score_type)
      ggplot(df_filtered, aes(x = factor(burden), y = Score_Value, fill = factor(burden))) +
        geom_violin(alpha = 0.3, color = "black") +
        geom_boxplot(width = 0.3, alpha = 0.7, outlier.shape = 21,
                     outlier.size = 2, outlier.color = "black") +
        scale_fill_manual(
          values = c("0" = "#0072B2", "1" = "#756bb1", "2" = "#cd4a33"),
          name = "Atherosclerosis Burden Level",
          labels = c(
            "0" = "0: No Atherosclerosis",
            "1" = "1: Single Vascular Bed Involvement",
            "2" = "2: Two or More Vascular Beds Involvement"
          )
        ) +
        labs(title = title, x = "Atherosclerotic burden level",
             y = "AtheroBurden Score (Z-score)") +
        theme_classic(base_size = 18) +
        theme(
          axis.line = element_line(size = 0.6),
          plot.title = element_text(size = 18, face = "bold", hjust = 0),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold", angle = 45, hjust = 1),
          axis.text.y = element_text(size = 18, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 18),
          panel.grid.major = element_line(size = 0.5, linetype = "dotted", color = "grey"),
          legend.position = "bottom"
        ) +
        geom_signif(
          comparisons = list(c("0", "1"), c("0", "2"), c("1", "2")),
          map_signif_level = TRUE,
          test = "t.test",
          step_increase = 0.1,
          tip_length = 0.01
        )
    }

    burden_plots <- lapply(seq_len(nrow(panel_info)), function(i) {
      create_box_violin_plot(panel_info$score_col[i], panel_info$signature_label[i])
    })
    burden_legend <- get_legend(burden_plots[[4]])
    burden_plots <- lapply(burden_plots, function(p) p + theme(legend.position = "none"))
    figure2e <- grid.arrange(arrangeGrob(grobs = burden_plots, ncol = 4), burden_legend,
                             nrow = 2, heights = c(5, 2))
    save_figure(figure2e, "Figure2E_burden_violin.pdf", width = 24, height = 7.5)
  }
}
