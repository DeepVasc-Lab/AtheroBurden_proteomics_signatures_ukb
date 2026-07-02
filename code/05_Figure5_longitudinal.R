# Input note:
# This script starts from anonymous longitudinal analysis-ready files.
#
# Main input files:
# - data/data_for_MACE_barplot.csv:
#   long-format AtheroBurden score trajectories stratified by later MACE status.
# - data/data_for_SCORE2_barplot.csv:
#   long-format AtheroBurden score trajectories stratified by baseline SCORE2 risk group.
# - data/restricted_score_correlation.csv:
#   paired full-panel and restricted-panel AtheroBurden scores for correlation plots.
#
# Required columns for the two longitudinal trajectory files:
# - eid: anonymous participant ID.
# - Instance: visit/assessment instance.
# - Years_since_baseline: time from baseline visit in years.
# - Score: AtheroBurden score value.
# - score_type: panel identifier, one of 248, 680, 402, or 2920.
# - score2_group: SCORE2 risk group, required for data_for_SCORE2_barplot.csv.
# - MACE: later MACE status, required for data_for_MACE_barplot.csv.
#
# Longitudinal model used for beta labels:
# Score ~ Years_since_baseline + (1 | eid)
# This estimates annual score change while allowing each participant to have
# their own baseline score through a random intercept.

library(data.table)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(cowplot)
library(export)

my_data <- fread("data/data_for_MACE_barplot.csv", data.table = FALSE)
data_for_SCORE2_barplot <- fread("data/data_for_SCORE2_barplot.csv", data.table = FALSE)
analysis_df <- fread("data/restricted_score_correlation.csv", data.table = FALSE)
# -------------------------------

common_theme <- theme_classic(base_size = 14) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    legend.position = "right",
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.subtitle = element_text(size = 14, face = "bold")
  )


bar_theme <- theme_classic(base_size = 14) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    legend.position = "right",
    text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.subtitle = element_text(size = 14, face = "bold")
  )

# -------------------------------

create_bar_plot_with_beta_SCORE2 <- function(data, score_type_value, subtitle_text = NULL, y_max = NULL) {
  data_filtered <- data %>%
    filter(score_type == score_type_value) %>%
    mutate(
      score2_group = factor(score2_group, levels = c("<2.5%", "2.5-5%", "5-7.5%", ">7.5%")),
      Instance = factor(Instance)
    )
  
  if(nrow(data_filtered) == 0) {
    message(paste("Warning: score_type =", score_type_value, "has no matching rows in the data"))
    return(ggplot() + labs(subtitle = paste0(subtitle_text, " (no data)")) + theme_void())
  }
  
  beta_results <- data_filtered %>%
    group_by(score2_group) %>%
    group_modify(~{
      if(nrow(.x) < 2 || all(is.na(.x$Score))) {
        return(tibble(beta = NA_real_, p_value = NA_real_))
      }
      model <- lmer(Score ~ Years_since_baseline + (1 | eid), data = .x)
      coefs <- summary(model)$coefficients
      tibble(
        beta    = coefs["Years_since_baseline", "Estimate"],
        p_value = coefs["Years_since_baseline", "Pr(>|t|)"]
      )
    })
  
  results <- data_filtered %>%
    group_by(score2_group, Instance) %>%
    summarise(
      mean = mean(Score, na.rm = TRUE),
      se   = sd(Score, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  if(nrow(results) == 0) {
    return(ggplot() + labs(subtitle = paste0(subtitle_text, " (no valid data)")) + theme_void())
  }
  
  local_y_max <- max(results$mean + results$se, na.rm = TRUE)
  y_max_to_use <- if (!is.null(y_max)) y_max else local_y_max
  y_range <- local_y_max - min(results$mean - results$se, na.rm = TRUE)
  
  beta_labels <- beta_results %>%
    mutate(
      p_str = if_else(!is.na(p_value) & p_value >= 0.001,
                      sprintf("%.3f", p_value),
                      format(p_value, scientific = TRUE, digits = 3)),
      label = sprintf("β = %.3f/year\np = %s", beta, p_str),
      y_pos = y_max_to_use + 0.15 * y_range
    )
  
  p <- ggplot(results, aes(x = score2_group, y = mean)) +
    geom_bar(aes(fill = Instance), stat = "identity",
             position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se, group = Instance),
                  position = position_dodge(width = 0.8), width = 0.25) +
    geom_text(data = beta_labels, aes(x = score2_group, y = y_pos, label = label),
              size = 4, inherit.aes = FALSE) +
    scale_fill_brewer(palette = "BuPu") +
    labs(
      subtitle = subtitle_text,
      x = "SCORE2 absolute risk",
      y = "Mean Score (±SE)",
      fill = "Instance"
    ) +
    scale_y_continuous(limits = c(-0.5, y_max_to_use + 0.2), expand = expansion(mult = c(0.1, 0.1))) +
    bar_theme
  
  return(p)
}

# -------------------------------

create_bar_plot_with_beta_MACE <- function(data, score_type_value, subtitle_text = NULL, y_max = NULL) {
  data_filtered <- data %>%
    filter(score_type == score_type_value) %>%
    mutate(
      MACE = factor(MACE, levels = c("0", "1"), labels = c("No", "Yes")),
      Instance = factor(Instance)
    )
  
  if(nrow(data_filtered) == 0) {
    message(paste("Warning: score_type =", score_type_value, "has no matching rows in the data"))
    return(ggplot() + labs(subtitle = paste0(subtitle_text, " (no data)")) + theme_void())
  }
  
  beta_results <- data_filtered %>%
    group_by(MACE) %>%
    group_modify(~{
      valid_data <- .x %>% filter(!is.na(Score))
      if(nrow(valid_data) < 2) {
        return(tibble(beta = NA_real_, p_value = NA_real_))
      }
      model <- lmer(Score ~ Years_since_baseline + (1 | eid), data = valid_data)
      coefs <- summary(model)$coefficients
      tibble(
        beta = coefs["Years_since_baseline", "Estimate"],
        p_value = coefs["Years_since_baseline", "Pr(>|t|)"]
      )
    })
  
  results <- data_filtered %>%
    group_by(MACE, Instance) %>%
    summarise(
      mean = mean(Score, na.rm = TRUE),
      se   = sd(Score, na.rm = TRUE) / sqrt(sum(!is.na(Score))),
      .groups = 'drop'
    )
  
  if(nrow(results) == 0) {
    return(ggplot() + labs(subtitle = paste0(subtitle_text, " (no valid data)")) + theme_void())
  }
  
  local_y_max <- max(results$mean + results$se, na.rm = TRUE)
  y_max_to_use <- if (!is.null(y_max)) y_max else local_y_max
  y_range <- local_y_max - min(results$mean - results$se, na.rm = TRUE)
  
  beta_labels <- beta_results %>%
    mutate(
      p_str = if_else(!is.na(p_value) & p_value >= 0.001,
                      sprintf("%.3f", p_value),
                      format(p_value, scientific = TRUE, digits = 3)),
      label = sprintf("β = %.3f/year\np = %s", beta, p_str),
      y_pos = y_max_to_use + 0.15 * y_range
    )
  
  p <- ggplot(results, aes(x = MACE, y = mean)) +
    geom_bar(aes(fill = Instance), stat = "identity",
             position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se, group = Instance),
                  position = position_dodge(width = 0.8), width = 0.25) +
    geom_text(data = beta_labels, aes(x = MACE, y = y_pos, label = label),
              size = 4, inherit.aes = FALSE) +
    scale_fill_brewer(palette = "BuPu") +
    labs(
      subtitle = subtitle_text,
      x = "MACE",
      y = "Mean Score (±SE)",
      fill = "Instance"
    ) +
    scale_y_continuous(limits = c(-0.5, y_max_to_use + 0.3 * y_range), expand = expansion(mult = c(0.1, 0.1))) +
    coord_cartesian(clip = "off") +
    bar_theme
  
  return(p)
}

# -------------------------------

score_types <- c("248", "680", "402", "2920")
score_subtitles <- c(
  "248"  = "AtheroBurden Arterial Signature",
  "680"  = "AtheroBurden Mechanistic Signature",
  "402"  = "AtheroBurden Genetic Signature",
  "2920" = "AtheroBurden WholeProteome Signature"
)

# -------------------------------

y_max_global_MACE <- 0
for (st in score_types) {
  tmp <- my_data %>%
    filter(score_type == st) %>%
    mutate(
      MACE = factor(MACE, levels = c("0", "1"), labels = c("No", "Yes")),
      Instance = factor(Instance)
    ) %>%
    group_by(MACE, Instance) %>%
    summarise(
      mean = mean(Score, na.rm = TRUE),
      se = sd(Score, na.rm = TRUE) / sqrt(sum(!is.na(Score))),
      .groups = "drop"
    )
  if(nrow(tmp) > 0) {
    local_max <- max(tmp$mean + tmp$se, na.rm = TRUE)
    y_max_global_MACE <- max(y_max_global_MACE, local_max, na.rm = TRUE)
  }
}

y_max_global_SCORE2 <- 0
for (st in score_types) {
  tmp <- data_for_SCORE2_barplot %>%
    filter(score_type == st) %>%
    mutate(
      score2_group = factor(score2_group, levels = c("<2.5%", "2.5-5%", "5-7.5%", ">7.5%")),
      Instance = factor(Instance)
    ) %>%
    group_by(score2_group, Instance) %>%
    summarise(
      mean = mean(Score, na.rm = TRUE),
      se = sd(Score, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  if(nrow(tmp) > 0) {
    local_max <- max(tmp$mean + tmp$se, na.rm = TRUE)
    y_max_global_SCORE2 <- max(y_max_global_SCORE2, local_max, na.rm = TRUE)
  }
}

# -------------------------------


mace_plots <- lapply(score_types, function(st) {
  create_bar_plot_with_beta_MACE(my_data, score_type_value = st,
                                 subtitle_text = score_subtitles[st],
                                 y_max = 2)
})


score2_plots <- lapply(score_types, function(st) {
  create_bar_plot_with_beta_SCORE2(data_for_SCORE2_barplot, score_type_value = st,
                                   subtitle_text = score_subtitles[st],
                                   y_max = y_max_global_SCORE2)
})

bar_plots_combined <- (wrap_plots(mace_plots, ncol = 4)) /
  (wrap_plots(score2_plots, ncol = 4)) +
  plot_layout(guides = "collect")

bar_plots_combined <- wrap_plots(
  list(
    wrap_plots(mace_plots, ncol = 4),
    plot_spacer(),
    wrap_plots(score2_plots, ncol = 4)
  ),
  ncol = 1,
  heights = c(0.9, 0.1, 1.1)
) + plot_layout(guides = "collect")

# -------------------------------

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



p_arterial <- plot_hex_lm_r2(
  data = analysis_df,
  xvar = "Artery_enriched_248",
  yvar = "pred_artery",
  xlab = "Full AtheroBurden Arterial Signature",
  ylab = "Restricted AtheroBurden Arterial Signature",
  bins = 40,
  maxcount = globalMaxCount,
  xlim = c(global_xmin, global_xmax)
)

p_mech <- plot_hex_lm_r2(
  data = analysis_df,
  xvar = "Atherosclerosis_680",
  yvar = "pred_arthero",
  xlab = "Full AtheroBurden Mechanistic Signature",
  ylab = "Restricted AtheroBurden Mechanistic Signature",
  bins = 40,
  maxcount = globalMaxCount,
  xlim = c(global_xmin, global_xmax)
)

p_genetic <- plot_hex_lm_r2(
  data = analysis_df,
  xvar = "MR_derived_402",
  yvar = "pred_mr",
  xlab = "Full AtheroBurden Genetic Signature",
  ylab = "Restricted AtheroBurden Genetic Signature",
  bins = 40,
  maxcount = globalMaxCount,
  xlim = c(global_xmin, global_xmax)
)

p_fullprot <- plot_hex_lm_r2(
  data = analysis_df,
  xvar = "Whole_proteome_2920",
  yvar = "pred_whole",
  xlab = "Full AtheroBurden WholeProteome Signature",
  ylab = "Restricted AtheroBurden WholeProteome Signature",
  bins = 40,
  maxcount = globalMaxCount,
  xlim = c(global_xmin, global_xmax)
)

hex_plots_combined <- (p_arterial | p_mech | p_genetic | p_fullprot) +
  plot_layout(guides = "collect")

# -------------------------------


final_combined_plot <- wrap_plots(
  list(
    hex_plots_combined,
    plot_spacer(),
    bar_plots_combined
  ),
  ncol = 1,
  heights = c(0.9, 0.1, 2)
) +
  plot_layout(guides = "collect")

print(final_combined_plot)

graph2ppt(final_combined_plot, file = "figures/Figure5.pptx", width = 24, height = 16)

ggsave("figures/Final_Combined_Plot.png", final_combined_plot, width = 24, height = 16, dpi = 300)
