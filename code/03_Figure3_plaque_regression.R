# Input note:
# This script starts from an anonymous plaque analysis-ready file.
# The file combines AtheroBurden scores, carotid ultrasound plaque outcomes,
# and covariates for the UK Biobank plaque subset.
#
# Main required columns:
# - predictions_248, predictions_680, predictions_402, predictions_2920:
#   AtheroBurden score columns from the four CatBoost panels.
# - presence: binary carotid plaque presence outcome.
# - plaque_number: carotid plaque burden/count outcome.
# - Age, Sex, ethnicity, Cholesterol, HDL, SBP, BMI, LDL, Triglycerides,
#   eGFRcr_cys, HbA1c, Diabetes_diag, Hypertension_diag:
#   clinical covariates used in adjusted models.
# - Current_smoke, Current_drink, Unknown_drink, fh_cvd:
#   smoking, drinking, and family-history covariates.

score_plaque_cli = fread('data/score_plaque_cli.csv',data.table = F)
colnames(score_plaque_cli)[1] = "eid_b"
head(score_plaque_cli)

map_b_to_151 <- bridge %>%
  distinct(eid_b, eid_m) %>%
  left_join(bridge_new %>% distinct(eid_m, eid_151281),
            by = "eid_m")
score_plaque_cli <- score_plaque_cli %>%
  left_join(map_b_to_151, by = "eid_b")
score_plaque_cli <- score_plaque_cli %>%
  mutate(eid_151281 = as.numeric(eid_151281)) %>%
  left_join(cov_add, by = "eid_151281")
head(score_plaque_cli)



library(dplyr)
library(purrr)
library(sandwich)
library(lmtest)

dat <- score_plaque_cli

#----------------------------

#----------------------------
dat <- dat %>%
  mutate(
    Sex = factor(Sex),
    ethnicity = as.character(ethnicity),
    ethnicity = ifelse(is.na(ethnicity) | ethnicity == "", "Unknown", ethnicity),
    ethnicity = factor(ethnicity, levels = c("White","Mixed","Asian","Black","Other","Unknown")),
    
    fh_cvd = as.character(fh_cvd),
    fh_cvd = ifelse(is.na(fh_cvd) | fh_cvd == "", "Unknown", fh_cvd),
    fh_cvd = factor(fh_cvd, levels = c("No","Yes","Unknown")),
    
    Current_smoke = ifelse(is.na(Current_smoke), 0L, as.integer(Current_smoke)),
    Current_drink = ifelse(is.na(Current_drink), 0L, as.integer(Current_drink)),
    Unknown_drink = ifelse(is.na(Unknown_drink), 0L, as.integer(Unknown_drink)),
    

    Diabetes_diag = factor(Diabetes_diag),
    Hypertension_diag = factor(Hypertension_diag)
  )


if (!("Unknown_smoke" %in% names(dat))) {
  if ("Smoking_status" %in% names(dat)) {

    dat <- dat %>% mutate(
      Unknown_smoke = as.integer(is.na(Smoking_status) | !(Smoking_status %in% c(1,2,3)))
    )
  } else {
    dat <- dat %>% mutate(Unknown_smoke = 0L)
  }
} else {
  dat <- dat %>% mutate(Unknown_smoke = ifelse(is.na(Unknown_smoke), 0L, as.integer(Unknown_smoke)))
}

#----------------------------

#----------------------------
score_map <- tibble::tribble(
  ~score_var,         ~prob_rs,
  "predictions_248",  "AtheroBurden Arterial signature",
  "predictions_680",  "AtheroBurden Mechanistic signature",
  "predictions_402",  "AtheroBurden Genetic signature",
  "predictions_2920", "AtheroBurden WholeProteome signature"
)


use_std <- TRUE
if (use_std) {
  for (v in score_map$score_var) {
    dat[[paste0(v, "_std")]] <- as.numeric(scale(dat[[v]]))
  }
  score_map <- score_map %>%
    mutate(score_var_std = paste0(score_var, "_std"))
} else {
  score_map <- score_map %>%
    mutate(score_var_std = score_var)
}

#----------------------------

#----------------------------
model_covars <- list(
  model1 = c("Age", "Sex", "ethnicity"),
  model2 = c("Age", "Sex", "Cholesterol", "HDL", "SBP", "Current_smoke", "Unknown_smoke"),
  model3 = c("Age", "Sex", "SBP", "BMI",
             "Current_smoke", "Unknown_smoke",
             "LDL", "Triglycerides", "eGFRcr_cys", "HbA1c",
             "Diabetes_diag", "Hypertension_diag",
             "ethnicity", "Current_drink", "Unknown_drink", "fh_cvd")
)

#----------------------------

#----------------------------
extract_glm <- function(d, outcome, score, covars, family, robust = FALSE) {
  vars_need <- c(outcome, score, covars)
  d2 <- d %>% filter(if_all(all_of(vars_need), ~ !is.na(.x)))
  
  fml <- as.formula(paste(outcome, "~", paste(c(score, covars), collapse = " + ")))
  fit <- glm(fml, data = d2, family = family)
  
  if (!robust) {
    sm <- summary(fit)$coefficients
    beta <- sm[score, "Estimate"]
    se   <- sm[score, "Std. Error"]
    pval <- sm[score, "Pr(>|z|)"]
  } else {
    V <- sandwich::vcovHC(fit, type = "HC0")
    ct <- lmtest::coeftest(fit, vcov. = V)
    beta <- ct[score, 1]
    se   <- ct[score, 2]
    pval <- ct[score, 4]
  }
  
  est <- exp(beta)
  lo  <- exp(beta - 1.96 * se)
  hi  <- exp(beta + 1.96 * se)
  
  list(estimate = est, ci_lower = lo, ci_upper = hi, p_value = pval)
}

#----------------------------

#----------------------------
run_block <- function(model_type, outcome, family, estimate_type, robust, model_covars, score_map, dat) {
  expand.grid(
    model = names(model_covars),
    score_idx = seq_len(nrow(score_map)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      model_num = as.numeric(sub("model", "", model)),
      model_type = model_type,
      estimate_type = estimate_type
    ) %>%
    pmap_dfr(function(model, score_idx, model_num, model_type, estimate_type) {
      sc <- score_map$score_var_std[score_idx]
      pr <- score_map$prob_rs[score_idx]
      covs <- model_covars[[model]]
      
      out <- extract_glm(
        d = dat,
        outcome = outcome,
        score = sc,
        covars = covs,
        family = family,
        robust = robust
      )
      
      tibble(
        model_type = factor(model_type),
        model = model,
        prob_rs = pr,
        estimate_type = estimate_type,
        estimate = out$estimate,
        ci_lower = out$ci_lower,
        ci_upper = out$ci_upper,
        p_value = out$p_value,
        model_num = model_num
      )
    })
}


res_logit <- run_block(
  model_type = "Logistic regression",
  outcome = "presence",
  family = binomial(link = "logit"),
  estimate_type = "OR",
  robust = FALSE,
  model_covars = model_covars,
  score_map = score_map,
  dat = dat
)


res_pois <- run_block(
  model_type = "Poisson regression",
  outcome = "plaque_number",
  family = poisson(link = "log"),
  estimate_type = "RR",
  robust = TRUE,
  model_covars = model_covars,
  score_map = score_map,
  dat = dat
)

final_results <- bind_rows(res_logit, res_pois) %>%
  group_by(model_type, model) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "fdr"),
    significance = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(model_type, model, prob_rs, estimate_type,
         estimate, ci_lower, ci_upper, p_value, p_adjusted,
         significance, model_num)

final_results
write.csv(final_results, "plaque_presence_logit_and_plaqueNumber_poisson_results.csv", row.names = FALSE)
final_results$model_type

data_Artery = final_results[final_results$prob_rs == "AtheroBurden Arterial signature",]
data_Atherosclerosis = final_results[final_results$prob_rs == "AtheroBurden Mechanistic signature",]
data_MR = final_results[final_results$prob_rs == "AtheroBurden Genetic signature",]
data_Whole = final_results[final_results$prob_rs == "AtheroBurden WholeProteome signature",]
library(ggplot2)
library(dplyr)
dput(colnames(final_results))
final_results <- final_results %>%
  mutate(model_type = ifelse(estimate_type == "OR", "Logistic regression", "Poisson regression"))
final_results$model_type <- factor(final_results$model_type, 
                                   levels = c("Poisson regression", "Logistic regression"))



create_forest_plot <- function(data, title) {
  # Ensure model_type is a factor and properly ordered
  data$model_type <- factor(data$model_type, levels = c("Poisson regression", "Logistic regression"))
  
  # Add y_position offset
  data$y_position <- as.numeric(data$model_type) + 
    ifelse(data$model == "model1", 0.3,
           ifelse(data$model == "model2", 0,
                  ifelse(data$model == "model3", -0.3, 0)))
  
  # Define background shadow data
  background_data <- data.frame(
    model_type = c("Poisson regression", "Logistic regression"),
    ymin = c(0.5, 1.5),
    ymax = c(1.5, 2.5)
  )
  
  # Build the plot
  p <- ggplot() +
    geom_rect(
      data = background_data,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = model_type),
      alpha = 0.1
    ) +  # Add background shadow
    geom_errorbarh(
      data = data,
      aes(xmin = ci_lower, xmax = ci_upper, y = y_position),
      height = 0.1, color = "black"
    ) +  # Add error bars
    geom_point(
      data = data,
      aes(x = estimate, y = y_position, color = model),
      size = 8
    ) +  # Add points
    scale_y_continuous(
      breaks = c(1, 2),
      labels = c("Plaque burden\n(Poisson)", "Plaque presence\n(Logistic)"),
      limits = c(0.4, 2.6)
    )+  # Y-axis settings
    scale_fill_manual(
      values = c("Logistic regression" = "#bdbdbd", "Poisson regression" = "#e6eaf2"),
      guide = "none"
    ) +  # Background fill
    scale_color_manual(
      values = c("model1" = "#e66101", "model2" = "#fdb863", "model3" = "#5e3c99"),
      labels = c("Adjusted for Age, Sex and ethnicity", "Adjusted for SCORE2 Variables", "Adjusted for VRFs")
    ) +  # Point colors
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
    scale_x_continuous(
      breaks = c(0.8, 1.0, 1.5),
      limits = c(0.8, 1.65)
    ) +
    theme_classic(base_size = 20) +
    labs(title = title, x = "OR (Logistic) / RR (Poisson) (95% CI)", y = "") +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 20),
      legend.position = "none",
      plot.title = element_text(size = 20, face = "bold", hjust = 0),
      axis.line = element_line(size = 1, color = "black"),
      axis.ticks = element_line(size = 1, color = "black"),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 20, face = "bold")
    )
  # Add significance labels on the right side of each point with corresponding model colors
  p <- p + geom_text(data = data, aes(x = ci_upper + 0.05, y = y_position, label = significance, 
                                      color = model), 
                     size = 8, hjust = 0)  # Use HR or another valid column
  
  
  return(p)
}

# Create individual forest plots for different risk scores
plot1 <- create_forest_plot(data_Artery, "AtheroBurden Arterial signature")
plot2 <- create_forest_plot(data_Atherosclerosis, "AtheroBurden Mechanistic signature")
plot3 <- create_forest_plot(data_MR, "AtheroBurden Genetic signature")
plot4 <- create_forest_plot(data_Whole, "AtheroBurden WholeProteome signature")

# Extract the legend from one of the plots
get_legend <- function(a_plot) {
  tmp <- ggplotGrob(a_plot + theme(legend.position = "bottom"))
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

legend_plot <- get_legend(plot1)

final_plot <- grid.arrange(
  arrangeGrob(plot4, plot3, plot2,plot1,ncol = 2),
  legend_plot,
  heights = c(10, 4)
)
final_plot
# 20*20 
