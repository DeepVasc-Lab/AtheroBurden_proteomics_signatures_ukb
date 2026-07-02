olink = fread("data/OlinkNDimpute.csv",data.table = F)
colnames(olink)[1] = "eid_b"
olink_44788 = olink[,c(1,2)]

map_b_to_151 <- bridge %>%
  distinct(eid_b, eid_m) %>%
  left_join(bridge_new %>% distinct(eid_m, eid_151281),
            by = "eid_m")


olink_44788 <- olink_44788 %>%
  left_join(map_b_to_151, by = "eid_b")
olink_44788 <- olink_44788 %>%
  mutate(eid_151281 = as.numeric(eid_151281)) %>%
  left_join(cov_add, by = "eid_151281")


head(olink_44788)
library(dplyr)
library(tidyr)

fmt_n_pct <- function(n, denom, digits = 0){
  pct <- 100 * n / denom
  sprintf(paste0("%d (%.", digits, "f%%)"), n, pct)
}

make_cat_block <- function(df, var, header, denom = nrow(df), na_label = "Missing"){
  tmp <- df %>%
    mutate(.v = as.character(.data[[var]]),
           .v = replace_na(.v, na_label)) %>%
    count(.v, name = "n") %>%
    mutate(value = fmt_n_pct(n, denom)) %>%
    arrange(.v)
  
  bind_rows(
    tibble(Characteristic = header, UKB = ""),
    tibble(Characteristic = paste0("  ", tmp$.v), UKB = tmp$value)
  )
}
N <- nrow(olink_44788)

tab_fh  <- make_cat_block(olink_44788, "fh_cvd",
                          "Family history of CVD, n (%)", denom = N)

tab_eth <- make_cat_block(olink_44788, "ethnicity",
                          "Ethnicity, n (%)", denom = N)


tab_alc_status <- make_cat_block(olink_44788, "alcohol_status",
                                 "Alcohol drinking status, n (%)", denom = N)


tab_alc_prevcur <- tibble(
  Characteristic = c("Previous drinking, n (%)", "Current drinking, n (%)"),
  UKB = c(
    fmt_n_pct(sum(olink_44788$Previous_drink == 1, na.rm = TRUE), N),
    fmt_n_pct(sum(olink_44788$Current_drink  == 1, na.rm = TRUE), N)
  )
)


tab_additional <- bind_rows(
  tab_fh,
  tab_eth,
  tab_alc_status,
  tab_alc_prevcur
)

tab_additional
tab_eth <- make_cat_block(cov_add, "ethnicity",
                          "Ethnicity, n (%)", denom = N)
colnames(rm_FAS_new_PSM_index)[1] = "eid_b"
discovery_data = rm_FAS_new_PSM_index[,c(1,2)]

map_b_to_151 <- bridge %>%
  distinct(eid_b, eid_m) %>%
  left_join(bridge_new %>% distinct(eid_m, eid_151281),
            by = "eid_m")


discovery_data <- discovery_data %>%
  left_join(map_b_to_151, by = "eid_b")
discovery_data <- discovery_data %>%
  mutate(eid_151281 = as.numeric(eid_151281)) %>%
  left_join(cov_add, by = "eid_151281")
table(discovery_data$group, useNA = "ifany")
library(dplyr)

discovery_data <- discovery_data %>%
  mutate(
    group = factor(group,
                   levels = c(0, 1),
                   labels = c("No Atherosclerotic events", "Atherosclerotic events"))
  )

library(dplyr)
library(tidyr)

fmt_n_pct <- function(n, denom, digits = 1){
  pct <- 100 * n / denom
  sprintf(paste0("%d (%.", digits, "f%%)"), n, pct)
}

fmt_p <- function(p){
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

p_cat <- function(g, v){
  tab <- table(g, v)

  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  if (all(dim(tab) == c(2,2)) && any(chi$expected < 5)) {
    return(fisher.test(tab)$p.value)
  }
  chi$p.value
}


cat_block <- function(df, by, var, header, na_label = "Missing", digits = 1){
  df2 <- df %>%
    mutate(
      .g = .data[[by]],
      .v = as.character(.data[[var]]),
      .v = ifelse(is.na(.v) | .v == "", na_label, .v)
    )
  
  pval <- fmt_p(p_cat(df2$.g, df2$.v))
  
  # overall
  overall <- df2 %>%
    count(.v, name="n") %>%
    mutate(Overall = fmt_n_pct(n, nrow(df2), digits)) %>%
    select(.v, Overall)
  
  # by-group
  bytab <- df2 %>%
    count(.g, .v, name="n") %>%
    group_by(.g) %>%
    mutate(denom = sum(n),
           value = fmt_n_pct(n, denom, digits)) %>%
    ungroup() %>%
    select(.g, .v, value) %>%
    pivot_wider(names_from = .g, values_from = value)
  
  out <- overall %>%
    left_join(bytab, by = ".v") %>%
    mutate(Characteristic = paste0("  ", .v),
           `p value` = "") %>%
    select(Characteristic, Overall, everything(), `p value`)
  

  header_row <- out[0, ]
  header_row[1, "Characteristic"] <- header
  header_row[1, "p value"] <- pval
  header_row[, setdiff(names(out), c("Characteristic","p value"))] <- ""
  
  bind_rows(header_row, out)
}


bin_row <- function(df, by, var, label, digits = 1, na_as0 = TRUE){
  x <- df[[var]]
  if (na_as0) x[is.na(x)] <- 0
  df2 <- df %>% mutate(.g = .data[[by]], .x = x)
  
  overall <- fmt_n_pct(sum(df2$.x == 1), nrow(df2), digits)
  
  bytab <- df2 %>%
    group_by(.g) %>%
    summarise(value = fmt_n_pct(sum(.x == 1), n(), digits), .groups="drop") %>%
    pivot_wider(names_from = .g, values_from = value)
  
  pval <- fmt_p(p_cat(df2$.g, factor(df2$.x, levels = c(0,1))))
  
  tibble(
    Characteristic = label,
    Overall = overall,
    !!!bytab,
    `p value` = pval
  )
}
N <- nrow(discovery_data)

tab_add_S6 <- bind_rows(
  cat_block(discovery_data, by = "group", var = "fh_cvd",
            header = "Family history of CVD, n (%)", digits = 1),
  
  cat_block(discovery_data, by = "group", var = "ethnicity",
            header = "Ethnicity, n (%)", digits = 1),
  

  bin_row(discovery_data, by = "group", var = "Previous_drink",
          label = "Previous drinking, n (%)", digits = 1),
  
  bin_row(discovery_data, by = "group", var = "Current_drink",
          label = "Current drinking, n (%)", digits = 1)
  

  # bin_row(discovery_data, by="group", var="Unknown_drink", label="Unknown drinking, n (%)", digits=1)
)

tab_add_S6
write.csv(tab_add_S6, "S6_add_family_ethnicity_alcohol.csv", row.names = FALSE)

df = df_final[,c(c("eid_b", "MACE1_event",
                   "fh_cvd", "ethnicity", "alcohol_status", "Previous_drink", "Current_drink", 
                   "Unknown_drink"))]


head(df)

library(dplyr)
library(tidyr)

df2 <- df_final %>%
  mutate(
    MACE_grp = factor(MACE1_event, levels = c(0, 1), labels = c("No MACE", "MACE")),
    
    # ---- Ethnicity: 5 groups, NA -> White, Other+Unknown -> "Other or Unknown"
    eth_chr = as.character(ethnicity),
    ethnicity5 = case_when(
      is.na(eth_chr) ~ "White",
      eth_chr %in% c("Asian", "Asian or Asian British", "Chinese") ~ "Asian",
      eth_chr %in% c("Black", "Black or Black British") ~ "Black",
      eth_chr == "White" ~ "White",
      eth_chr == "Mixed" ~ "Mixed",
      TRUE ~ "Other or Unknown"
    ),
    ethnicity5 = factor(ethnicity5, levels = c("Asian","Black","White","Mixed","Other or Unknown")),
    
    # ---- Family history: missing -> Unknown
    fh_chr = as.character(fh_cvd),
    fh3 = case_when(
      is.na(fh_chr) | fh_chr == "" ~ "Unknown",
      fh_chr %in% c("Yes","No","Unknown") ~ fh_chr,
      TRUE ~ "Unknown"
    ),
    fh3 = factor(fh3, levels = c("Yes","No","Unknown"))
  ) %>%
  select(-eth_chr, -fh_chr)
fmt_n_pct <- function(n, denom, digits = 0){
  pct <- 100 * n / denom
  sprintf(paste0("%d (%.", digits, "f%%)"), n, pct)
}

fmt_p <- function(p){
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

p_cat <- function(g, v){
  tab <- table(g, v)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  if (all(dim(tab) == c(2,2)) && any(chi$expected < 5)) fisher.test(tab)$p.value else chi$p.value
}

cat_block <- function(df, by, var, header, digits = 0){
  df2 <- df %>% mutate(.g = .data[[by]], .v = .data[[var]])
  
  pval <- fmt_p(p_cat(df2$.g, df2$.v))
  
  overall <- df2 %>%
    count(.v, name="n") %>%
    mutate(Overall = fmt_n_pct(n, nrow(df2), digits)) %>%
    select(.v, Overall)
  
  bytab <- df2 %>%
    count(.g, .v, name="n") %>%
    group_by(.g) %>%
    mutate(value = fmt_n_pct(n, sum(n), digits)) %>%
    ungroup() %>%
    select(.g, .v, value) %>%
    pivot_wider(names_from = .g, values_from = value)
  
  out <- overall %>%
    left_join(bytab, by = ".v") %>%
    mutate(Characteristic = paste0("  ", .v),
           `p value` = "") %>%
    select(Characteristic, Overall, everything(), `p value`)
  
  header_row <- out[0, ]
  header_row[1, "Characteristic"] <- header
  header_row[1, "p value"] <- pval
  header_row[, setdiff(names(out), c("Characteristic","p value"))] <- ""
  
  bind_rows(header_row, out)
}
tab_eth <- cat_block(df2, by = "MACE_grp", var = "ethnicity5",
                     header = "Ethnicity, n (%)", digits = 0)

tab_fh  <- cat_block(df2, by = "MACE_grp", var = "fh3",
                     header = "Family history of CVD, n (%)", digits = 0)

tab_add <- bind_rows(tab_eth, tab_fh)

tab_add
write.csv(tab_add, "Table_add_ethnicity_famhx_byMACE.csv", row.names = FALSE)
library(dplyr)

df2 <- df2 %>%
  mutate(
    alc4 = case_when(
      is.na(alcohol_status) | as.character(alcohol_status) == "" ~ "Unknown",
      as.character(alcohol_status) %in% c("Never","Previous","Current","Unknown") ~ as.character(alcohol_status),
      TRUE ~ "Unknown"
    ),
    alc4 = factor(alc4, levels = c("Never","Previous","Current","Unknown")),
    
    Previous_drink2 = as.integer(alc4 == "Previous"),
    Current_drink2  = as.integer(alc4 == "Current"),
    Unknown_drink2  = as.integer(alc4 == "Unknown")
  )
library(dplyr)
library(tidyr)

bin_row <- function(df, by, var, label, digits = 0){
  x <- df[[var]]
  x[is.na(x)] <- 0
  df2 <- df %>% mutate(.g = .data[[by]], .x = x)
  
  overall <- fmt_n_pct(sum(df2$.x == 1), nrow(df2), digits)
  
  bytab <- df2 %>%
    group_by(.g) %>%
    summarise(value = fmt_n_pct(sum(.x == 1), n(), digits), .groups="drop") %>%
    pivot_wider(names_from = .g, values_from = value)
  
  pval <- fmt_p(p_cat(df2$.g, factor(df2$.x, levels = c(0,1))))
  
  tibble(
    Characteristic = label,
    Overall = overall,
    !!!bytab,
    `p value` = pval
  )
}

tab_drink <- bind_rows(
  bin_row(df2, by = "MACE_grp", var = "Previous_drink2", label = "Previous drinking, n (%)"),
  bin_row(df2, by = "MACE_grp", var = "Current_drink2",  label = "Current drinking, n (%)"),
  bin_row(df2, by = "MACE_grp", var = "Unknown_drink2",  label = "Unknown drinking, n (%)")
)

tab_drink
tab_add <- bind_rows(
  tab_eth,
  tab_fh,
  tab_drink
)

write.csv(tab_add, "Table_add_eth_fh_drink_byMACE.csv", row.names = FALSE)


colnames(score_plaque_cli)[1] = "eid_b"
score_plaque_cli = score_plaque_cli[,c(1,5)]

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
library(tidyr)

# -------------------------

# -------------------------
plaque_df <- score_plaque_cli %>%
  mutate(
    plaque_grp = factor(presence, levels = c(0, 1), labels = c("No plaque", "Plaque present")),
    

    eth_chr = as.character(ethnicity),
    ethnicity5 = case_when(
      is.na(eth_chr) ~ "White",
      eth_chr %in% c("Asian", "Asian or Asian British", "Chinese") ~ "Asian",
      eth_chr %in% c("Black", "Black or Black British") ~ "Black",
      eth_chr == "White" ~ "White",
      eth_chr == "Mixed" ~ "Mixed",
      TRUE ~ "Other or Unknown"
    ),
    ethnicity5 = factor(ethnicity5, levels = c("Asian","Black","White","Mixed","Other or Unknown")),
    
    # Family history: missing -> Unknown
    fh_chr = as.character(fh_cvd),
    fh3 = case_when(
      is.na(fh_chr) | fh_chr == "" ~ "Unknown",
      fh_chr %in% c("Yes","No","Unknown") ~ fh_chr,
      TRUE ~ "Unknown"
    ),
    fh3 = factor(fh3, levels = c("Yes","No","Unknown")),
    

    alc_chr = as.character(alcohol_status),
    alc4 = case_when(
      is.na(alc_chr) | alc_chr == "" ~ "Unknown",
      alc_chr %in% c("Never","Previous","Current","Unknown") ~ alc_chr,
      TRUE ~ "Unknown"
    ),
    alc4 = factor(alc4, levels = c("Never","Previous","Current","Unknown")),
    
    Previous_drink2 = as.integer(alc4 == "Previous"),
    Current_drink2  = as.integer(alc4 == "Current"),
    Unknown_drink2  = as.integer(alc4 == "Unknown")
  ) %>%
  select(-eth_chr, -fh_chr, -alc_chr)

# -------------------------

# -------------------------
fmt_n_pct <- function(n, denom, digits = 0){
  pct <- 100 * n / denom
  sprintf(paste0("%d (%.", digits, "f%%)"), n, pct)
}
fmt_p <- function(p){
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
p_cat_pearson <- function(g, v){
  suppressWarnings(chisq.test(table(g, v), correct = FALSE)$p.value)
}

cat_block <- function(df, by, var, header, digits = 0){
  df2 <- df %>% mutate(.g = .data[[by]], .v = .data[[var]])
  pval <- fmt_p(p_cat_pearson(df2$.g, df2$.v))
  
  overall <- df2 %>%
    count(.v, name="n") %>%
    mutate(Overall = fmt_n_pct(n, nrow(df2), digits)) %>%
    select(.v, Overall)
  
  bytab <- df2 %>%
    count(.g, .v, name="n") %>%
    group_by(.g) %>%
    mutate(value = fmt_n_pct(n, sum(n), digits)) %>%
    ungroup() %>%
    select(.g, .v, value) %>%
    pivot_wider(names_from = .g, values_from = value)
  
  out <- overall %>%
    left_join(bytab, by = ".v") %>%
    mutate(Characteristic = paste0("  ", .v),
           `p value` = "") %>%
    select(Characteristic, Overall, everything(), `p value`)
  
  header_row <- out[0, ]
  header_row[1, "Characteristic"] <- header
  header_row[1, "p value"] <- pval
  header_row[, setdiff(names(out), c("Characteristic","p value"))] <- ""
  
  bind_rows(header_row, out)
}

bin_row <- function(df, by, var, label, digits = 0){
  x <- df[[var]]
  x[is.na(x)] <- 0
  df2 <- df %>% mutate(.g = .data[[by]], .x = x)
  
  overall <- fmt_n_pct(sum(df2$.x == 1), nrow(df2), digits)
  
  bytab <- df2 %>%
    group_by(.g) %>%
    summarise(value = fmt_n_pct(sum(.x == 1), n(), digits), .groups="drop") %>%
    pivot_wider(names_from = .g, values_from = value)
  
  pval <- fmt_p(p_cat_pearson(df2$.g, factor(df2$.x, levels = c(0,1))))
  
  tibble(
    Characteristic = label,
    Overall = overall,
    !!!bytab,
    `p value` = pval
  )
}

# -------------------------

# -------------------------
tab_plaque_add <- bind_rows(
  cat_block(plaque_df, by = "plaque_grp", var = "ethnicity5", header = "Ethnicity, n (%)", digits = 0),
  cat_block(plaque_df, by = "plaque_grp", var = "fh3",       header = "Family history of CVD, n (%)", digits = 0),
  bin_row(plaque_df,  by = "plaque_grp", var = "Previous_drink2", label = "Previous drinking, n (%)", digits = 0),
  bin_row(plaque_df,  by = "plaque_grp", var = "Current_drink2",  label = "Current drinking, n (%)", digits = 0)

  # ,bin_row(plaque_df, by="plaque_grp", var="Unknown_drink2", label="Unknown drinking, n (%)", digits=0)
)

tab_plaque_add
write.csv(tab_plaque_add, "PlaqueTable_add_eth_fh_drink.csv", row.names = FALSE)
library(dplyr)
library(tidyr)

ukb_base2 <- ukb_base2 %>%
  mutate(
    # family history：NA -> Unknown
    fh3 = case_when(
      is.na(fh_cvd) | as.character(fh_cvd) == "" ~ "Unknown",
      TRUE ~ as.character(fh_cvd)
    ),
    fh3 = factor(fh3, levels = c("Yes","No","Unknown")),
    

    eth_chr = as.character(ethnicity),
    ethnicity5 = case_when(
      is.na(eth_chr) ~ "Other or Unknown",
      eth_chr %in% c("Asian", "Asian or Asian British", "Chinese") ~ "Asian",
      eth_chr %in% c("Black", "Black or Black British") ~ "Black",
      eth_chr == "White" ~ "White",
      eth_chr == "Mixed" ~ "Mixed",
      TRUE ~ "Other or Unknown"
    ),
    ethnicity5 = factor(ethnicity5, levels = c("Asian","Black","White","Mixed","Other or Unknown")),
    
    # alcohol：NA -> Unknown
    alc_chr = as.character(alcohol_status),
    alc4 = case_when(
      is.na(alc_chr) | alc_chr == "" ~ "Unknown",
      alc_chr %in% c("Never","Previous","Current","Unknown") ~ alc_chr,
      TRUE ~ "Unknown"
    ),
    alc4 = factor(alc4, levels = c("Never","Previous","Current","Unknown")),
    
    Previous_drink2 = as.integer(alc4 == "Previous"),
    Current_drink2  = as.integer(alc4 == "Current")
  ) %>%
  select(-eth_chr, -alc_chr)
fmt_n_pct <- function(n, denom, digits = 1){
  pct <- 100 * n / denom
  sprintf(paste0("%d (%.", digits, "f%%)"), n, pct)
}
fmt_p <- function(p){
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
p_pearson <- function(g, v){
  suppressWarnings(chisq.test(table(g, v), correct = FALSE)$p.value)
}

cat_block <- function(df, by, var, header, digits = 1){
  df2 <- df %>% mutate(.g = .data[[by]], .v = .data[[var]])
  pval <- fmt_p(p_pearson(df2$.g, df2$.v))
  
  overall <- df2 %>%
    count(.v, name="n") %>%
    mutate(Overall = fmt_n_pct(n, nrow(df2), digits)) %>%
    select(.v, Overall)
  
  bytab <- df2 %>%
    count(.g, .v, name="n") %>%
    group_by(.g) %>%
    mutate(value = fmt_n_pct(n, sum(n), digits)) %>%
    ungroup() %>%
    select(.g, .v, value) %>%
    pivot_wider(names_from = .g, values_from = value)
  
  out <- overall %>%
    left_join(bytab, by = ".v") %>%
    mutate(Characteristic = paste0("  ", .v),
           `p value` = "") %>%
    select(Characteristic, Overall, everything(), `p value`)
  
  header_row <- out[0, ]
  header_row[1, "Characteristic"] <- header
  header_row[1, "p value"] <- pval
  header_row[, setdiff(names(out), c("Characteristic","p value"))] <- ""
  
  bind_rows(header_row, out)
}

bin_row <- function(df, by, var, label, digits = 1){
  x <- df[[var]]; x[is.na(x)] <- 0
  df2 <- df %>% mutate(.g = .data[[by]], .x = x)
  
  overall <- fmt_n_pct(sum(df2$.x == 1), nrow(df2), digits)
  
  bytab <- df2 %>%
    group_by(.g) %>%
    summarise(value = fmt_n_pct(sum(.x == 1), n(), digits), .groups="drop") %>%
    pivot_wider(names_from = .g, values_from = value)
  
  pval <- fmt_p(p_pearson(df2$.g, factor(df2$.x, levels = c(0,1))))
  
  tibble(Characteristic = label, Overall = overall, !!!bytab, `p value` = pval)
}
tab_S14_add <- bind_rows(
  cat_block(ukb_base2, by = "group", var = "fh3",
            header = "Family history of CVD, n (%)", digits = 1),
  
  cat_block(ukb_base2, by = "group", var = "ethnicity5",
            header = "Ethnicity, n (%)", digits = 1),
  
  bin_row(ukb_base2, by = "group", var = "Previous_drink2",
          label = "Previous drinking, n (%)", digits = 1),
  
  bin_row(ukb_base2, by = "group", var = "Current_drink2",
          label = "Current drinking, n (%)", digits = 1)
)

tab_S14_add
write.csv(tab_S14_add, "S14_add_fh_eth_drink_single_vs_longitudinal.csv", row.names = FALSE)