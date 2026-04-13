rm(list = ls())

# ---- packages ----
library(dplyr)
library(readr)
library(synthdid)
library(purrr)
library(tidyr)

# ------------------------------------------------------------
# 1) Read panel data
# ------------------------------------------------------------
panel <- readRDS(
  "C:/Users/Owner/Desktop/GreenSpace/data/panel_16_24_R300_DTO_SEV.rds"
)

# Drop the index column if present (...1)
if ("...1" %in% names(panel)) {
  panel <- panel %>% select(-`...1`)
}

# Basic types
panel <- panel %>%
  mutate(
    person_id = as.character(person_id),
    year      = as.integer(year),
    label_300 = as.character(label_300)
  )

# Choose outcome explicitly: BMI_DTO (change to "bmi_sev" if needed)
outcome_var <- "fbg_dto"

# ------------------------------------------------------------
# 2) Restrict to Improved vs Stable, outcome available, year window
# ------------------------------------------------------------
panel_bmi <- panel %>%
  filter(
    label_300 %in% c("improved", "stable"),
    year >= 2016, year <= 2024,
    !is.na(.data[[outcome_var]])    # ensures no NA for outcome
  ) %>%
  mutate(
    pre  = year <= 2020,
    post = year >= 2021
  ) %>%
  group_by(person_id) %>%
  # require at least 3 pre and 2 post outcome observations
  filter(
    sum(pre  & !is.na(.data[[outcome_var]])) >= 3,
    sum(post & !is.na(.data[[outcome_var]])) >= 2
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 3) Define time-varying treatment indicator
#    treat_time = 1 for Improved units in post years (>= 2021)
# ------------------------------------------------------------
panel_bmi <- panel_bmi %>%
  mutate(
    treat_time = if_else(label_300 == "improved" & year >= 2021, 1L, 0L)
  ) %>%
  arrange(person_id, year)

# Sanity checks
table(panel_bmi$label_300, useNA = "ifany")
range(panel_bmi$year, na.rm = TRUE)

# ------------------------------------------------------------
# 4) Baseline age/sex (use latest pre-treatment year <= 2020)
# ------------------------------------------------------------
baseline <- panel_bmi %>%
  filter(year <= 2020) %>%
  arrange(person_id, year) %>%
  group_by(person_id) %>%
  # last pre-treatment observation
  slice_max(order_by = year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    person_id,
    age_baseline = age,
    sex_baseline = sex
  )

panel_bmi <- panel_bmi %>%
  left_join(baseline, by = "person_id")

# ------------------------------------------------------------
# 5) Define age/sex strata
#    (adjust age cutpoints and sex coding as you like)
# ------------------------------------------------------------

panel_bmi <- panel_bmi %>%
  mutate(
    # Here we just treat sex_baseline as a factor "0","1",...
    # Replace labels with "male"/"female" if you know the coding.
    sex_strata = factor(sex_baseline),
    
    # Age groups: <40, 40-59, 60+
    age_group = cut(
      age_baseline,
      breaks = c(-Inf, 40, 60, Inf),
      labels = c("<40", "40-59", "60+"),
      right  = FALSE
    ),
    
    # Combined stratum: e.g. "0.<40", "1.60+", ...
    stratum = interaction(sex_strata, age_group, drop = TRUE)
  )

# Look at sample sizes by stratum and treatment group
table(panel_bmi$stratum, panel_bmi$label_300, useNA = "ifany")

# ------------------------------------------------------------
# 6) Helper: run SDID within a single stratum
# ------------------------------------------------------------

run_sdid_stratum <- function(df, outcome_var = "fpg_dto") {
  # df: data for a single stratum (all years, all persons in that stratum)
  
  stratum_level <- df$stratum[1]
  
  # Ensure outcome not missing
  df <- df %>% filter(!is.na(.data[[outcome_var]]))
  
  # Check treated vs control units in this stratum
  treat_units <- df %>%
    group_by(person_id) %>%
    summarise(any_treat = any(treat_time == 1L), .groups = "drop")
  
  n_treated_units <- sum(treat_units$any_treat)
  n_control_units <- nrow(treat_units) - n_treated_units
  
  if (n_treated_units == 0 || n_control_units == 0) {
    message("Skipping stratum ", stratum_level,
            ": needs both treated and control units.")
    return(NULL)
  }
  
  # Build panel matrices for this stratum
  pm <- panel.matrices(
    panel     = as.data.frame(df),
    unit      = "person_id",
    time      = "year",
    outcome   = outcome_var,
    treatment = "treat_time"
  )
  
  Y  <- pm$Y
  N0 <- pm$N0
  T0 <- pm$T0
  
  # Fit SDID
  est <- synthdid_estimate(Y, N0, T0)
  
  # Standard error (jackknife or placebo; choose what you prefer)
  se <- sqrt(vcov(est, method = "jackknife"))
  
  # Weight for aggregation: number of treated post observations in this stratum
  w <- df %>%
    summarise(w = sum(treat_time == 1L & year >= 2021 &
                        !is.na(.data[[outcome_var]]))) %>%
    pull(w)
  
  tibble(
    stratum          = stratum_level,
    n_units          = nrow(Y),
    n_control_units  = N0,
    n_treated_units  = nrow(Y) - N0,
    weight_raw       = w,
    att_hat          = as.numeric(est),          # SDID point estimate
    se_hat           = as.numeric(se)            # its SE
  )
}

# ------------------------------------------------------------
# 7) Run SDID in each age/sex stratum
# ------------------------------------------------------------

sdid_strata <- panel_bmi %>%
  filter(!is.na(stratum)) %>%
  group_by(stratum) %>%
  group_modify(~{
    res <- run_sdid_stratum(.x, outcome_var = outcome_var)
    if (is.null(res)) tibble() else res
  }) %>%
  ungroup()

sdid_strata

# ------------------------------------------------------------
# 8) Aggregate stratum-specific ATTs to an overall ATT
#    (weighted by # of treated post observations)
# ------------------------------------------------------------

if (nrow(sdid_strata) > 0) {
  
  sdid_strata <- sdid_strata %>%
    mutate(
      weight_norm = weight_raw / sum(weight_raw),
      contrib     = weight_norm * att_hat
    )
  
  overall_att <- sum(sdid_strata$contrib)
  
  # Approximate SE for the weighted average assuming independence across strata
  overall_var <- sum((sdid_strata$weight_norm^2) * (sdid_strata$se_hat^2))
  overall_se  <- sqrt(overall_var)
  
  overall_ci  <- overall_att + c(-1.96, 1.96) * overall_se
  
  list(
    strata_results = sdid_strata,
    overall = list(
      att = overall_att,
      se  = overall_se,
      ci95 = overall_ci
    )
  ) -> sdid_strata_summary
  
  sdid_strata_summary$overall
}

library(dplyr)
library(ggplot2)

# Prepare data for plotting
plot_df <- sdid_strata %>%
  mutate(
    ci_low  = att_hat - 1.96 * se_hat,
    ci_high = att_hat + 1.96 * se_hat,
    stratum = as.factor(stratum)
  )

# If you computed an overall ATT:
overall_att <- sdid_strata_summary$overall$att

ggplot(plot_df, aes(x = att_hat, y = stratum)) +
  # zero line
  geom_vline(xintercept = 0, linetype = "dashed") +
  # overall ATT (optional)
  geom_vline(xintercept = overall_att,
             linetype = "dotted") +
  # CI
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.2) +
  # point size = normalized weight
  geom_point(aes(size = weight_norm)) +
  scale_size_continuous(name = "Weight (treated post share)") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Stratum-specific SDID effects on BMI_DTO",
    subtitle = "Age × sex strata; error bars: 95% CI",
    x = "ATT on BMI_DTO (treated − synthetic)",
    y = "Age × sex stratum"
  )
