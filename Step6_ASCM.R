# ============================================================
# Pooled ASCM (no strata) for NDVI Improved vs Stable
# Outcome: fbg (raw fasting glucose) or fbg_dto, etc.
# Treated: label_300 == "improved"
# Controls: label_300 == "stable"
# Post-treatment starts in 2021
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
})

# Install/load augsynth if needed
if (!requireNamespace("augsynth", quietly = TRUE)) {
  install.packages("augsynth", repos = "https://cloud.r-project.org")
}
library(augsynth)

# ------------------------------------------------------------
# 1) Read your panel with DTO/SEV variables
# ------------------------------------------------------------
panel <- readRDS("C:/Users/Owner/Desktop/GreenSpace/data/panel_16_24_R300_DTO_SEV.rds")

panel <- panel %>%
  mutate(
    person_id = as.character(person_id),
    year      = as.integer(year),
    label_300 = as.character(label_300)
  ) %>%
  filter(label_300 %in% c("improved", "stable"),
         year >= 2016, year <= 2024)

# Define time-invariant treatment (ever-treated) indicator
panel <- panel %>%
  mutate(
    treat = if_else(label_300 == "improved", 1L, 0L)
  )

# ------------------------------------------------------------
# 2) Choose outcome
#    - for raw fasting glucose: outcome_var <- "fbg"
#    - for distance-to-optimal: outcome_var <- "fbg_dto"
#    - you can later run again with "bmi_dto", "map_dto", "pp_dto"
# ------------------------------------------------------------
outcome_var <- "fbg"   # <-- change to "fbg_dto", "bmi_dto", "map_dto", "pp_dto" as needed

panel <- panel %>%
  filter(!is.na(.data[[outcome_var]]))

# ------------------------------------------------------------
# 3) Build baseline covariates (last pre-treatment year ≤ 2020)
# ------------------------------------------------------------
baseline <- panel %>%
  filter(year <= 2020) %>%
  arrange(person_id, year) %>%
  group_by(person_id) %>%
  slice_max(order_by = year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    person_id,
    age_baseline = as.numeric(age),
    sex_baseline = as.factor(sex)
  )

panel <- panel %>%
  left_join(baseline, by = "person_id")

# ------------------------------------------------------------
# 4) Common support trimming for donors (optional but recommended)
#    Trim controls whose pre-period mean outcome is far from treated mean
# ------------------------------------------------------------
trim_common_support <- function(df, outcome, band_sd = 2) {
  pre_df <- df %>% filter(year <= 2020)
  mu_i <- pre_df %>%
    group_by(person_id) %>%
    summarise(pre_mean = mean(.data[[outcome]], na.rm = TRUE),
              .groups = "drop")
  
  df2 <- df %>% left_join(mu_i, by = "person_id")
  
  mu_t <- mean(df2$pre_mean[df2$treat == 1], na.rm = TRUE)
  sd_t <- sd(df2$pre_mean[df2$treat == 1],   na.rm = TRUE)
  if (!is.finite(sd_t) || sd_t == 0) sd_t <- 1e-6
  
  df2 %>%
    filter(
      treat == 1 |
        (treat == 0 & is.finite(pre_mean) &
           abs(pre_mean - mu_t) <= band_sd * sd_t)
    ) %>%
    select(-pre_mean)
}

# Apply trimming
panel_trim <- trim_common_support(panel, outcome = outcome_var, band_sd = 2)

# ------------------------------------------------------------
# 5) Ensure sufficient pre/post coverage for each unit
#    (≥3 pre years, ≥2 post years with outcome)
# ------------------------------------------------------------
panel_trim <- panel_trim %>%
  group_by(person_id) %>%
  mutate(
    pre  = year <= 2020,
    post = year >= 2021
  ) %>%
  filter(
    sum(pre  & !is.na(.data[[outcome_var]])) >= 3,
    sum(post & !is.na(.data[[outcome_var]])) >= 2
  ) %>%
  ungroup()

# Quick check: treated & controls
unit_treat <- panel_trim %>%
  group_by(person_id) %>%
  summarise(treat_unit = max(treat), .groups = "drop")

cat("Units after trimming & coverage filters:\n")
cat("  Treated units: ", sum(unit_treat$treat_unit == 1), "\n")
cat("  Control units: ", sum(unit_treat$treat_unit == 0), "\n")

# ------------------------------------------------------------
# 6) Fit ASCM (pooled)
# ------------------------------------------------------------
# Formula: outcome ~ treat | covariates
use_covariates <- TRUE

if (use_covariates && all(c("age_baseline","sex_baseline") %in% names(panel_trim))) {
  form <- as.formula(paste0(outcome_var, " ~ treat | age_baseline + sex_baseline"))
} else {
  form <- as.formula(paste0(outcome_var, " ~ treat"))
}

t_int <- 2021  # first post-treatment year

fit_ascm <- augsynth(
  form,
  unit   = person_id,
  time   = year,
  t_int  = t_int,
  data   = panel_trim,
  progfunc = "Ridge",
  scm      = TRUE
)

sm <- summary(fit_ascm)
print(sm)

# Extract ATT & SE robustly
att_hat <- if (!is.null(sm$att)) sm$att else if (!is.null(sm$Estimate)) sm$Estimate else fit_ascm$att
se_hat  <- if (!is.null(sm$se))  sm$se  else if (!is.null(sm$Std..Error)) sm$Std..Error else fit_ascm$se

ci95 <- att_hat + c(-1.96, 1.96)*se_hat

cat("\n=== ASCM pooled results for", outcome_var, "===\n")
cat("ATT:", att_hat, "\n")
cat("SE :", se_hat, "\n")
cat("95% CI:", paste(round(ci95, 3), collapse = " to "), "\n")

# ------------------------------------------------------------
# 7) (Optional) Plot treated vs synthetic mean trajectories
# ------------------------------------------------------------
# This plot does NOT draw all unit lines; it's just group means.
try({
  plot(fit_ascm)  # if your augsynth version supports plotting
}, silent = TRUE)
