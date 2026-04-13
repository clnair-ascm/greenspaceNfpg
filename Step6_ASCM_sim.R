# ============================================================
# FAST pooled ASCM / SCM-only for NDVI Improved vs Stable
# Treated:  label_300 == "improved"
# Controls: label_300 == "stable"
# Post-treatment starts in 2021
# Covariates: baseline hypertension, diabetes, smoking (pre 2021)
# Outcome: fbg_dto (distance-to-optimal fasting glucose) by default
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
})

if (!requireNamespace("augsynth", quietly = TRUE)) {
  install.packages("augsynth", repos = "https://cloud.r-project.org")
}
library(augsynth)

# ------------------------------------------------------------
# 1) Read panel with DTO/SEV variables
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

# Define time-invariant treatment indicator
panel <- panel %>%
  mutate(
    treat = if_else(label_300 == "improved", 1L, 0L)
  )

# ------------------------------------------------------------
# 2) Choose outcome
#    Recommended: outcome_var <- "fbg_dto"
#    Other options: "fbg", "bmi_dto", "map_dto", "pp_dto"
# ------------------------------------------------------------
outcome_var <- "fbg"   # <<< change here if needed

panel <- panel %>%
  filter(!is.na(.data[[outcome_var]]))

# ------------------------------------------------------------
# 3) Build baseline covariates (ever hypertensive/diabetic/smoker before 2021)
# ------------------------------------------------------------
baseline_cov <- panel %>%
  filter(year <= 2020) %>%
  group_by(person_id) %>%
  summarise(
    hyp_base = as.integer(any(hypertension == 1, na.rm = TRUE)),
    dm_base  = as.integer(any(diabetes    == 1, na.rm = TRUE)),
    smk_base = as.integer(any(smoking     == 1, na.rm = TRUE)),
    .groups  = "drop"
  )

panel <- panel %>%
  left_join(baseline_cov, by = "person_id") %>%
  mutate(
    hyp_base = if_else(is.na(hyp_base), 0L, hyp_base),
    dm_base  = if_else(is.na(dm_base),  0L, dm_base),
    smk_base = if_else(is.na(smk_base), 0L, smk_base)
  )

# ------------------------------------------------------------
# 4) Sample a small set of controls for speed
#    - keep all treated units
# ------------------------------------------------------------
set.seed(123)

treated_ids <- panel %>%
  filter(treat == 1) %>%
  distinct(person_id) %>%
  pull(person_id)

control_ids_all <- panel %>%
  filter(treat == 0) %>%
  distinct(person_id) %>%
  pull(person_id)

n_controls_target <- 1000L   # very small = fast; can increase later (e.g., 3000–5000)
control_ids <- sample(
  control_ids_all,
  size    = min(length(control_ids_all), n_controls_target),
  replace = FALSE
)

panel_fast <- panel %>%
  filter(person_id %in% c(treated_ids, control_ids))

cat("Fast sample: treated =",
    length(unique(panel_fast$person_id[panel_fast$treat == 1])),
    ", controls =",
    length(unique(panel_fast$person_id[panel_fast$treat == 0])), "\n")

# ------------------------------------------------------------
# 5) Ensure each unit has enough pre/post data (>=3 pre, >=2 post)
# ------------------------------------------------------------
panel_fast <- panel_fast %>%
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

unit_treat <- panel_fast %>%
  group_by(person_id) %>%
  summarise(treat_unit = max(treat), .groups = "drop")

cat("After coverage: treated =",
    sum(unit_treat$treat_unit == 1),
    ", controls =",
    sum(unit_treat$treat_unit == 0), "\n")

# ------------------------------------------------------------
# 6) Fit SCM-only (fast ASCM without Ridge augmentation)
# ------------------------------------------------------------
form <- as.formula(paste0(outcome_var, " ~ treat | hyp_base + dm_base + smk_base"))
t_int <- 2021

fit_fast <- augsynth(
  form,
  unit   = person_id,
  time   = year,
  t_int  = t_int,
  data   = panel_fast,
  progfunc = "None",   # no Ridge augmentation
  scm      = TRUE      # SCM weights only
)

# Extract ATT & SE as numeric scalars (no summary() needed)
att_hat <- unname(as.numeric(fit_fast$att)[1])
se_hat  <- unname(as.numeric(fit_fast$se)[1])

ci95 <- att_hat + c(-1.96, 1.96) * se_hat

cat("\n=== FAST SCM-only results for", outcome_var, "===\n")
cat("ATT:", att_hat, "\n")
cat("SE :", se_hat, "\n")
cat("95% CI:", paste(round(ci95, 3), collapse = " to "), "\n")

# ------------------------------------------------------------
# 7) Optional: plot treated vs synthetic mean trajectories
# ------------------------------------------------------------
try(plot(fit_fast), silent = TRUE)

# ------------------------------------------------------------
# 8) Save results to RDS
# ------------------------------------------------------------
ascm_output <- list(
  outcome    = outcome_var,
  att        = att_hat,
  se         = se_hat,
  ci95       = ci95,
  n_treated  = sum(unit_treat$treat_unit == 1),
  n_controls = sum(unit_treat$treat_unit == 0),
  model      = fit_fast,
  data_used  = panel_fast
)

saveRDS(
  ascm_output,
  file = paste0("C:/Users/Owner/Desktop/GreenSpace/data/ASCM_fast_SCM_", outcome_var, ".rds")
)

cat("Saved to: ",
    paste0("C:/Users/Owner/Desktop/GreenSpace/data/ASCM_fast_SCM_", outcome_var, ".rds"), "\n")
