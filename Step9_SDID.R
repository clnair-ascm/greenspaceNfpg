# ============================================================
# Incident diabetes in high-risk subgroup (baseline FBG >= 5.6)
# Clean discrete-time hazard DID
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(sandwich)
  library(lmtest)
  library(tibble)
})

# ------------------------------------------------------------
# 0) Settings
# ------------------------------------------------------------
fp <- "/Users/yuqing/Downloads/GreenSpace/data/panel_16_24_R300_DTO_SEV.rds"
intervention_year <- 2021

# ------------------------------------------------------------
# 1) Read panel
# ------------------------------------------------------------
panel <- readRDS(fp)

if ("...1" %in% names(panel)) {
  panel <- panel %>% select(-`...1`)
}

panel <- panel %>%
  mutate(
    person_id = as.character(person_id),
    year = as.integer(year),
    label_300 = as.character(label_300)
  ) %>%
  distinct(person_id, year, .keep_all = TRUE)

# 如果你的数据没有 2023，这里会自动按实际年份处理
years_use <- sort(unique(panel$year))

# ------------------------------------------------------------
# 2) Keep improved / stable only + complete annual diabetes records
# ------------------------------------------------------------
panel2 <- panel %>%
  filter(
    label_300 %in% c("improved", "stable"),
    year %in% years_use
  ) %>%
  group_by(person_id) %>%
  filter(
    n_distinct(year) == length(years_use),
    all(!is.na(diabetes))
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 3) Baseline table (2016), rename clearly to avoid .x/.y
# ------------------------------------------------------------
base <- panel2 %>%
  filter(year == 2016) %>%
  transmute(
    person_id,
    label_base = label_300,
    improved_base = as.integer(label_300 == "improved"),
    diabetes_base = diabetes,
    age_base = age,
    sex_base = sex,
    bmi_base = bmi,
    smoking_base = smoking,
    hypertension_base = hypertension,
    fbg_base = fbg
  ) %>%
  filter(
    diabetes_base == 0,
    !is.na(fbg_base),
    fbg_base >= 5.6,
    complete.cases(age_base, sex_base, bmi_base, smoking_base, hypertension_base)
  )

cat("High-risk baseline n =", nrow(base), "\n")

# 只保留高风险基线人群
df <- panel2 %>%
  semi_join(base, by = "person_id")

# ------------------------------------------------------------
# 4) First diabetes year
# ------------------------------------------------------------
first_event <- df %>%
  group_by(person_id) %>%
  summarise(
    first_diab_year = if (any(diabetes == 1, na.rm = TRUE)) {
      min(year[diabetes == 1], na.rm = TRUE)
    } else {
      Inf
    },
    .groups = "drop"
  )

# ------------------------------------------------------------
# 5) Build person-year risk set
# ------------------------------------------------------------
risk <- df %>%
  left_join(base, by = "person_id") %>%
  left_join(first_event, by = "person_id") %>%
  filter(year > 2016) %>%   # follow-up starts after baseline
  mutate(
    event = as.integer(is.finite(first_diab_year) & year == first_diab_year),
    in_risk = year <= first_diab_year,
    post = as.integer(year >= intervention_year),
    treat_time = improved_base * post,
    year_f = relevel(factor(year), ref = "2020")
  ) %>%
  filter(in_risk)

cat("Person-years in risk set =", nrow(risk), "\n")
cat("Total incident diabetes events =", sum(risk$event), "\n")

# ------------------------------------------------------------
# 6) Main model: discrete-time hazard DID
# ------------------------------------------------------------
model <- glm(
  event ~ improved_base + treat_time + year_f +
    age_base + sex_base + bmi_base + smoking_base + hypertension_base + fbg_base,
  family = binomial(),
  data = risk
)

# Cluster-robust SE by person
vc <- vcovCL(model, cluster = risk$person_id)
robust_res <- coeftest(model, vcov. = vc)
print(robust_res)

# ------------------------------------------------------------
# 7) Robust OR table
# ------------------------------------------------------------
beta <- coef(model)
se <- sqrt(diag(vc))

or_tab <- tibble(
  term = names(beta),
  OR = exp(beta),
  conf.low = exp(beta - 1.96 * se),
  conf.high = exp(beta + 1.96 * se),
  p.value = 2 * pnorm(-abs(beta / se))
)

print(or_tab)

# ------------------------------------------------------------
# 8) Simple cumulative incidence summary
# ------------------------------------------------------------
person_event <- risk %>%
  group_by(person_id, improved_base) %>%
  summarise(
    incident = as.integer(any(event == 1)),
    .groups = "drop"
  )

inc_summary <- person_event %>%
  group_by(improved_base) %>%
  summarise(
    n = n(),
    cases = sum(incident),
    risk = cases / n,
    .groups = "drop"
  ) %>%
  mutate(group = if_else(improved_base == 1, "Improved", "Stable")) %>%
  select(group, n, cases, risk)

print(inc_summary)

# ------------------------------------------------------------
# 9) Forest plot: ONLY key coefficients
# ------------------------------------------------------------
plot_df <- or_tab %>%
  filter(term %in% c(
    "improved_base", "treat_time",
    "age_base", "sex_base", "bmi_base",
    "smoking_base", "hypertension_base", "fbg_base"
  )) %>%
  mutate(
    label = recode(
      term,
      improved_base = "Improved vs Stable\n(baseline difference)",
      treat_time = "Greening effect\n(Improved × Post2021)",
      age_base = "Age (baseline)",
      sex_base = "Sex (baseline)",
      bmi_base = "BMI (baseline)",
      smoking_base = "Smoking (baseline)",
      hypertension_base = "Hypertension (baseline)",
      fbg_base = "FBG (baseline)"
    ),
    label = factor(label, levels = rev(label))
  )

p <- ggplot(plot_df, aes(x = label, y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_y_log10() +
  theme_classic(base_size = 12) +
  labs(
    x = "",
    y = "Odds ratio (log scale, cluster-robust 95% CI)",
    title = "Incident diabetes in high-risk subgroup (baseline FBG ≥ 5.6)"
  )

print(p)