rm(list = ls())

library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# Suppose you have a list of yearly data frames as before:
panel <- readRDS("/Users/yqdai/Downloads/GreenSpace/data/D_all_SI_greenT/panel_16_24_R1000.rds")

panel_filtered_ids <- panel %>% 
  mutate(
    year = as.integer(year),
    age  = as.numeric(age)
  ) %>%
  group_by(person_id) %>%
  # 标记：这个人是否在 2016 年及以后年龄 > 60
  mutate(eligible = any(year >= 2016 & age > 60, na.rm = TRUE)) %>%
  ungroup() %>%
  # 去掉不符合的 person_id
  filter(eligible) %>%
  select(-eligible)

expected_years <- setdiff(2016:2024, 2023)
expected_years
# [1] 2016 2017 2018 2019 2020 2021 2022 2024

library(dplyr)

fbg_continuity_r300 <- panel_filtered_ids %>%
  mutate(
    year = as.integer(year)
  ) %>%
  # 只看我们关心的年份（2016–2024，不含 2023）
  filter(year %in% expected_years) %>%
  group_by(person_id) %>%
  summarise(
    n_years_with_fbg = n_distinct(year[!is.na(fbg)]),
    years_have_fbg   = paste(sort(unique(year[!is.na(fbg)])), collapse = ","),
    years_missing_fbg = paste(setdiff(expected_years, sort(unique(year[!is.na(fbg)]))), collapse = ","),
    is_continuous    = n_years_with_fbg == length(expected_years),
    .groups = "drop"
  )

# 快速检查：是不是所有人都连续？
all(fbg_continuity_r300$is_continuous)
# TRUE 表示所有 person_id 都有完整 fbg（2016–2022 + 2024）
# FALSE 表示至少有一个缺
saveRDS(panel_filtered_ids, "/Users/yqdai/Downloads/GreenSpace/data/D_all_SI_greenT/panel_16_24_R1000.rds")