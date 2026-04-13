# ============================================================
# SDID for FBG: Improved vs Stable (NDVI 300m, 2016–2024)
# Stratified by sex / hypertension / diabetes / smoking
# ============================================================

rm(list = ls())

# ---- packages ----
library(dplyr)
library(readr)
library(synthdid)
library(rlang)   # 用于 .data[[strat_var]]

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

# ------------------------------------------------------------
# 2) Restrict to Improved vs Stable, FBG available, year window
# ------------------------------------------------------------
panel_fbg <- panel %>%
  filter(
    label_300 %in% c("improved", "stable"),
    year >= 2016, year <= 2024,
    !is.na(fbg)              # ← 用 fbg
  ) %>%
  mutate(
    pre  = year <= 2020,
    post = year >= 2021
  ) %>%
  group_by(person_id) %>%
  # 要求：每个人至少 3 个 pre + 2 个 post 有 FBG
  filter(
    sum(pre  & !is.na(fbg)) >= 3,
    sum(post & !is.na(fbg)) >= 2
  ) %>%
  ungroup() %>%
  arrange(person_id, year)

# ------------------------------------------------------------
# 3) Define time-varying treatment indicator
#    treat_time = 1 for Improved units in post years (>= 2021)
# ------------------------------------------------------------
panel_fbg <- panel_fbg %>%
  mutate(
    treat_time = if_else(label_300 == "improved" & year >= 2021, 1L, 0L)
  ) %>%
  arrange(person_id, year)

# ------------------------------------------------------------
# 4) 简单检查：整体 SDID 是否能跑（可选）
# ------------------------------------------------------------
pm_fbg <- panel.matrices(
  panel     = as.data.frame(panel_fbg),
  unit      = "person_id",
  time      = "year",
  outcome   = "fbg",       # ← outcome 用 fbg
  treatment = "treat_time"
)

Y  <- pm_fbg$Y
N0 <- pm_fbg$N0
T0 <- pm_fbg$T0

cat("Overall Y dim:", dim(Y), "\n")
cat("Overall N0 (controls):", N0, "\n")
cat("Overall T0 (pre-periods):", T0, "\n")

sdid_fbg_overall <- synthdid_estimate(Y, N0, T0)
# plot(sdid_fbg_overall)

# ------------------------------------------------------------
# 5) 封装：在某个分层 + 取值下，构造“平衡子面板”并做 SDID
# ------------------------------------------------------------
run_sdid_stratum <- function(panel_dat, strat_var, level) {
  
  # ---- 5.1 先按分层变量筛选 ----
  subdat <- panel_dat %>%
    filter(
      !is.na(.data[[strat_var]]),
      .data[[strat_var]] == level
    )
  
  if (nrow(subdat) == 0) {
    warning(paste0("No data for ", strat_var, " = ", level))
    return(NULL)
  }
  
  # ---- 5.2 保证每个 unit-year 只有一行 ----
  subdat <- subdat %>%
    group_by(person_id, year) %>%
    slice(1) %>%
    ungroup()
  
  # ---- 5.3 在这个分层内部构造“平衡面板” ----
  all_years <- sort(unique(subdat$year))
  
  subdat_bal <- subdat %>%
    group_by(person_id) %>%
    # 要求：这个人在该分层的所有年份 all_years 都有记录
    filter(n_distinct(year) == length(all_years)) %>%
    ungroup()
  
  if (nrow(subdat_bal) == 0) {
    warning(paste0("No balanced units for ", strat_var, " = ", level))
    return(NULL)
  }
  
  # ---- 5.4 检查 treated / control 是否都存在 ----
  treated_units  <- unique(subdat_bal$person_id[subdat_bal$treat_time == 1])
  control_units  <- unique(subdat_bal$person_id[subdat_bal$treat_time == 0])
  
  if (length(treated_units) == 0L || length(control_units) == 0L) {
    warning(paste0("Skipping ", strat_var, " = ", level,
                   ": no treated or no controls in balanced subset."))
    return(NULL)
  }
  
  # ---- 5.5 丢给 panel.matrices & synthdid ----
  pm <- panel.matrices(
    panel     = as.data.frame(subdat_bal),
    unit      = "person_id",
    time      = "year",
    outcome   = "fbg",     # ← 这里同样用 fbg
    treatment = "treat_time"
  )
  
  est <- synthdid_estimate(pm$Y, pm$N0, pm$T0)
  est
}

# ------------------------------------------------------------
# 6) 按 sex / hypertension / diabetes / smoking 分层，画 8 张图
# ------------------------------------------------------------
library(ggplot2)

strata_vars  <- c("sex", "hypertension", "diabetes", "smoking")
sdid_strata  <- list()
p_list <- list()   # ← 用来存 8 张图

for (v in strata_vars) {
  for (lvl in 0:1) {
    
    est  <- run_sdid_stratum(panel_fbg, v, lvl)
    name <- paste0(v, "_", lvl)
    sdid_strata[[name]] <- est
    
    if (!is.null(est)) {
      p <- synthdid::synthdid_plot(est)
      p <- p + ggplot2::ggtitle(paste0("SDID: FBG\n", v, " = ", lvl))
    } else {
      # 如果没法做 SDID，生成一个空白图
      p <- ggplot2::ggplot() + 
        ggplot2::theme_void() +
        ggplot2::ggtitle(paste0("No valid SDID\n", v, " = ", lvl))
    }
    
    p_list[[name]] <- p
  }
}

library(patchwork)

wrap_plots(p_list, ncol = 2)
