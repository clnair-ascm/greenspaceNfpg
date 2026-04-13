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
  "/Users/yqdai/Downloads/GreenSpace/data/panel_16_24_R300_DTO_SEV.rds"
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

combined_plot <- wrap_plots(p_list, ncol = 2)
# Option B: High-res TIFF (often preferred by journals)
tiff(
  filename = "/Users/yqdai/Downloads/GreenSpace/sdid_fbg_strata_fbg.tiff",
  width    = 10,
  height   = 16,
  units    = "in",
  res      = 300
)
print(combined_plot)
dev.off()


library(dplyr)
library(purrr)
library(synthdid)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 1) Run SDID by stratum and store ATT + SE
#    (这里用你已经准备好的 panel_fbg 和 run_sdid_stratum)
# ------------------------------------------------------------

# 想分层的变量；现在只做 sex，如果之后要扩展可以再加
strata_vars <- "sex"
# strata_vars <- c("sex", "hypertension", "diabetes", "smoking")

# 每一行对应一个 (var, level) 组合
strata_specs <- expand.grid(
  var   = strata_vars,
  level = 0:1,
  stringsAsFactors = FALSE
)

# 保存 est 对象（如果以后要画轨迹图）
sdid_strata <- list()

# 一次性把 ATT 和 SE 存下来
effects_df <- purrr::pmap_dfr(
  strata_specs,
  function(var, level) {
    est  <- run_sdid_stratum(panel_fbg, strat_var = var, level = level)
    name <- paste0(var, "_", level)
    sdid_strata[[name]] <<- est   # 存到外层 list
    
    if (is.null(est)) {
      return(tibble(
        stratum = name,
        att     = NA_real_,
        se      = NA_real_
      ))
    }
    
    tibble(
      stratum = name,
      att     = as.numeric(est),
      # 关键改动：用 HC0（解析型、非重抽样）而不是 placebo，加快很多
      se      = sqrt(as.numeric(vcov(est, method = "HC0")))
    )
  }
)

# ------------------------------------------------------------
# 2) 整理 effects_df
# ------------------------------------------------------------

effects_df <- effects_df %>%
  filter(!is.na(att)) %>%
  mutate(
    ci_low  = att - 1.96 * se,
    ci_high = att + 1.96 * se
  ) %>%
  # 把 "sex_0" 拆成 var = sex, level = 0/1
  tidyr::separate(stratum, into = c("var", "level"), sep = "_", remove = FALSE) %>%
  mutate(
    level = ifelse(level == "1", "Yes", "No"),
    var = factor(
      var,
      levels = c("sex", "hypertension", "diabetes", "smoking"),
      labels = c("Male (sex = 1)", "Hypertension", "Diabetes", "Current smoking")
    ),
    level = factor(level, levels = c("No", "Yes"))
  )

# ------------------------------------------------------------
# 3) PNAS 风格的 forest plot
# ------------------------------------------------------------

gg_sdid_strata <- ggplot(effects_df,
                         aes(x = att, y = level)) +
  # 零效应线
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour   = "grey60") +
  # 95% CI 横线
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height    = 0.2,
                 linewidth = 0.8) +
  # 点估计
  geom_point(size = 2.8) +
  # 分面：每个分层变量一个 panel（现在只有 sex）
  facet_wrap(~ var, ncol = 2, scales = "free_y") +
  labs(
    title    = "Stratum-specific effects of greenness improvement on fasting glucose",
    subtitle = "Synthetic difference-in-differences estimates (ATT, NDVI 300 m)",
    x        = "ATT on fasting glucose (mmol/L)",
    y        = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text         = element_text(size = 12, colour = "grey20"),
    axis.title.x      = element_text(size = 14),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 13, face = "bold"),
    plot.title        = element_text(size = 18, face = "bold"),
    plot.subtitle     = element_text(size = 13),
    plot.margin       = margin(t = 5, r = 10, b = 5, l = 10)
  )

gg_sdid_strata

# ------------------------------------------------------------
# 4) 高分辨率输出（PNAS style）
# ------------------------------------------------------------

ggsave(
  filename    = "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/fig/FBG_sdid_strata_forest.tiff",
  plot        = gg_sdid_strata,
  width       = 7.5,
  height      = 5.5,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)
