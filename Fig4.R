# ============================================================
# SDID for FBG: Improved vs Stable (NDVI 300m, 2016–2024)
# Stratified by sex / hypertension / diabetes / smoking
# Add overall estimate into diabetes panel for comparison
# ============================================================

rm(list = ls())

# ---- packages ----
library(dplyr)
library(readr)
library(synthdid)
library(rlang)      # for .data[[strat_var]]
library(ggplot2)
library(patchwork)
library(tibble)

# ------------------------------------------------------------
# 1) Read panel data
# ------------------------------------------------------------
panel <- readRDS(
  "/Users/yuqing/Downloads/GreenSpace/data/panel_16_24_R1000_DTO_SEV.rds"
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
    !is.na(fbg)
  ) %>%
  mutate(
    pre  = year <= 2020,
    post = year >= 2021
  ) %>%
  group_by(person_id) %>%
  # Require at least 3 pre + 2 post FBG observations per person
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
# 4) Overall SDID
# ------------------------------------------------------------
pm_fbg <- panel.matrices(
  panel     = as.data.frame(panel_fbg),
  unit      = "person_id",
  time      = "year",
  outcome   = "fbg",
  treatment = "treat_time"
)

Y  <- pm_fbg$Y
N0 <- pm_fbg$N0
T0 <- pm_fbg$T0

cat("Overall Y dim:", dim(Y), "\n")
cat("Overall N0 (controls):", N0, "\n")
cat("Overall T0 (pre-periods):", T0, "\n")

sdid_fbg_overall <- synthdid_estimate(Y, N0, T0)

# Optional:
# plot(sdid_fbg_overall)

# ------------------------------------------------------------
# 5) Helper function:
#    In a given stratum, build balanced sub-panel and run SDID
# ------------------------------------------------------------
run_sdid_stratum <- function(panel_dat, strat_var, level) {
  
  # ---- 5.1 Filter by stratum ----
  subdat <- panel_dat %>%
    filter(
      !is.na(.data[[strat_var]]),
      .data[[strat_var]] == level
    )
  
  if (nrow(subdat) == 0) {
    warning(paste0("No data for ", strat_var, " = ", level))
    return(NULL)
  }
  
  # ---- 5.2 Ensure one row per unit-year ----
  subdat <- subdat %>%
    group_by(person_id, year) %>%
    slice(1) %>%
    ungroup()
  
  # ---- 5.3 Build balanced panel within this stratum ----
  all_years <- sort(unique(subdat$year))
  
  subdat_bal <- subdat %>%
    group_by(person_id) %>%
    filter(n_distinct(year) == length(all_years)) %>%
    ungroup()
  
  if (nrow(subdat_bal) == 0) {
    warning(paste0("No balanced units for ", strat_var, " = ", level))
    return(NULL)
  }
  
  # ---- 5.4 Check both treated and controls exist ----
  treated_units <- unique(subdat_bal$person_id[subdat_bal$treat_time == 1])
  
  # Better to define controls as units never treated in the balanced subset
  unit_treat_status <- subdat_bal %>%
    group_by(person_id) %>%
    summarise(ever_treated = max(treat_time, na.rm = TRUE), .groups = "drop")
  
  control_units <- unit_treat_status$person_id[unit_treat_status$ever_treated == 0]
  
  if (length(treated_units) == 0L || length(control_units) == 0L) {
    warning(paste0(
      "Skipping ", strat_var, " = ", level,
      ": no treated or no controls in balanced subset."
    ))
    return(NULL)
  }
  
  # ---- 5.5 Run panel.matrices & synthdid ----
  pm <- panel.matrices(
    panel     = as.data.frame(subdat_bal),
    unit      = "person_id",
    time      = "year",
    outcome   = "fbg",
    treatment = "treat_time"
  )
  
  est <- synthdid_estimate(pm$Y, pm$N0, pm$T0)
  return(est)
}

# ------------------------------------------------------------
# 6) Stratified SDID plots: sex / hypertension / diabetes / smoking
# ------------------------------------------------------------
strata_vars <- c("sex", "hypertension", "diabetes", "smoking")
sdid_strata <- list()
p_list <- list()

for (v in strata_vars) {
  for (lvl in 0:1) {
    
    est  <- run_sdid_stratum(panel_fbg, v, lvl)
    name <- paste0(v, "_", lvl)
    sdid_strata[[name]] <- est
    
    if (!is.null(est)) {
      p <- synthdid::synthdid_plot(est)
      p <- p + ggplot2::ggtitle(paste0("SDID: FBG\n", v, " = ", lvl))
    } else {
      p <- ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle(paste0("No valid SDID\n", v, " = ", lvl))
    }
    
    p_list[[name]] <- p
  }
}

combined_plot <- wrap_plots(p_list, ncol = 2)

# Optional save:
# ggsave(
#   filename = "/Users/yuqing/Downloads/GreenSpace/sdid_fbg_strata_combined.pdf",
#   plot     = combined_plot,
#   width    = 10,
#   height   = 16,
#   units    = "in"
# )

# Optional TIFF:
# tiff(
#   filename = "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/fig/sdid_fbg_strata_fbg.tiff",
#   width    = 10,
#   height   = 16,
#   units    = "in",
#   res      = 300
# )
# print(combined_plot)
# dev.off()

# ------------------------------------------------------------
# 7) Build summary table (ATE + 95% CI)
#    Re-use existing estimates; do NOT re-estimate
# ------------------------------------------------------------
alpha     <- 0.05
z_crit    <- qnorm(1 - alpha / 2)
se_method <- "jackknife"

res_list <- list()

# ---- 7.1 Overall row: put into diabetes panel for comparison ----
tau_overall <- as.numeric(sdid_fbg_overall)
vc_overall  <- stats::vcov(sdid_fbg_overall, method = se_method)
se_overall  <- sqrt(as.numeric(vc_overall))

if (is.na(se_overall)) {
  warning("Overall SE is NA - overall CI will be NA")
}

overall_row <- tibble(
  strat_var = "diabetes",   # place overall into diabetes facet
  level     = -1,
  label     = "Overall",
  tau       = tau_overall,
  ci_low    = tau_overall - z_crit * se_overall,
  ci_high   = tau_overall + z_crit * se_overall
)

# ---- 7.2 Stratified rows ----
for (v in strata_vars) {
  for (lvl in 0:1) {
    
    nm  <- paste0(v, "_", lvl)
    est <- sdid_strata[[nm]]
    
    if (is.null(est)) next
    
    if (!inherits(est, "synthdid_estimate")) {
      warning(paste("Object", nm, "is not a synthdid_estimate; skipping"))
      next
    }
    
    tau <- as.numeric(est)
    vc  <- stats::vcov(est, method = se_method)
    se  <- sqrt(as.numeric(vc))
    
    if (is.na(se)) {
      warning(paste("SE is NA for", nm, "- CI will be NA"))
    }
    
    lab <- dplyr::case_when(
      v == "sex"          & lvl == 0 ~ "sex = 0",
      v == "sex"          & lvl == 1 ~ "sex = 1",
      v == "hypertension" & lvl == 0 ~ "hypertension = 0",
      v == "hypertension" & lvl == 1 ~ "hypertension = 1",
      v == "diabetes"     & lvl == 0 ~ "diabetes = 0",
      v == "diabetes"     & lvl == 1 ~ "diabetes = 1",
      v == "smoking"      & lvl == 0 ~ "smoking = 0",
      v == "smoking"      & lvl == 1 ~ "smoking = 1",
      TRUE ~ paste0(v, " = ", lvl)
    )
    
    res_list[[nm]] <- tibble(
      strat_var = v,
      level     = lvl,
      label     = lab,
      tau       = tau,
      ci_low    = tau - z_crit * se,
      ci_high   = tau + z_crit * se
    )
  }
}

forest_df <- dplyr::bind_rows(res_list) %>%
  dplyr::bind_rows(overall_row)

# ------------------------------------------------------------
# 8) Order facets and labels
#    Put Overall at the TOP of diabetes panel
# ------------------------------------------------------------
forest_df <- forest_df %>%
  mutate(
    strat_var = factor(
      strat_var,
      levels = c("diabetes", "sex", "hypertension", "smoking")
    ),
    label = factor(
      label,
      levels = c(
        "smoking = 1", "smoking = 0",
        "diabetes = 1", "diabetes = 0", "Overall",
        "hypertension = 1", "hypertension = 0",
        "sex = 1", "sex = 0"
      )
    )
  )

# Global label order across all facets
# Because ggplot places higher factor levels higher on the y-axis,
# "Overall" is placed after diabetes labels in the levels below,
# then reversed.
# label_levels <- rev(c(
#   "sex = 0", "sex = 1",
#   "hypertension = 0", "hypertension = 1",
#   "diabetes = 0", "diabetes = 1", "Overall",
#   "smoking = 0", "smoking = 1"
# ))
# 
# forest_df <- forest_df %>%
#   mutate(
#     label = factor(label, levels = label_levels)
#   )
# 
# print(forest_df)

# ------------------------------------------------------------
# 9) PNAS-style forest plot
# ------------------------------------------------------------
forest_plot <- ggplot(forest_df, aes(x = tau, y = label)) +
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.3,
    colour     = "grey40"
  ) +
  geom_pointrange(
    aes(xmin = ci_low, xmax = ci_high),
    linewidth = 0.3,
    fatten    = 2
  ) +
  facet_grid(
    strat_var ~ .,
    scales = "free_y",
    space  = "free_y"
  ) +
  labs(
    x = "SDID effect on FBG",
    y = NULL
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme_classic(base_size = 10) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    panel.grid   = element_blank(),
    
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
    strip.text       = element_text(face = "bold", size = 9),
    
    axis.title.x = element_text(face = "bold"),
    axis.text    = element_text(colour = "black", size = 8),
    axis.ticks   = element_line(colour = "black", linewidth = 0.3),
    
    legend.position = "none",
    plot.title      = element_blank(),
    plot.subtitle   = element_blank()
  )

forest_plot <- forest_plot +
  coord_cartesian(xlim = c(-0.15, 0.05))

print(forest_plot)

# ------------------------------------------------------------
# 10) Save PNAS-style figure
# ------------------------------------------------------------
ggsave(
  filename = "fig_fbg_sdid_forest_pnas.pdf",
  plot     = forest_plot,
  width    = 6.5,
  height   = 6,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename = "/Users/yuqing/Downloads/GreenSpace/fig_fbg_sdid_forest_pnas300.tiff",
  plot     = forest_plot,
  width    = 3,
  height   = 4.5,
  units    = "in",
  dpi      = 600
)
