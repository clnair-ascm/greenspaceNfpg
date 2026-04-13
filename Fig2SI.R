# ============================================================
# SDID for FBG: Improved vs Stable (NDVI 300m, 2016–2024)
# Whole population (no strata)
# ============================================================

rm(list = ls())

# ---- packages ----
library(dplyr)
library(readr)
library(synthdid)
library(ggplot2)

# ------------------------------------------------------------
# 1) Read panel data
# ------------------------------------------------------------
panel <- readRDS(
  "/Users/yqdai/Downloads/GreenSpace/data/panel_14_24_R300_DTO_SEV.rds"
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
    !is.na(fbg)              # outcome: fasting blood glucose
  ) %>%
  mutate(
    pre  = year <= 2020,
    post = year >= 2021
  ) %>%
  group_by(person_id) %>%
  # require: each person has at least 3 pre + 2 post FBG measures
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
# 4) Build synthdid panel matrices and fit overall SDID
# ------------------------------------------------------------
pm_fbg <- panel.matrices(
  panel     = as.data.frame(panel_fbg),
  unit      = "person_id",
  time      = "year",
  outcome   = "fbg_dto",       # outcome: fasting blood glucose
  treatment = "treat_time"
)

Y  <- pm_fbg$Y
N0 <- pm_fbg$N0
T0 <- pm_fbg$T0

cat("Overall Y dim:", dim(Y), "\n")
cat("Overall N0 (controls):", N0, "\n")
cat("Overall T0 (pre-periods):", T0, "\n")

sdid_fbg_overall <- synthdid_estimate(Y, N0, T0)

# ============================================================
# Figure 1: Raw FBG trajectories (Improved vs Stable)
# ============================================================

fig1_df <- panel_fbg %>%
  group_by(label_300, year) %>%
  summarise(
    mean_fbg = mean(fbg, na.rm = TRUE),
    se_fbg   = sd(fbg, na.rm = TRUE) / sqrt(sum(!is.na(fbg))),
    .groups  = "drop"
  )

# ============================================================
# Figure 1: Raw FBG trajectories (Improved vs Stable) – PNAS style
# ============================================================

gg_fig1 <- ggplot(
  fig1_df,
  aes(x = year, y = mean_fbg, colour = label_300, group = label_300)
) +
  # Intervention line
  geom_vline(xintercept = 2021, linetype = "dashed", colour = "grey60") +
  # 95% CI ribbon
  geom_ribbon(
    aes(ymin = mean_fbg - 1.96 * se_fbg,
        ymax = mean_fbg + 1.96 * se_fbg,
        fill = label_300),
    alpha = 0.15, colour = NA, show.legend = FALSE
  ) +
  # Mean trajectories
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  # Colour-blind safe palette, distinct from NDVI figure
  scale_colour_manual(
    values = c("stable" = "#E69F00", "improved" = "#009E73"),
    labels = c("Stable", "Improved"),
    name   = "Greenness group (300 m)"
  ) +
  scale_fill_manual(
    values = c("stable" = "#F6C96B", "improved" = "#8FD3BF")
  ) +
  labs(
    title    = "Mean fasting glucose by greenness change, 2016–2024",
    subtitle = "NDVI 300 m: Improved vs Stable; dashed line marks 2021 (post period)",
    x        = "Year",
    y        = "Fasting glucose (mmol/L)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text         = element_text(size = 12, colour = "grey20"),
    axis.title        = element_text(size = 14),
    legend.position   = "top",
    legend.title      = element_text(size = 13),
    legend.text       = element_text(size = 12),
    plot.title        = element_text(face = "bold", size = 18),
    plot.subtitle     = element_text(size = 13),
    plot.margin       = margin(t = 5, r = 10, b = 5, l = 10)
  )

gg_fig1

# ============================================================
# Figure 2: Overall SDID trajectories (FBG, whole population) – PNAS style
# ============================================================

gg_fig2 <- synthdid::synthdid_plot(
  sdid_fbg_overall,
  facet.vertical = FALSE
) +
  ggtitle(
    "Synthetic difference-in-differences for fasting glucose",
    subtitle = "Improved NDVI (treated) vs synthetic control, 2016–2024"
  ) +
  xlab("Year") +
  ylab("Fasting glucose (mmol/L)") +
  theme_classic(base_size = 16) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text         = element_text(size = 12, colour = "grey20"),
    axis.title        = element_text(size = 14),
    legend.position   = "top",
    legend.title      = element_text(size = 13),
    legend.text       = element_text(size = 12),
    plot.title        = element_text(face = "bold", size = 18),
    plot.subtitle     = element_text(size = 13),
    plot.margin       = margin(t = 5, r = 10, b = 5, l = 10)
  )

gg_fig2





# ============================================================
# (Optional) Save high-resolution versions for PNAS
# ============================================================

# Figure 1: raw trajectories
ggsave(
  filename = "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/fig/FBG_raw_traj_improved_vs_stable.tiff",
  plot     = gg_fig1,
  width    = 6,
  height   = 4.5,
  units    = "in",
  dpi      = 600
)

# Figure 2: SDID trajectories
ggsave(
  filename = "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/fig/FBG_sdid_overall.tiff",
  plot     = gg_fig2,
  width    = 6,
  height   = 4.5,
  units    = "in",
  dpi      = 600
)






# ============================================================
# SDID for FBG: Improved vs Stable (NDVI 300m, 2016–2024)
# Whole population (no strata)
# ============================================================

rm(list = ls())

# ---- packages ----
library(dplyr)
library(readr)
library(synthdid)
library(ggplot2)

# ------------------------------------------------------------
# 1) Read panel data
# ------------------------------------------------------------
panel <- readRDS(
  "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/data/panel_16_24_R300_DTO_SEV.rds"
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
    label_1000 = as.character(label_1000)
  )

# ------------------------------------------------------------
# 2) Restrict to Improved vs Stable, FBG available, year window
# ------------------------------------------------------------
panel_fbg <- panel %>%
  filter(
    label_1000 %in% c("improved", "stable"),
    year >= 2016, year <= 2024,
    !is.na(fbg)              # outcome: fasting blood glucose
  ) %>%
  mutate(
    pre  = year <= 2020,
    post = year >= 2021
  ) %>%
  group_by(person_id) %>%
  # require: each person has at least 3 pre + 2 post FBG measures
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
    treat_time = if_else(label_1000 == "improved" & year >= 2021, 1L, 0L)
  ) %>%
  arrange(person_id, year)

# ------------------------------------------------------------
# 4) Build synthdid panel matrices and fit overall SDID
# ------------------------------------------------------------
pm_fbg <- panel.matrices(
  panel     = as.data.frame(panel_fbg),
  unit      = "person_id",
  time      = "year",
  outcome   = "fbg",       # outcome: fasting blood glucose
  treatment = "treat_time"
)

Y  <- pm_fbg$Y
N0 <- pm_fbg$N0
T0 <- pm_fbg$T0

cat("Overall Y dim:", dim(Y), "\n")
cat("Overall N0 (controls):", N0, "\n")
cat("Overall T0 (pre-periods):", T0, "\n")

sdid_fbg_overall <- synthdid_estimate(Y, N0, T0)

# ============================================================
# Figure 1: Raw FBG trajectories (Improved vs Stable)
# ============================================================

fig1_df <- panel_fbg %>%
  group_by(label_1000, year) %>%
  summarise(
    mean_fbg = mean(fbg, na.rm = TRUE),
    se_fbg   = sd(fbg, na.rm = TRUE) / sqrt(sum(!is.na(fbg))),
    .groups  = "drop"
  )

# ============================================================
# Figure 1: Raw FBG trajectories (Improved vs Stable) – PNAS style
# ============================================================

gg_fig1 <- ggplot(
  fig1_df,
  aes(x = year, y = mean_fbg, colour = label_1000, group = label_1000)
) +
  # Intervention line
  geom_vline(xintercept = 2021, linetype = "dashed", colour = "grey60") +
  # 95% CI ribbon
  geom_ribbon(
    aes(ymin = mean_fbg - 1.96 * se_fbg,
        ymax = mean_fbg + 1.96 * se_fbg,
        fill = label_1000),
    alpha = 0.15, colour = NA, show.legend = FALSE
  ) +
  # Mean trajectories
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  # Colour-blind safe palette, distinct from NDVI figure
  scale_colour_manual(
    values = c("stable" = "#E69F00", "improved" = "#009E73"),
    labels = c("Stable", "Improved"),
    name   = "Greenness group (300 m)"
  ) +
  scale_fill_manual(
    values = c("stable" = "#F6C96B", "improved" = "#8FD3BF")
  ) +
  labs(
    title    = "Mean fasting glucose by greenness change, 2016–2024",
    subtitle = "NDVI 1000 m: Improved vs Stable; dashed line marks 2021 (post period)",
    x        = "Year",
    y        = "Fasting glucose (mmol/L)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text         = element_text(size = 12, colour = "grey20"),
    axis.title        = element_text(size = 14),
    legend.position   = "top",
    legend.title      = element_text(size = 13),
    legend.text       = element_text(size = 12),
    plot.title        = element_text(face = "bold", size = 18),
    plot.subtitle     = element_text(size = 13),
    plot.margin       = margin(t = 5, r = 10, b = 5, l = 10)
  )

gg_fig1

# ============================================================
# Figure 2: Overall SDID trajectories (FBG, whole population) – PNAS style
# ============================================================

gg_fig2 <- synthdid::synthdid_plot(
  sdid_fbg_overall,
  facet.vertical = FALSE
) +
  ggtitle(
    "Synthetic difference-in-differences for fasting glucose",
    subtitle = "Improved NDVI (treated) vs synthetic control, 2016–2024"
  ) +
  xlab("Year") +
  ylab("Fasting glucose (mmol/L)") +
  theme_classic(base_size = 16) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.text         = element_text(size = 12, colour = "grey20"),
    axis.title        = element_text(size = 14),
    legend.position   = "top",
    legend.title      = element_text(size = 13),
    legend.text       = element_text(size = 12),
    plot.title        = element_text(face = "bold", size = 18),
    plot.subtitle     = element_text(size = 13),
    plot.margin       = margin(t = 5, r = 10, b = 5, l = 10)
  )

gg_fig2





# ============================================================
# (Optional) Save high-resolution versions for PNAS
# ============================================================

# Figure 1: raw trajectories
ggsave(
  filename = "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/fig/FBG_raw_traj_improved_vs_stable.tiff",
  plot     = gg_fig1,
  width    = 6,
  height   = 4.5,
  units    = "in",
  dpi      = 600
)

# Figure 2: SDID trajectories
ggsave(
  filename = "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/fig/FBG_sdid_overall.tiff",
  plot     = gg_fig2,
  width    = 6,
  height   = 4.5,
  units    = "in",
  dpi      = 600
)