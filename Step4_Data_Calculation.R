# ============================================================
# Build DTO & SEV for BMI, FBG, MAP, PP on panel_16_24_R300
# ============================================================

rm(list = ls())
library(dplyr)

# 1) Read your panel (RDS or CSV) ----
# If it's an RDS:
panel <- readRDS("/Users/yqdai/Downloads/GreenSpace/data/D_all_SI/panel_14_24_R300.rds")

# If instead it's a CSV, use:
# panel <- read.csv("C:/path/to/panel_16_24_R300.csv")

# Ensure lower-case names (your example already looks like this)
panel <- panel %>% 
  rename_with(tolower)

# Quick check
str(panel)
head(panel)

# 2) Helper functions ----
dto_above <- function(x, thr_opt) {
  # distance to optimal: penalty only when above healthy/target threshold
  pmax(x - thr_opt, 0)
}

sev_above <- function(x, thr_sev) {
  # distance to severity: penalty only when above a severe threshold
  pmax(x - thr_sev, 0)
}

# 3) Clinical thresholds (you can tune these) ----
thr <- list(
  # BMI (kg/m^2), Chinese guideline cutoffs:
  # Normal: <24; Overweight: 24–27.9; Obesity: >=28
  bmi_opt = 24,   # "optimal" upper limit
  bmi_sev = 28,   # obesity
  
  # FBG (mmol/L)
  # Normal: <5.6; Diabetes: >=7.0
  fbg_opt = 5.6,  # healthy upper bound you requested
  fbg_sev = 7.0,  # diabetes threshold
  
  # MAP (mmHg) – approximate targets
  map_opt = 100,  # upper limit of acceptable MAP
  map_sev = 110,  # clearly high MAP
  
  # Pulse Pressure (mmHg)
  pp_opt  = 40,   # ideal PP
  pp_sev  = 60    # clearly high PP
)

# 4) Create DTO & SEV variables ----
panel <- panel %>%
  mutate(
    # Ensure these are numeric
    bmi = as.numeric(bmi),
    fbg = as.numeric(fbg),
    map = as.numeric(map),
    pp  = as.numeric(pp),
    
    # BMI
    bmi_dto = dto_above(bmi, thr$bmi_opt),
    bmi_sev = sev_above(bmi, thr$bmi_sev),
    
    # FBG (fasting glucose)
    fbg_dto = dto_above(fbg, thr$fbg_opt),
    fbg_sev = sev_above(fbg, thr$fbg_sev),
    
    # MAP
    map_dto = dto_above(map, thr$map_opt),
    map_sev = sev_above(map, thr$map_sev),
    
    # Pulse pressure
    pp_dto  = dto_above(pp, thr$pp_opt),
    pp_sev  = sev_above(pp, thr$pp_sev)
  )

# 5) Quick sanity checks ----
panel %>%
  summarise(
    mean_bmi      = mean(bmi, na.rm = TRUE),
    mean_bmi_dto  = mean(bmi_dto, na.rm = TRUE),
    prop_bmi_high = mean(bmi_dto > 0, na.rm = TRUE),
    
    mean_fbg      = mean(fbg, na.rm = TRUE),
    mean_fbg_dto  = mean(fbg_dto, na.rm = TRUE),
    prop_fbg_high = mean(fbg_dto > 0, na.rm = TRUE)
  )

# 6) Build SDID/ASCM-ready panel ----
#   - keep only improved/stable 300m greenness
#   - define treated indicator (improved)
#   - keep core covariates and outcomes

panel_ready <- panel %>%
  filter(label_300 %in% c("improved", "stable")) %>%
  mutate(
    treat = if_else(label_300 == "improved", 1L, 0L)
  ) %>%
  select(
    person_id, year,
    label_300, label_1000, treat,
    bmi, fbg, map, pp,
    bmi_dto, bmi_sev,
    fbg_dto, fbg_sev,
    map_dto, map_sev,
    pp_dto, pp_sev,
    hypertension, diabetes, smoking, age, sex,
    lat_gcj, lon_gcj
  )

# (Optional) save this for SDID/ASCM scripts
saveRDS(
  panel_ready,
  "/Users/yqdai/Downloads/GreenSpace/data/D_all_SI/panel_14_24_R300_DTO_SEV.rds"
)

# Or CSV if you prefer
# write.csv(panel_ready,
#           "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/data/panel_16_24_R300_DTO_SEV.csv",
#           row.names = FALSE)
