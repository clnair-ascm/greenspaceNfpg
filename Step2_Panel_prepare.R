rm(list = ls())

library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# Suppose you have a list of yearly data frames as before:
my_list_labeled <- readRDS("/Users/yqdai/Downloads/GreenSpace/data/D_all_SI_greenT/mydat_16T24NO23_loose.rds")

norm_names <- function(df) {
  nms <- names(df)
  nms2 <- nms %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
  names(df) <- nms2
  df
}

detect_person_id <- function(df) {
  candidates <- c("person_id","personid","pid","id","subject_id","subjectid",
                  "subject","participant_id","participantid","participant","case_id","caseid")
  nms <- names(df)
  found <- intersect(candidates, nms)
  if(length(found)) return(found[1])
  regex_match <- nms[str_detect(nms, "\\b(id|pid|person)\\b")]
  if(length(regex_match)) return(regex_match[1])
  NA_character_
}

get_year_from_name <- function(nm) {
  yr <- str_extract(nm, "\\d{4}")
  if(!is.na(yr)) as.integer(yr) else NA_integer_
}

panel <- imap_dfr(my_list_labeled, function(df, nm) {
  df2 <- as.data.frame(df) |> norm_names()
  id_col <- detect_person_id(df2)
  if (is.na(id_col)) stop("No person id column in: ", nm)
  
  # year: from column if exists, otherwise from name
  year_col <- dplyr::case_when(
    "year"      %in% names(df2) ~ "year",
    "exam_year" %in% names(df2) ~ "exam_year",
    TRUE                         ~ NA_character_
  )
  year_vec <- if (!is.na(year_col)) as.integer(df2[[year_col]]) else rep(get_year_from_name(nm), nrow(df2))
  
  tibble(
    person_id = as.character(df2[[id_col]]),
    year      = year_vec,
    label_300 = as.character(df2$label_300),
    label_1000 = as.character(df2$labels_1000),
    bmi   = as.numeric(df2$bmi),
    fbg   = as.numeric(df2$fasting_glucose),
    map   = as.numeric(df2$map),
    pp    = as.numeric(df2$pulse_pressure),
    hypertension = as.numeric(df2$hypertension),
    diabetes = as.numeric(df2$diabetes),
    smoking = as.numeric(df2$smoking),
    bmi_dto   = as.numeric(df2$bmi_dto),
    fbg_dto   = as.numeric(df2$fpg_dto),
    map_dto   = as.numeric(df2$map_dto),
    bmi_sev   = as.numeric(df2$bmi_sev),
    fbg_sev   = as.numeric(df2$fpg_sev),
    map_sev   = as.numeric(df2$map_sev),
    waist_dto = as.numeric(df2$waist_dto),
    tc_dto    = as.numeric(df2$tc_dto),
    pp_dto    = as.numeric(df2$pp_dto),
    age       = as.numeric(df2$age),
    sex       = as.numeric(df2$sex),
    lat_gcj   = as.numeric(df2$lat_gcj),
    lon_gcj   = as.numeric(df2$lon_gcj)
  )
})

panel_r300 <- panel %>%
  filter(year >= 2014, year <= 2024,
         label_300 %in% c("improved", "stable"))
panel_r1000 <- panel %>%
  filter(year >= 2014, year <= 2024,
         label_1000 %in% c("improved", "stable"))

# saveRDS(panel_r300,
#         "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/data/panel_16_24_R300.rds")
# 
# saveRDS(panel_r1000,
#         "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/data/panel_16_24_R1000.rds")

# library(openxlsx)
# 
# write.xlsx(
#   panel_r300,
#   "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/data/panel_16_24_R300.xlsx"
# )
# 
# write.xlsx(
#   panel_r1000,
#   "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/data/panel_16_24_R1000.xlsx"
# )

saveRDS(panel_r300,  "/Users/yqdai/Downloads/GreenSpace/data/D_all_SI_greenT/panel_16_24_R300.rds")

saveRDS(panel_r1000, "/Users/yqdai/Downloads/GreenSpace/data/D_all_SI_greenT/panel_16_24_R1000.rds")
# 
# library(dplyr)
# 
# # 1. Restrict to 2020 only
# panel_2020 <- panel %>%
#   filter(year == 2020)
# 
# # Quick check
# table(panel_2020$label_300, useNA = "ifany")
# 
# # 2. Summaries by group
# age_sex_summary <- panel_2020 %>%
#   group_by(label_300) %>%
#   summarise(
#     n        = n(),
#     age_mean = mean(age, na.rm = TRUE),
#     age_sd   = sd(age, na.rm = TRUE),
#     sex_mean = mean(sex, na.rm = TRUE),  # if sex is coded 0/1, this is proportion of "1"
#     .groups  = "drop"
#   )
# 
# age_sex_summary
