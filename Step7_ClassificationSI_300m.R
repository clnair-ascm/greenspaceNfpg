# ============================================================
# NDVI cross-year cleanup:
# - Remove rows where NDVI is NA/blank or < threshold in that year
# - Remove those person_id rows from ALL year files
# - Save cleaned files + a removal log to a new folder
# ============================================================

rm(list = ls())

# tidy read: read all Excel files into a named list
# update this path if needed
data_dir <- "/Users/yqdai/Downloads/GreenSpace/data/D_all_SI_greenT/"

my_list <- readRDS("/Users/yqdai/Downloads/GreenSpace/data/D_all_SI_greenT/mydat_16T24NO23_300.rds")
# my_list$complete_2014_cleaned$label_300 <- NULL
# my_list$complete_2014_cleaned$labels_1000 <- NULL
# my_list$complete_2015_cleaned$label_300 <- NULL
# my_list$complete_2015_cleaned$labels_1000 <- NULL
# ---- packages ----
pkgs <- c("dplyr", "purrr", "stringr", "tidyr")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
library(dplyr); library(purrr); library(stringr); library(tidyr)

# ---- helpers ----
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

# ---- 0) Safety check ----
if(!exists("my_list") || !is.list(my_list)) stop("Please ensure 'my_list' is present and is a list of data.frames.")

# ---- 1) Build long NDVI table from my_list ----
# iterate over list elements; names(my_list) used to infer year if needed
list_names <- names(my_list)
if(is.null(list_names)) list_names <- rep(NA_character_, length(my_list))

ndvi_long <- imap_dfr(my_list, function(df, nm) {
  df2 <- as.data.frame(df) |> norm_names()
  id_col <- detect_person_id(df2)
  if(is.na(id_col)) stop("No person id column detected in list element: '", nm, "'.")
  # find ndvi 300m columns (flexible)
  ndvi_cols <- names(df2)[str_detect(names(df2), "ndvi") & str_detect(names(df2), "300")]
  # also accept columns like "ndvi_2016" without 300, as fallback
  if(length(ndvi_cols) == 0) {
    ndvi_cols <- names(df2)[str_detect(names(df2), "^ndvi[_0-9]*")]
  }
  if(length(ndvi_cols) == 0) {
    # If no NDVI columns, return empty tibble for this element
    return(tibble(person_id = character(), year = integer(), ndvi_300m = numeric()))
  }
  # For each NDVI column, attempt to extract year from column name, otherwise use list element name
  rows <- map_dfr(ndvi_cols, function(col) {
    yr_col <- str_extract(col, "\\d{4}") %>% as.integer()
    if(is.na(yr_col)) yr_col <- get_year_from_name(nm)
    tibble(
      person_id = as.character(df2[[id_col]]),
      year      = yr_col,
      ndvi_300m = as.numeric(df2[[col]])
    )
  })
  # if multiple ndvi_cols yielded multiple rows per person-year, keep the first non-NA per person-year
  rows %>%
    group_by(person_id, year) %>%
    summarise(ndvi_300m = {
      v <- ndvi_300m
      v[!is.na(v)][1]  # first non-NA
    }, .groups = "drop")
}, .id = "source")

# keep only years in target range and non-missing person_id
ndvi_long <- ndvi_long %>%
  filter(!is.na(person_id)) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2016 & year <= 2024)

# ---- 2) Pivot wide (one row per person) if you want ----
ndvi_wide <- ndvi_long %>%
  select(person_id, year, ndvi_300m) %>%    # drop 'source' and any other unwanted cols
  pivot_wider(
    names_from  = year,
    values_from = ndvi_300m,
    names_prefix = "ndvi_"
  )

# ---- packages ----
library(dplyr); library(tidyr); library(stringr); library(purrr)

ts_slope <- function(y, t) {
  ok <- is.finite(y) & is.finite(t)
  y <- as.numeric(y[ok]); t <- as.numeric(t[ok])
  n <- length(y); if (n < 2) return(NA_real_)
  idx <- utils::combn(n, 2)
  slopes <- (y[idx[2, ]] - y[idx[1, ]]) / (t[idx[2, ]] - t[idx[1, ]])
  stats::median(slopes[is.finite(slopes)])
}

classify_ndvi_groups_v2 <- function(ndvi_wide,
                                    # cleaning
                                    min_ndvi = 0.05, strict_clean = TRUE, clamp_01 = FALSE,
                                    # windows
                                    years_pre = 2016:2020, years_post = 2021:2024,
                                    min_pre = 3, min_post = 2,
                                    # thresholds
                                    q_noise = 0.90, q_imp = 0.85,
                                    floor_improve = 0.02,
                                    slope_quant = 0.80, sd_quant = 0.80,
                                    slope_floor = 0.01, sd_floor = 0.03,
                                    # persistence
                                    persistence_share = 0.50,  # “half the post years”
                                    donor_band_sd = 2) {
  
  if (!"person_id" %in% names(ndvi_wide)) stop("ndvi_wide must contain 'person_id'.")
  yr_cols <- names(ndvi_wide)[grepl("^ndvi_\\d{4}$", names(ndvi_wide))]
  if (!length(yr_cols)) stop("No ndvi_YYYY columns found.")
  avail_years <- sort(unique(as.integer(str_extract(yr_cols, "\\d{4}"))))
  pre_years  <- intersect(years_pre,  avail_years)
  post_years <- intersect(years_post, avail_years)
  if (length(post_years) < min_post) stop("Insufficient post years.")
  
  nw <- ndvi_wide
  if (clamp_01) nw[yr_cols] <- lapply(nw[yr_cols], function(x) pmax(pmin(as.numeric(x),1),0))
  
  if (strict_clean) {
    nw <- nw %>%
      filter(if_all(all_of(yr_cols), ~ !is.na(.x) & .x >= min_ndvi))
  }
  
  ndvi_long <- nw %>%
    pivot_longer(all_of(yr_cols), names_to="var", values_to="ndvi_300m") %>%
    mutate(year = as.integer(str_extract(var, "\\d{4}")),
           person_id = as.character(person_id)) %>%
    select(person_id, year, ndvi_300m) %>%
    filter(year %in% c(pre_years, post_years))
  
  # summaries (include post_n here)
  ndvi_pre <- ndvi_long %>%
    filter(year %in% pre_years) %>%
    arrange(person_id, year) %>%
    group_by(person_id) %>%
    summarise(pre_mean = mean(ndvi_300m, na.rm=TRUE),
              pre_sd   = sd(ndvi_300m,   na.rm=TRUE),
              pre_slope= ts_slope(ndvi_300m, year),
              pre_n    = sum(!is.na(ndvi_300m)), .groups="drop")
  
  ndvi_post <- ndvi_long %>%
    filter(year %in% post_years) %>%
    group_by(person_id) %>%
    summarise(post_mean = mean(ndvi_300m, na.rm=TRUE),
              post_n    = sum(!is.na(ndvi_300m)), .groups="drop")
  
  ndvi_sum <- ndvi_pre %>% full_join(ndvi_post, by="person_id") %>%
    mutate(delta = post_mean - pre_mean)
  
  # placebo
  placebo <- ndvi_long %>%
    filter(year %in% pre_years) %>%
    arrange(person_id, year) %>%
    group_by(person_id) %>%
    summarise(d_pre = {
      yi <- ndvi_300m; ti <- year
      ok <- is.finite(yi); yi <- yi[ok]; ti <- ti[ok]
      n <- length(yi)
      if (n < 2) NA_real_ else {
        K <- min(3L, floor(n/2))
        m1 <- mean(yi[seq_len(K)], na.rm=TRUE)
        m2 <- mean(yi[(n-K+1):n],   na.rm=TRUE)
        m2 - m1
      }
    }, .groups="drop")
  
  placebo_vals <- placebo$d_pre[is.finite(placebo$d_pre)]
  noise_thr <- if (length(placebo_vals)>=10) quantile(abs(placebo_vals), q_noise, na.rm=TRUE) else 0.02
  imp_thr   <- if (length(placebo_vals)>=10) max(quantile(placebo_vals, q_imp, na.rm=TRUE), floor_improve) else floor_improve
  
  # persistence margin = improvement_thr (key modelling choice)
  persist_margin <- imp_thr
  
  # pre-stability thresholds (data-driven with floors)
  slope_vals <- ndvi_sum$pre_slope[is.finite(ndvi_sum$pre_slope)]
  sd_vals    <- ndvi_sum$pre_sd[is.finite(ndvi_sum$pre_sd)]
  slope_thr  <- if (length(slope_vals)>=10) max(quantile(abs(slope_vals), slope_quant, na.rm=TRUE), slope_floor) else slope_floor
  sd_thr     <- if (length(sd_vals)>=10)    max(quantile(sd_vals,       sd_quant,  na.rm=TRUE), sd_floor)       else sd_floor
  
  # persistence (count-based) — NOTE: rename to avoid collision
  pers <- ndvi_long %>%
    left_join(ndvi_pre %>% select(person_id, pre_mean), by="person_id") %>%
    filter(year %in% post_years) %>%
    mutate(above = ndvi_300m >= (pre_mean + persist_margin),
           below = ndvi_300m <= (pre_mean - persist_margin)) %>%
    group_by(person_id) %>%
    summarise(
      post_n_pers = n(),                  # <-- renamed
      n_up   = sum(above, na.rm=TRUE),
      n_down = sum(below, na.rm=TRUE),
      post_share_up   = ifelse(post_n_pers>0, n_up/post_n_pers, NA_real_),
      post_share_down = ifelse(post_n_pers>0, n_down/post_n_pers, NA_real_),
      required_up   = pmax(1L, ceiling(persistence_share * post_n_pers)),
      required_down = pmax(1L, ceiling(persistence_share * post_n_pers)),
      .groups="drop"
    )
  
  labels <- ndvi_sum %>%
    # keep left names; any duplicates from pers get ".pers"
    left_join(pers, by="person_id", suffix = c("", ".pers")) %>%
    mutate(
      pre_stable = (pre_n >= min_pre) & (post_n >= min_post) &
        is.finite(pre_slope) & is.finite(pre_sd) &
        (abs(pre_slope) <= slope_thr) & (pre_sd <= sd_thr),
      
      improved = pre_stable & (delta >= imp_thr) & (n_up   >= required_up),
      declined = pre_stable & (delta <= -imp_thr) & (n_down >= required_down),
      stable   = pre_stable & (abs(delta) <= noise_thr),
      
      ndvi_group  = case_when(improved ~ "improved",
                              stable   ~ "stable",
                              TRUE     ~ "ambiguous_unstable"),
      ndvi_group4 = case_when(improved ~ "improved",
                              declined ~ "declined",
                              stable   ~ "stable",
                              TRUE     ~ "ambiguous_unstable"),
      treated = ndvi_group == "improved"
    )
  
  # donor common support among controls
  if (any(labels$treated, na.rm=TRUE)) {
    mu_t <- mean(labels$pre_mean[labels$treated], na.rm=TRUE)
    sd_t <- sd(labels$pre_mean[labels$treated], na.rm=TRUE); if (is.finite(sd_t) && sd_t==0) sd_t <- 1e-6
    labels <- labels %>%
      mutate(donor_common_support = (!treated) &
               is.finite(pre_mean) & abs(pre_mean - mu_t) <= donor_band_sd * sd_t)
  } else {
    labels$donor_common_support <- NA
  }
  
  thresholds <- list(
    years_pre = pre_years, years_post = post_years,
    min_ndvi = min_ndvi,
    noise_threshold = noise_thr,
    improvement_thr = imp_thr,
    persistence_margin = persist_margin,
    slope_threshold = slope_thr,
    sd_threshold = sd_thr,
    min_pre = min_pre, min_post = min_post,
    persistence_req = persistence_share,
    donor_band_sd = donor_band_sd
  )
  diagnostics <- list(
    counts_3grp = table(labels$ndvi_group,  useNA="ifany"),
    counts_4grp = table(labels$ndvi_group4, useNA="ifany")
  )
  list(labels = labels, thresholds = thresholds, diagnostics = diagnostics,
       ndvi_long = ndvi_long, ndvi_wide_clean = nw)
}

# ---- Run with your current settings (same as before, but v2 margin/persistence) ----
res <- classify_ndvi_groups_v2(ndvi_wide, strict_clean = TRUE, min_ndvi = 0.05)
table(res$labels$ndvi_group)
res$thresholds


# ndvi_wide_labeled <- ndvi_wide %>% dplyr::left_join(ndvi_labels, by = "person_id")
labels_300 <- res$labels %>%
  mutate(
    label_300 = case_when(
      ndvi_group == "improved"            ~ "improved",
      ndvi_group == "stable"              ~ "stable",
      ndvi_group == "ambiguous_unstable"  ~ "ambiguous",
      TRUE                                ~ NA_character_
    )
  ) %>%
  select(person_id, label_300)

# labels_300 <- res$labels %>%
#   mutate(
#     label_300 = case_when(
#       ndvi_group4 == "improved"           ~ "improved",
#       ndvi_group4 == "declined"           ~ "declined",
#       ndvi_group4 == "stable"             ~ "stable",
#       ndvi_group4 == "ambiguous_unstable" ~ "ambiguous",
#       TRUE                                ~ NA_character_
#     )
#   ) %>%
#   select(person_id, label_300)


# (Optional) label every data frame in your original my_list
add_label_300_to_list <- function(my_list, labels_300, id_col = "person_id") {
  lapply(my_list, function(df) {
    df2 <- df
    
    # detect ID column if `id_col` not present
    if (!(id_col %in% names(df2))) {
      candidates <- c("person_id","personid","pid","id","subject_id","subject",
                      "participant_id","case_id")
      hit <- candidates[candidates %in% names(df2)]
      if (!length(hit)) stop("No ID column found in one list element.")
      id_col <- hit[1]
    }
    
    df2[[id_col]] <- as.character(df2[[id_col]])
    
    df2 <- df2 %>%
      left_join(
        labels_300 %>% rename(!!id_col := person_id),
        by = id_col
      )
    
    df2
  })
}

# apply
my_list_labeled <- add_label_300_to_list(my_list, labels_300)

saveRDS(
  my_list_labeled,
  file = "/Users/yqdai/Downloads/GreenSpace/data/D_all_SI_greenT/mydat_16T24NO23_300.rds"
)

