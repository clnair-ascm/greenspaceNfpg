## ============================================================
## Full pipeline:
## - Read NDVI data list (mydat_16T24NO23.rds)
## - Build 300 m & 1000 m NDVI panels
## - Classify improved vs stable using time-series rules
## - Plot trajectories (300 m & 1000 m) in two PNAS-style panels
## ============================================================

rm(list = ls())

## ---- packages ----
pkgs <- c("dplyr", "purrr", "stringr", "tidyr", "ggplot2")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)

## ---- paths ----
data_dir <- "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/data/"
rds_path <- file.path(data_dir, "mydat_16T24NO23.rds")

## ---- helpers ----
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
  if (length(found)) return(found[1])
  regex_match <- nms[str_detect(nms, "\\b(id|pid|person)\\b")]
  if (length(regex_match)) return(regex_match[1])
  NA_character_
}

get_year_from_name <- function(nm) {
  yr <- str_extract(nm, "\\d{4}")
  if (!is.na(yr)) as.integer(yr) else NA_integer_
}

ts_slope <- function(y, t) {
  ok <- is.finite(y) & is.finite(t)
  y <- as.numeric(y[ok]); t <- as.numeric(t[ok])
  n <- length(y); if (n < 2) return(NA_real_)
  idx <- utils::combn(n, 2)
  slopes <- (y[idx[2, ]] - y[idx[1, ]]) / (t[idx[2, ]] - t[idx[1, ]])
  stats::median(slopes[is.finite(slopes)])
}

## ============================================================
## PART A: 300 m NDVI
## ============================================================

my_list_300 <- readRDS(rds_path)
if (!is.list(my_list_300)) stop("my_list_300 must be a list of data frames.")

# 1) Build long NDVI 300 m table
list_names_300 <- names(my_list_300)
if (is.null(list_names_300)) list_names_300 <- rep(NA_character_, length(my_list_300))

ndvi_long_300 <- imap_dfr(my_list_300, function(df, nm) {
  df2 <- as.data.frame(df) |> norm_names()
  id_col <- detect_person_id(df2)
  if (is.na(id_col)) stop("No person id column detected in list element (300 m): '", nm, "'.")
  
  # find ndvi 300m columns (flexible)
  ndvi_cols <- names(df2)[str_detect(names(df2), "ndvi") & str_detect(names(df2), "300")]
  # also accept columns like "ndvi_2016" without 300, as fallback
  if (length(ndvi_cols) == 0) {
    ndvi_cols <- names(df2)[str_detect(names(df2), "^ndvi[_0-9]*")]
  }
  if (length(ndvi_cols) == 0) {
    return(tibble(person_id = character(), year = integer(), ndvi_300m = numeric()))
  }
  
  # For each NDVI column, attempt to extract year from column name, otherwise use list element name
  rows <- map_dfr(ndvi_cols, function(col) {
    yr_col <- str_extract(col, "\\d{4}") %>% as.integer()
    if (is.na(yr_col)) yr_col <- get_year_from_name(nm)
    tibble(
      person_id = as.character(df2[[id_col]]),
      year      = yr_col,
      ndvi_300m = as.numeric(df2[[col]])
    )
  })
  
  rows %>%
    group_by(person_id, year) %>%
    summarise(ndvi_300m = {
      v <- ndvi_300m
      v[!is.na(v)][1]
    }, .groups = "drop")
}, .id = "source")

ndvi_long_300 <- ndvi_long_300 %>%
  filter(!is.na(person_id)) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2014 & year <= 2024)

# 2) Pivot wide
ndvi_wide_300 <- ndvi_long_300 %>%
  select(person_id, year, ndvi_300m) %>%
  pivot_wider(
    names_from  = year,
    values_from = ndvi_300m,
    names_prefix = "ndvi_"
  )

# 3) Classification function (300 m)
classify_ndvi_groups_v2_300 <- function(ndvi_wide,
                                        min_ndvi = 0.05, strict_clean = TRUE, clamp_01 = FALSE,
                                        years_pre = 2014:2020, years_post = 2021:2024,
                                        min_pre = 3, min_post = 2,
                                        q_noise = 0.90, q_imp = 0.95,
                                        floor_improve = 0.03,
                                        slope_quant = 0.80, sd_quant = 0.80,
                                        slope_floor = 0.01, sd_floor = 0.03,
                                        persistence_share = 0.50,
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
  
  persist_margin <- imp_thr
  
  slope_vals <- ndvi_sum$pre_slope[is.finite(ndvi_sum$pre_slope)]
  sd_vals    <- ndvi_sum$pre_sd[is.finite(ndvi_sum$pre_sd)]
  slope_thr  <- if (length(slope_vals)>=10) max(quantile(abs(slope_vals), slope_quant, na.rm=TRUE), slope_floor) else slope_floor
  sd_thr     <- if (length(sd_vals)>=10)    max(quantile(sd_vals,       sd_quant,  na.rm=TRUE), sd_floor)       else sd_floor
  
  pers <- ndvi_long %>%
    left_join(ndvi_pre %>% select(person_id, pre_mean), by="person_id") %>%
    filter(year %in% post_years) %>%
    mutate(above = ndvi_300m >= (pre_mean + persist_margin),
           below = ndvi_300m <= (pre_mean - persist_margin)) %>%
    group_by(person_id) %>%
    summarise(
      post_n_pers = n(),
      n_up   = sum(above, na.rm=TRUE),
      n_down = sum(below, na.rm=TRUE),
      post_share_up   = ifelse(post_n_pers>0, n_up/post_n_pers, NA_real_),
      post_share_down = ifelse(post_n_pers>0, n_down/post_n_pers, NA_real_),
      required_up   = pmax(1L, ceiling(persistence_share * post_n_pers)),
      required_down = pmax(1L, ceiling(persistence_share * post_n_pers)),
      .groups="drop"
    )
  
  labels <- ndvi_sum %>%
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

res_300 <- classify_ndvi_groups_v2_300(ndvi_wide_300, strict_clean = TRUE, min_ndvi = 0.05)
table(res_300$labels$ndvi_group)

## ============================================================
## PART B: 1000 m NDVI
## ============================================================

my_list_1000 <- readRDS(rds_path)
if (!is.list(my_list_1000)) stop("my_list_1000 must be a list of data frames.")

list_names_1000 <- names(my_list_1000)
if (is.null(list_names_1000)) list_names_1000 <- rep(NA_character_, length(my_list_1000))

ndvi_long_1000 <- imap_dfr(my_list_1000, function(df, nm) {
  df2 <- as.data.frame(df) |> norm_names()
  
  id_col <- detect_person_id(df2)
  if (is.na(id_col)) stop("No person id column detected in list element (1000 m): '", nm, "'.")
  if (!"year" %in% names(df2)) stop("No 'year' column detected in list element (1000 m): '", nm, "'.")
  
  ndvi_cols <- grep("^ndvi_1000m_", names(df2), value = TRUE)
  
  if (length(ndvi_cols) == 0L) {
    message("No ndvi_1000m_* column in ", nm,
            "; NDVI columns present: ",
            paste(grep("ndvi", names(df2), value = TRUE), collapse = ", "))
    return(tibble(
      source    = character(),
      person_id = character(),
      year      = integer(),
      ndvi_1000m = numeric()
    ))
  }
  
  if (length(ndvi_cols) > 1L) {
    message("Multiple ndvi_1000m_* columns in ", nm, 
            ": ", paste(ndvi_cols, collapse = ", "),
            ". Using the first.")
  }
  
  col_use <- ndvi_cols[1L]
  
  tibble(
    source     = nm,
    person_id  = as.character(df2[[id_col]]),
    year       = as.integer(df2$year),
    ndvi_1000m = as.numeric(df2[[col_use]])
  )
}, .id = "source_id")

ndvi_long_1000 <- ndvi_long_1000 %>%
  filter(!is.na(person_id)) %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2014 & year <= 2024)

ndvi_wide_1000 <- ndvi_long_1000 %>%
  select(person_id, year, ndvi_1000m) %>%
  pivot_wider(
    names_from  = year,
    values_from = ndvi_1000m,
    names_prefix = "ndvi_"
  )

# Classification function (1000 m)
classify_ndvi_groups_v2_1000 <- function(ndvi_wide,
                                         min_ndvi = 0.05, strict_clean = TRUE, clamp_01 = FALSE,
                                         years_pre = 2014:2020, years_post = 2021:2024,
                                         min_pre = 3, min_post = 2,
                                         q_noise = 0.90, q_imp = 0.95,
                                         floor_improve = 0.03,
                                         slope_quant = 0.80, sd_quant = 0.80,
                                         slope_floor = 0.01, sd_floor = 0.03,
                                         persistence_share = 0.50,
                                         donor_band_sd = 2) {
  
  if (!"person_id" %in% names(ndvi_wide))
    stop("ndvi_wide must contain 'person_id'.")
  
  yr_cols <- names(ndvi_wide)[grepl("^ndvi_\\d{4}$", names(ndvi_wide))]
  if (!length(yr_cols)) stop("No ndvi_YYYY columns found.")
  
  avail_years <- sort(unique(as.integer(str_extract(yr_cols, "\\d{4}"))))
  pre_years  <- intersect(years_pre,  avail_years)
  post_years <- intersect(years_post, avail_years)
  if (length(post_years) < min_post) stop("Insufficient post years.")
  
  nw <- ndvi_wide
  
  if (clamp_01) {
    nw[yr_cols] <- lapply(nw[yr_cols],
                          function(x) pmax(pmin(as.numeric(x), 1), 0))
  }
  
  if (strict_clean) {
    nw <- nw %>%
      filter(if_all(all_of(yr_cols), ~ !is.na(.x) & .x >= min_ndvi))
  }
  
  ndvi_long <- nw %>%
    pivot_longer(all_of(yr_cols),
                 names_to  = "var",
                 values_to = "ndvi_1000m") %>%
    mutate(
      year      = as.integer(str_extract(var, "\\d{4}")),
      person_id = as.character(person_id)
    ) %>%
    select(person_id, year, ndvi_1000m) %>%
    filter(year %in% c(pre_years, post_years))
  
  ndvi_pre <- ndvi_long %>%
    filter(year %in% pre_years) %>%
    arrange(person_id, year) %>%
    group_by(person_id) %>%
    summarise(
      pre_mean  = mean(ndvi_1000m, na.rm = TRUE),
      pre_sd    = sd(ndvi_1000m,   na.rm = TRUE),
      pre_slope = ts_slope(ndvi_1000m, year),
      pre_n     = sum(!is.na(ndvi_1000m)),
      .groups   = "drop"
    )
  
  ndvi_post <- ndvi_long %>%
    filter(year %in% post_years) %>%
    group_by(person_id) %>%
    summarise(
      post_mean = mean(ndvi_1000m, na.rm = TRUE),
      post_n    = sum(!is.na(ndvi_1000m)),
      .groups   = "drop"
    )
  
  ndvi_sum <- ndvi_pre %>%
    full_join(ndvi_post, by = "person_id") %>%
    mutate(delta = post_mean - pre_mean)
  
  placebo <- ndvi_long %>%
    filter(year %in% pre_years) %>%
    arrange(person_id, year) %>%
    group_by(person_id) %>%
    summarise(
      d_pre = {
        yi <- ndvi_1000m; ti <- year
        ok <- is.finite(yi); yi <- yi[ok]; ti <- ti[ok]
        n  <- length(yi)
        if (n < 2) NA_real_ else {
          K  <- min(3L, floor(n/2))
          m1 <- mean(yi[seq_len(K)], na.rm = TRUE)
          m2 <- mean(yi[(n - K + 1):n], na.rm = TRUE)
          m2 - m1
        }
      },
      .groups = "drop"
    )
  
  placebo_vals <- placebo$d_pre[is.finite(placebo$d_pre)]
  noise_thr <- if (length(placebo_vals) >= 10)
    quantile(abs(placebo_vals), q_noise, na.rm = TRUE) else 0.02
  imp_thr   <- if (length(placebo_vals) >= 10)
    max(quantile(placebo_vals, q_imp, na.rm = TRUE), floor_improve) else floor_improve
  
  persist_margin <- imp_thr
  
  slope_vals <- ndvi_sum$pre_slope[is.finite(ndvi_sum$pre_slope)]
  sd_vals    <- ndvi_sum$pre_sd[is.finite(ndvi_sum$pre_sd)]
  slope_thr  <- if (length(slope_vals) >= 10)
    max(quantile(abs(slope_vals), slope_quant, na.rm = TRUE), slope_floor) else slope_floor
  sd_thr     <- if (length(sd_vals) >= 10)
    max(quantile(sd_vals, sd_quant, na.rm = TRUE), sd_floor) else sd_floor
  
  pers <- ndvi_long %>%
    left_join(ndvi_pre %>% select(person_id, pre_mean), by = "person_id") %>%
    filter(year %in% post_years) %>%
    mutate(
      above = ndvi_1000m >= (pre_mean + persist_margin),
      below = ndvi_1000m <= (pre_mean - persist_margin)
    ) %>%
    group_by(person_id) %>%
    summarise(
      post_n_pers      = n(),
      n_up             = sum(above, na.rm = TRUE),
      n_down           = sum(below, na.rm = TRUE),
      post_share_up    = ifelse(post_n_pers > 0, n_up   / post_n_pers, NA_real_),
      post_share_down  = ifelse(post_n_pers > 0, n_down / post_n_pers, NA_real_),
      required_up      = pmax(1L, ceiling(persistence_share * post_n_pers)),
      required_down    = pmax(1L, ceiling(persistence_share * post_n_pers)),
      .groups = "drop"
    )
  
  labels <- ndvi_sum %>%
    left_join(pers, by = "person_id", suffix = c("", ".pers")) %>%
    mutate(
      pre_stable = (pre_n >= min_pre) & (post_n >= min_post) &
        is.finite(pre_slope) & is.finite(pre_sd) &
        (abs(pre_slope) <= slope_thr) & (pre_sd <= sd_thr),
      
      improved = pre_stable & (delta >= imp_thr) & (n_up   >= required_up),
      declined = pre_stable & (delta <= -imp_thr) & (n_down >= required_down),
      stable   = pre_stable & (abs(delta) <= noise_thr),
      
      ndvi_group  = case_when(
        improved ~ "improved",
        stable   ~ "stable",
        TRUE     ~ "ambiguous_unstable"
      ),
      ndvi_group4 = case_when(
        improved ~ "improved",
        declined ~ "declined",
        stable   ~ "stable",
        TRUE     ~ "ambiguous_unstable"
      ),
      treated = ndvi_group == "improved"
    )
  
  if (any(labels$treated, na.rm = TRUE)) {
    mu_t <- mean(labels$pre_mean[labels$treated], na.rm = TRUE)
    sd_t <- sd(labels$pre_mean[labels$treated], na.rm = TRUE)
    if (is.finite(sd_t) && sd_t == 0) sd_t <- 1e-6
    labels <- labels %>%
      mutate(
        donor_common_support = (!treated) &
          is.finite(pre_mean) &
          abs(pre_mean - mu_t) <= donor_band_sd * sd_t
      )
  } else {
    labels$donor_common_support <- NA
  }
  
  thresholds <- list(
    years_pre         = pre_years,
    years_post        = post_years,
    min_ndvi          = min_ndvi,
    noise_threshold   = noise_thr,
    improvement_thr   = imp_thr,
    persistence_margin = persist_margin,
    slope_threshold   = slope_thr,
    sd_threshold      = sd_thr,
    min_pre           = min_pre,
    min_post          = min_post,
    persistence_req   = persistence_share,
    donor_band_sd     = donor_band_sd
  )
  
  diagnostics <- list(
    counts_3grp = table(labels$ndvi_group,  useNA = "ifany"),
    counts_4grp = table(labels$ndvi_group4, useNA = "ifany")
  )
  
  list(
    labels          = labels,
    thresholds      = thresholds,
    diagnostics     = diagnostics,
    ndvi_long       = ndvi_long,
    ndvi_wide_clean = nw
  )
}

res_1000 <- classify_ndvi_groups_v2_1000(ndvi_wide_1000, strict_clean = TRUE, min_ndvi = 0.05)
table(res_1000$labels$ndvi_group)

## ============================================================
## PART C: Combined PNAS-style plot (300 m & 1000 m)
## ============================================================
## ---------------------------------
## 300 m summary
## ---------------------------------
labels_300_plot <- res_300$labels %>%
  dplyr::mutate(
    group = dplyr::case_when(
      ndvi_group == "improved" ~ "Improved",
      ndvi_group == "stable"   ~ "Stable",
      TRUE                     ~ NA_character_
    )
  ) %>%
  dplyr::select(person_id, group)

ndvi_summary_300 <- res_300$ndvi_long %>%
  dplyr::left_join(labels_300_plot, by = "person_id") %>%
  dplyr::filter(group %in% c("Improved", "Stable")) %>%
  dplyr::group_by(year, group) %>%
  dplyr::summarise(
    mean_ndvi = mean(ndvi_300m, na.rm = TRUE),
    sd_ndvi   = sd(ndvi_300m,   na.rm = TRUE),
    n         = dplyr::n(),
    se_ndvi   = sd_ndvi / sqrt(n),
    ci_lower  = mean_ndvi - 1.96 * se_ndvi,
    ci_upper  = mean_ndvi + 1.96 * se_ndvi,
    .groups   = "drop"
  ) %>%
  dplyr::mutate(buffer = "300 m")

## ---------------------------------
## 1000 m summary
## ---------------------------------
labels_1000_plot <- res_1000$labels %>%
  dplyr::mutate(
    group = dplyr::case_when(
      ndvi_group == "improved" ~ "Improved",
      ndvi_group == "stable"   ~ "Stable",
      TRUE                     ~ NA_character_
    )
  ) %>%
  dplyr::select(person_id, group)

ndvi_summary_1000 <- res_1000$ndvi_long %>%
  dplyr::left_join(labels_1000_plot, by = "person_id") %>%
  dplyr::filter(group %in% c("Improved", "Stable")) %>%
  dplyr::group_by(year, group) %>%
  dplyr::summarise(
    mean_ndvi = mean(ndvi_1000m, na.rm = TRUE),
    sd_ndvi   = sd(ndvi_1000m,   na.rm = TRUE),
    n         = dplyr::n(),
    se_ndvi   = sd_ndvi / sqrt(n),
    ci_lower  = mean_ndvi - 1.96 * se_ndvi,
    ci_upper  = mean_ndvi + 1.96 * se_ndvi,
    .groups   = "drop"
  ) %>%
  dplyr::mutate(buffer = "1000 m")

## ---------------------------------
## Combine & plot (PNAS-style)
## ---------------------------------
ndvi_summary_both <- dplyr::bind_rows(ndvi_summary_300, ndvi_summary_1000) %>%
  dplyr::mutate(
    buffer = factor(buffer, levels = c("300 m", "1000 m")),
    group  = factor(group,  levels = c("Stable", "Improved"))
  )

y_max <- max(ndvi_summary_both$mean_ndvi, na.rm = TRUE)

p_ndvi <- ggplot(ndvi_summary_both,
                 aes(x = year, y = mean_ndvi,
                     color = group, fill = group)) +
  # 95% CI ribbon
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = 0.15, colour = NA, show.legend = FALSE) +
  # Lines + points
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  # Intervention line
  geom_vline(xintercept = 2020,
             linetype = "dashed",
             linewidth = 0.7,
             colour = "grey30") +
  annotate(
    "text",
    x = 2020 + 0.1,
    y = y_max * 0.9,
    label = "2020\n(intervention)",
    hjust = 0, vjust = 0,
    size = 4
  ) +
  # Facets for buffer size
  facet_wrap(~ buffer, ncol = 2) +
  # Elegant, colour-blind safe palette
  scale_color_manual(
    values = c(
      "Stable"   = "#7F7F7F",  # soft grey
      "Improved" = "#1B73B4"   # muted blue
    )
  ) +
  scale_fill_manual(
    values = c(
      "Stable"   = "#BFBFBF",
      "Improved" = "#9FC4E3"
    )
  ) +
  labs(
    x = "Year",
    y = "Mean NDVI",
    color = "NDVI group",
    fill  = "NDVI group",
    title = "Trajectories of Neighbourhood Greenness",
    subtitle = "Improved vs stable residential greenness at 300 m and 1000 m buffers"
  ) +
  scale_x_continuous(breaks = sort(unique(ndvi_summary_both$year))) +
  # PNAS-like theme: larger text
  theme_classic(base_size = 16) +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid         = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 15, face = "bold"),
    legend.position    = "top",
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 13),
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(size = 15),
    axis.text          = element_text(size = 13),
    axis.title         = element_text(size = 15)
  )

p_ndvi

# High-resolution TIFF output (PNAS style)
tiff("/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/fig/NDVI_trajectories_300_1000m.tiff",
     width = 10, height = 4.5, units = "in", res = 600)

y_max <- max(ndvi_summary_both$mean_ndvi, na.rm = TRUE)

print(
  ggplot(ndvi_summary_both,
         aes(x = year, y = mean_ndvi,
             color = group, fill = group)) +
    # 95% CI ribbon
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                alpha = 0.15, colour = NA, show.legend = FALSE) +
    # Lines + points
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.8) +
    # Intervention line
    geom_vline(xintercept = 2020, linetype = "dashed",
               linewidth = 0.7, colour = "grey30") +
    annotate("text",
             x = 2020 + 0.1,
             y = y_max * 0.9,
             label = "2020\n(intervention)",
             hjust = 0, vjust = 0,
             size = 4) +
    # Facets for buffer size
    facet_wrap(~ buffer, ncol = 2) +
    # Elegant, colour-blind safe palette
    scale_color_manual(
      values = c("Stable" = "#7F7F7F", "Improved" = "#1B73B4")
    ) +
    scale_fill_manual(
      values = c("Stable" = "#BFBFBF", "Improved" = "#9FC4E3")
    ) +
    labs(
      x = NULL,  # Year is obvious from axis ticks; keep figure cleaner
      y = "Mean NDVI",
      color = "NDVI group",
      fill  = "NDVI group",
      title = "Trajectories of Neighbourhood Greenness",
      subtitle = "Improved vs stable residential greenness at 300 m and 1000 m buffers"
    ) +
    scale_x_continuous(breaks = sort(unique(ndvi_summary_both$year))) +
    # PNAS-like theme: slightly larger, readable text
    theme_classic(base_size = 16) +
    theme(
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid         = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(size = 15, face = "bold"),
      # Legend closer to plot (reduce extra gap above panels)
      legend.position    = "top",
      legend.title       = element_text(size = 14),
      legend.text        = element_text(size = 13),
      legend.box.margin  = margin(b = -5, t = 0, l = 0, r = 0),
      legend.spacing.y   = unit(0.1, "lines"),
      # Titles and axes
      plot.title         = element_text(face = "bold", size = 20),
      plot.subtitle      = element_text(size = 14),
      axis.text          = element_text(size = 13),
      axis.title.y       = element_text(size = 15),
      axis.title.x       = element_blank(),
      plot.margin        = margin(t = 5, r = 10, b = 5, l = 10)
    )
)

dev.off()

# Period-specific mean NDVI (2016–2020 and 2021–2024) by buffer & group
ndvi_period_lines <- ndvi_summary_both %>%
  dplyr::filter(year >= 2016 & year <= 2024) %>%
  dplyr::mutate(
    period = dplyr::case_when(
      year <= 2020 ~ "2016–2020",
      year >= 2021 ~ "2021–2024"
    )
  ) %>%
  dplyr::group_by(buffer, group, period) %>%
  dplyr::summarise(
    mean_period = mean(mean_ndvi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    x_start = dplyr::case_when(
      period == "2016–2020" ~ 2016,
      period == "2021–2024" ~ 2021
    ),
    x_end   = dplyr::case_when(
      period == "2016–2020" ~ 2020,
      period == "2021–2024" ~ 2024
    )
  )

y_max <- max(ndvi_summary_both$mean_ndvi, na.rm = TRUE)

p_ndvi <- ggplot(ndvi_summary_both,
                 aes(x = year, y = mean_ndvi,
                     color = group, fill = group)) +
  # 95% CI ribbon
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = 0.15, colour = NA, show.legend = FALSE) +
  # Lines + points
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  # NEW: horizontal dashed lines for period means (per buffer × group)
  geom_segment(
    data = ndvi_period_lines,
    aes(x = x_start, xend = x_end,
        y = mean_period, yend = mean_period,
        color = group),
    linetype = "dashed",
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  # Intervention line
  geom_vline(xintercept = 2020,
             linetype = "dashed",
             linewidth = 0.7,
             colour = "grey30") +
  annotate(
    "text",
    x = 2020 + 0.1,
    y = y_max * 0.9,
    label = "2020\n(intervention)",
    hjust = 0, vjust = 0,
    size = 4
  ) +
  facet_wrap(~ buffer, ncol = 2) +
  scale_color_manual(
    values = c(
      "Stable"   = "#7F7F7F",
      "Improved" = "#1B73B4"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Stable"   = "#BFBFBF",
      "Improved" = "#9FC4E3"
    )
  ) +
  labs(
    x = "Year",
    y = "Mean NDVI",
    color = "NDVI group",
    fill  = "NDVI group",
    title = "Trajectories of Neighbourhood Greenness",
    subtitle = "Improved vs stable residential greenness at 300 m and 1000 m buffers"
  ) +
  scale_x_continuous(breaks = sort(unique(ndvi_summary_both$year))) +
  theme_classic(base_size = 16) +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid         = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 15, face = "bold"),
    legend.position    = "top",
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 13),
    plot.title         = element_text(face = "bold", size = 20),
    plot.subtitle      = element_text(size = 15),
    axis.text          = element_text(size = 13),
    axis.title         = element_text(size = 15)
  )

p_ndvi

tiff("/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/fig/NDVI_trajectories_300_1000m.tiff",
     width = 10, height = 4.5, units = "in", res = 600)

y_max <- max(ndvi_summary_both$mean_ndvi, na.rm = TRUE)

print(
  ggplot(ndvi_summary_both,
         aes(x = year, y = mean_ndvi,
             color = group, fill = group)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
                alpha = 0.15, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.8) +
    # NEW: horizontal dashed lines for period means
    geom_segment(
      data = ndvi_period_lines,
      aes(x = x_start, xend = x_end,
          y = mean_period, yend = mean_period,
          color = group),
      linetype = "dashed",
      linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 2020, linetype = "dashed",
               linewidth = 0.7, colour = "grey30") +
    annotate("text",
             x = 2020 + 0.1,
             y = y_max * 0.9,
             label = "2020\n(intervention)",
             hjust = 0, vjust = 0,
             size = 4) +
    facet_wrap(~ buffer, ncol = 2) +
    scale_color_manual(
      values = c("Stable" = "#7F7F7F", "Improved" = "#1B73B4")
    ) +
    scale_fill_manual(
      values = c("Stable" = "#BFBFBF", "Improved" = "#9FC4E3")
    ) +
    labs(
      x = NULL,
      y = "Mean NDVI",
      color = "NDVI group",
      fill  = "NDVI group",
      title = "Trajectories of Neighbourhood Greenness",
      subtitle = "Improved vs stable residential greenness at 300 m and 1000 m buffers"
    ) +
    scale_x_continuous(breaks = sort(unique(ndvi_summary_both$year))) +
    theme_classic(base_size = 16) +
    theme(
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid         = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(size = 15, face = "bold"),
      legend.position    = "top",
      legend.title       = element_text(size = 14),
      legend.text        = element_text(size = 13),
      legend.box.margin  = margin(b = -5, t = 0, l = 0, r = 0),
      legend.spacing.y   = unit(0.1, "lines"),
      plot.title         = element_text(face = "bold", size = 20),
      plot.subtitle      = element_text(size = 14),
      axis.text          = element_text(size = 13),
      axis.title.y       = element_text(size = 15),
      axis.title.x       = element_blank(),
      plot.margin        = margin(t = 5, r = 10, b = 5, l = 10)
    )
)

dev.off()

