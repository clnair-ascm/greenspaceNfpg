# ============================================================
# Keep only selected variables in each Excel file (drop others)
# ============================================================

library(readxl)
library(writexl)
library(dplyr)

# --- 1. Directory containing your Excel files ---
dir_path <- "/Volumes/shiz-wm-netzero/users/yuqing/daiyy/GreenSpace/treated/data/data_14T24NO23"

# --- 2. List of columns to keep (case-insensitive) ---
cols_keep <- c(
  "age","waist","hemoglobin","fasting_glucose","creatinine",
  "cholesterol_total","hypertension","diabetes","sex","smoking","bmi",
  "person_id","year","lat_gcj","lon_gcj","pulse_pressure","map",
  "waist_u","hb_l","hb_u","cr_u","fpg_sev","waist_sev","hb_sev","cr_sev",
  "tc_sev","bmi_sev","pp_sev","map_sev","fpg_dto","waist_t","waist_dto",
  "hb_t","hb_dto","cr_t","cr_dto","tc_dto","bmi_dto","pp_dto","map_dto"
)

# --- 3. List Excel files (recursively if needed) ---
files <- list.files(dir_path, pattern = "\\.xlsx$", full.names = TRUE)
message("Found ", length(files), " Excel file(s).")

# --- 4. Loop through each file ---
for (f in files) {
  message("Processing: ", basename(f))
  
  df <- tryCatch(readxl::read_excel(f), error = function(e) {
    message("  ⚠️ Failed to read: ", conditionMessage(e)); return(NULL)
  })
  if (is.null(df)) next
  
  # make names lowercase and trim spaces
  names(df) <- tolower(trimws(names(df)))
  
  # check which columns exist
  keep_cols <- intersect(cols_keep, names(df))
  missing_cols <- setdiff(cols_keep, names(df))
  
  if (length(keep_cols) == 0) {
    message("  ⚠️ No matching columns in ", basename(f), ", skipping.")
    next
  }
  
  if (length(missing_cols) > 0) {
    message("  ⚠️ Missing ", length(missing_cols), " expected columns (skipped): ",
            paste(missing_cols, collapse = ", "))
  }
  
  # subset dataframe
  df_clean <- df %>% select(all_of(keep_cols))
  
  # overwrite file (or create _clean version)
  writexl::write_xlsx(df_clean, f)
  message("  ✅ Saved cleaned file: ", basename(f), "\n")
}

message("🎯 All Excel files processed successfully.")
