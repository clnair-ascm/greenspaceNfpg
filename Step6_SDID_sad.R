# ============================================================
# Clear workspace
rm(list = ls())

# ---- packages ----
library(openair)
library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(ggplot2)
library(leaflet)

# ============================================================
# 1. Define sites and common settings
# ============================================================
years      <- 2019:2024
pollutants <- c("pm2.5", "pm10", "o3", "no2", "nox")

site_codes <- c(
  Ladywood                    = "BMLD",
  A4540_Roadside              = "BIRR",
  West_Bromwich               = "WBKP",
  Oldbury_Birmingham_Road     = "BOLD",
  Coventry_Allesley           = "COAL",
  Coventry_Binley_Road        = "COBR"
)

# Safe wrapper: if a site/year has no data, return NULL instead of stopping
safe_importAURN <- purrr::possibly(importAURN, otherwise = NULL)

# ============================================================
# 2. Download AURN data for all sites
# ============================================================

aurn_list <- imap(site_codes, ~{
  message("Downloading ", .y, " (", .x, ") ...")
  safe_importAURN(
    site      = .x,
    year      = years,
    pollutant = pollutants,
    meta      = TRUE
  )
})

# Drop sites that returned NULL (no data)
aurn_list <- aurn_list[!vapply(aurn_list, is.null, logical(1))]

# Combine into one data frame with site_name column
aurn_all <- bind_rows(aurn_list, .id = "site_name")

# Quick check
dplyr::glimpse(aurn_all)

# ============================================================
# 3. Long-format data for ggplot (pm25, pm10, o3, no2, nox)
# ============================================================

aurn_long <- aurn_all %>%
  select(site_name, date, all_of(pollutants)) %>%
  pivot_longer(
    cols      = all_of(pollutants),
    names_to  = "pollutant",
    values_to = "value"
  )

# --- Optional: filter out completely missing rows ---
aurn_long <- aurn_long %>% filter(!is.na(value))

pollutant_order <- c("pm25", "pm10", "no2", "nox", "o3")

weekly_long <- aurn_long %>%
  filter(pollutant %in% pollutant_order) %>%
  mutate(
    week      = as.Date(floor_date(date, unit = "week", week_start = 1)),  # <- Date
    pollutant = factor(pollutant, levels = pollutant_order)
  ) %>%
  group_by(site_name, pollutant, week) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    value = ifelse(is.nan(value), NA_real_, value)
  )

pollutant_labels <- c(
  pm25 = "PM2.5 (µg/m³)",
  pm10 = "PM10 (µg/m³)",
  no2  = "NO2 (µg/m³)",
  nox  = "NOₓ (µg/m³)",
  o3   = "O₃ (µg/m³)"
)

gg_weekly <- ggplot(weekly_long,
                    aes(x = week, y = value, colour = site_name)) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  facet_wrap(
    ~ pollutant,
    scales   = "free_y",
    ncol     = 1,
    labeller = as_labeller(pollutant_labels)
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  labs(
    title  = "Weekly mean concentrations at selected AURN sites (2019–2024)",
    x      = NULL,
    y      = "Weekly mean concentration",
    colour = "Site"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position      = "bottom",
    legend.title         = element_text(face = "bold"),
    panel.grid.minor     = element_blank(),
    strip.background     = element_rect(fill = "grey90", colour = NA),
    strip.text           = element_text(face = "bold"),
    plot.title           = element_text(face = "bold", hjust = 0),
    axis.title.y         = element_text(margin = margin(r = 5))
  )

print(gg_weekly)
ggsave("/Users/admin-macpro/Downloads/TimeV.jpeg", gg_weekly, width = 12, height = 8, dpi = 300)

# ============================================================
# 5. Openair timeVariation example (all sites together)
# ============================================================

# This gives diurnal/weekday/annual patterns for selected pollutants
# timeVariation(
#   aurn_all,
#   pollutant = c("pm2.5", "pm10", "no2", "o3"),
#   type      = "site_name"   # separate panels by site
# )

# ============================================================
# 6. Map of site locations
# ============================================================
library(sf)
library(leaflet)
library(dplyr)

# ------------------------------------------------------------
# 1. Read West Midlands CA wards boundary and reproject
# ------------------------------------------------------------

wm_wards <- st_read(
  "/Users/admin-macpro/Downloads/cn_shp/boundaries-wards-2023-wmca/boundaries-wards-2023-wmca.shp"
)

# Make sure it's in WGS84 (lat/lon) for leaflet
wm_wards <- st_transform(wm_wards, 4326)

# (Optional) make a single West Midlands polygon if you prefer outline only
wm_boundary <- wm_wards %>%
  summarise(geometry = st_union(geometry))

# ------------------------------------------------------------
# 2. AURN site metadata (as you had)
# ------------------------------------------------------------

meta_aurn <- importMeta(source = "aurn")

meta_sites <- meta_aurn %>%
  filter(code %in% site_codes) %>%
  left_join(
    tibble(site_name = names(site_codes), code = as.character(site_codes)),
    by = "code"
  )

meta_sites %>% select(site_name, site, code, latitude, longitude)

# ------------------------------------------------------------
# 3. Leaflet map: West Midlands boundary + AURN sites
# ------------------------------------------------------------

library(leaflet)

leaflet() %>%
  addTiles() %>%
  
  # ---- West Midlands boundary ----
addPolygons(
  data        = wm_boundary,
  fill        = FALSE,
  weight      = 2,
  color       = "black",
  opacity     = 0.8,
  label       = ~"West Midlands CA"
) %>%
  
  # ---- Circle markers with hover labels + popups ----
addCircleMarkers(
  data   = meta_sites,
  lng    = ~longitude,
  lat    = ~latitude,
  radius = 6,
  stroke = TRUE,
  fillOpacity = 0.9,
  label = ~paste0(site_name, " (", site_type, ")"),  # hover label
  popup = ~paste0(
    "<b>", site_name, "</b><br>",
    site, " (", code, ")<br>",
    "Type: ", site_type, "<br>",
    "Lat: ", round(latitude, 4), ", Lon: ", round(longitude, 4)
  ),
  labelOptions = labelOptions(
    direction = "auto",
    textsize  = "12px"
  )
)


library(sf)
library(dplyr)
library(ggplot2)
library(ggrepel)

# -----------------------------------------------------------
# 1. Read West Midlands CA boundary (your local shapefile)
# -----------------------------------------------------------
wm_wards <- st_read(
  "/Users/admin-macpro/Downloads/cn_shp/boundaries-wards-2023-wmca/boundaries-wards-2023-wmca.shp"
)

# Combine all wards into one region boundary for a clean outline
wm_boundary <- wm_wards %>%
  summarise(geometry = st_union(geometry))

# -----------------------------------------------------------
# 2. Convert AURN site metadata to sf
# -----------------------------------------------------------
meta_sites_sf <- st_as_sf(
  meta_sites,
  coords = c("longitude", "latitude"),
  crs = 4326
)

# -----------------------------------------------------------
# 3. Static map using ggplot2 + ggrepel
# -----------------------------------------------------------

gg_map <- ggplot() +
  
  # West Midlands outline
  geom_sf(
    data = wm_boundary,
    fill = "grey95",
    color = "black",
    linewidth = 0.7
  ) +
  
  # AURN site points
  geom_sf(
    data = meta_sites_sf,
    aes(color = site_type),
    size = 3
  ) +
  
  # Smart non-overlapping labels
  geom_text_repel(
    data = meta_sites_sf,
    aes(
      geometry = geometry,
      label = paste0(site_name, "\n(", site_type, ")")
    ),
    stat = "sf_coordinates",
    size = 3.5,
    nudge_x = 0.03,
    nudge_y = 0.03,
    min.segment.length = 0,
    segment.size = 0.4
  ) +
  
  scale_color_brewer(palette = "Dark2") +
  
  labs(
    title = "AURN Monitoring Sites in the West Midlands",
    color = "Site Type"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    plot.title = element_text(face = "bold", size = 14)
  )

# Show the map
print(gg_map)

# You can export to PPT using ggsave:
ggsave("/Users/admin-macpro/Downloads/AURN_WestMidlands_Map.jpeg", gg_map, width = 14, height = 10, dpi = 300)
