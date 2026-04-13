## =========================================================
## China (with province borders + Taiwan) and Tianjin highlighted
## For Nature-style figure panel
## =========================================================

rm(list = ls())

# ---- 1. Packages ----
library(geodata)   # for gadm() online download
library(sf)
library(dplyr)
library(ggplot2)

# ---- 2. Download level-1 admin units (provinces) ----
gadm_dir <- tempdir()  # or set to a permanent folder

# China provinces / municipalities
chn1  <- gadm(country = "China",  level = 1,
              path = gadm_dir, version = "latest")

# Taiwan counties/cities (level 1 so that internal borders show)
twn1  <- gadm(country = "Taiwan", level = 1,
              path = gadm_dir, version = "latest")

chn1_sf <- st_as_sf(chn1)
twn1_sf <- st_as_sf(twn1)

# ---- 3. Combine China + Taiwan at province level ----
china_prov_all <- bind_rows(
  chn1_sf |> mutate(country = "China"),
  twn1_sf |> mutate(country = "Taiwan")
)

# ---- 4. Extract Tianjin polygon from China provinces ----
# You can check names if needed:
# sort(unique(chn1_sf$NAME_1))
tianjin_sf <- chn1_sf |>
  filter(NAME_1 %in% c("Tianjin", "Tianjin Shi"))

# ---- 5. CRS: keep WGS84 (4326) for now ----
china_prov_all <- st_transform(china_prov_all, 4326)
tianjin_sf     <- st_transform(tianjin_sf, 4326)

# ---- 6. Plot: province borders + Tianjin highlight ----
p_china_tj <- ggplot() +
  # background: all provinces (China + Taiwan)
  geom_sf(data = china_prov_all,
          fill      = "grey92",
          color     = "grey60",   # province boundaries visible
          linewidth = 0.25) +
  # highlight: Tianjin municipality
  geom_sf(data = tianjin_sf,
          fill      = "#d73027",  # highlight colour
          color     = "black",
          linewidth = 0.35) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.title       = element_blank(),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.margin      = margin(3, 3, 3, 3)
  ) +
  ggtitle("Tianjin in China")

print(p_china_tj)

# ---- 7. Save with a "proper" Nature-style size ----
# Nature single-column width ≈ 8.9 cm; adjust height slightly
out_width  <- 8.9   # cm
out_height <- 7.0   # cm

ggsave("china_tianjin_province_borders.pdf", p_china_tj,
       width = out_width, height = out_height,
       units = "cm", dpi = 600)

ggsave("china_tianjin_province_borders.tiff", p_china_tj,
       width = out_width, height = out_height,
       units = "cm", dpi = 600, compression = "lzw", bg = "white")
