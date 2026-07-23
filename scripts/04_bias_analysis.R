R.version.string
install.packages(c("dplyr", "readr", "ggplot2"))

# 04_bias_analysis.R
#
# Computes bias, RMSE, and correlation between observed and modelled SWE
# for each station, joins with station metadata (lat/lon/elevation), and
# makes a quick diagnostic plot of bias vs. elevation.
#
# Usage (from RStudio, with working directory set to the repo root):
#   source("scripts/04_bias_analysis.R")
#
# Inputs:
#   results/tables/modelled_vs_observed_swe.csv.gz
#   results/tables/station_metadata.csv
#
# Outputs:
#   results/tables/station_bias_metrics.csv
#   results/figures/bias_vs_elevation.png

library(dplyr)
library(readr)
library(ggplot2)

# --- Paths -------------------------------------------------------------
# Assumes RStudio's working directory is the repo root (aleyna-swe-thesis).
# If not, run setwd() first to point there.

combined_path <- "results/tables/modelled_vs_observed_swe.csv.gz"
metadata_path <- "results/tables/station_metadata.csv"
out_metrics_path <- "results/tables/station_bias_metrics.csv"
out_fig_path <- "results/figures/bias_vs_elevation.png"

# --- Load data -----------------------------------------------------------
cat("Loading combined observed/modelled data...\n")
combined <- read_csv(combined_path, show_col_types = FALSE)
cat("Rows loaded:", nrow(combined), "\n")

cat("Loading station metadata...\n")
metadata <- read_csv(metadata_path, show_col_types = FALSE)
cat("Stations in metadata:", nrow(metadata), "\n")

# --- Compute per-station bias metrics ------------------------------------
cat("\nComputing bias metrics per station...\n")

combined_clean <- combined %>%
  filter(!is.na(observed_swe_mm), !is.na(modelled_swe_mm))

station_metrics <- combined_clean %>%
  group_by(station_id) %>%
  summarise(
    n_obs = n(),
    bias_mm = mean(modelled_swe_mm - observed_swe_mm),
    rmse_mm = sqrt(mean((modelled_swe_mm - observed_swe_mm)^2)),
    correlation = if (n() > 1) cor(modelled_swe_mm, observed_swe_mm) else NA_real_,
    mean_observed_mm = mean(observed_swe_mm),
    mean_modelled_mm = mean(modelled_swe_mm),
    .groups = "drop"
  )

cat("Stations with computed metrics:", nrow(station_metrics), "\n")

# --- Join with metadata (lat/lon/elevation/name) -------------------------
station_metrics <- station_metrics %>%
  left_join(
    metadata %>% select(station_id, station_name, latitude, longitude, elevation_m),
    by = "station_id"
  )

# --- Save the result -------------------------------------------------------
dir.create(dirname(out_metrics_path), showWarnings = FALSE, recursive = TRUE)
write_csv(station_metrics, out_metrics_path)
cat("\nSaved station bias metrics to:", out_metrics_path, "\n")

# --- Quick summary printed to console ------------------------------------
cat("\n--- Overall summary across all stations ---\n")
cat("Mean bias (mm):    ", round(mean(station_metrics$bias_mm, na.rm = TRUE), 2), "\n")
cat("Mean RMSE (mm):    ", round(mean(station_metrics$rmse_mm, na.rm = TRUE), 2), "\n")
cat("Mean correlation:  ", round(mean(station_metrics$correlation, na.rm = TRUE), 3), "\n")

# --- Diagnostic plot: bias vs elevation -----------------------------------
cat("\nCreating bias vs elevation plot...\n")

dir.create(dirname(out_fig_path), showWarnings = FALSE, recursive = TRUE)

p <- ggplot(station_metrics, aes(x = elevation_m, y = bias_mm)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, color = "darkorange") +
  labs(
    title = "mHM SWE bias vs. station elevation",
    subtitle = "Bias = mean(modelled - observed); positive = model overpredicts",
    x = "Elevation (m)",
    y = "Bias (mm)"
  ) +
  theme_minimal()

ggsave(out_fig_path, plot = p, width = 8, height = 6, dpi = 150)
cat("Saved plot to:", out_fig_path, "\n")

cat("\nDone.\n")

source("scripts/04_bias_analysis.R")
