############################################################
# Correlations between CP frequency and climate indices (Spearman)
#
# This script:
# 1) Reads global climate indices from an Excel file (one sheet per season).
# 2) Reads daily CP classification (K-means / PCA+kmeans / UV-based variants).
# 3) Computes CP frequencies by season and by hydrological year (año_hidro).
#    - DJF is treated using hydrological years (July–June) to avoid splitting DJF
#      across calendar years.
#    - MAM/JJA/SON use calendar year.
#    - ANNUAL uses hydrological year.
# 4) Detrends CP-frequency series and climate-index series using linear detrending:
#      residuals(lm(x ~ seq_along(x)))
# 5) Computes Spearman correlation matrices (CPs vs indices) and p-values.
# 6) Produces correlation heatmaps, labeling only significant cells (p <= 0.05).
#
# Inputs (required):
# - Excel file with indices:
#   "modos_variabilidad_global.xlsx"
#   Sheets must be named "Todas_<SEASON>" where SEASON in:
#     DJF, MAM, JJA, SON, ANNUAL
#   Must include a column 'año' plus index columns (e.g., SAM, ENSO, PDO, etc.)
# - CP classification CSV:
#   dias_clasificados_<tecnica>.csv
#   Must include columns: año, mes, dia, CP
#
# User-defined parameters (must be set before running):
# - año1, año2: start/end year for analysis (inclusive)
#
# Outputs:
# - Example time-series plots saved as "series.png" (overwritten each run)
# - Correlation heatmap saved as "cor_CPs_SIN_LAG_ALT.png"
#
############################################################

rm(list=ls()) ; graphics.off() ; gc()

# -------------------------
# Libraries
# -------------------------
library(readxl)   # read Excel indices
library(ggplot2)  # plotting
library(reshape2) # melt for matrices
library(corrplot)
library(writexl)
library(tidyr)
library(dplyr)
library(trend)
library(gridExtra)
library(scales)
library(zoo)

# -------------------------
# User parameters (SET THESE)
# -------------------------
año1 <- 1979
año2 <- 2023

# Paths
ruta_excel <- "C:/Users/Franco/Documents/Doctorado/seasonal_precipitation_NWA/modos_variabilidad_global.xlsx"
ruta_cps   <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures/dias_clasificados_UV_PCA_k-means7.csv")
out_dir    <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures"

# -------------------------
# Helper functions
# -------------------------

# Assign season from month (DJF/MAM/JJA/SON)
month_to_season <- function(mes) {
  out <- rep("DJF", length(mes))
  out[mes %in% c(3,4,5)]   <- "MAM"
  out[mes %in% c(6,7,8)]   <- "JJA"
  out[mes %in% c(9,10,11)] <- "SON"
  out
}

# Hydrological year definition used here:
# año_hidro = año, but if month >= 7 then año_hidro = año + 1
# (i.e., hydro-year runs Jul–Jun and is labeled by the ending year)
calc_año_hidro <- function(año, mes) {
  ah <- año
  ah[mes >= 7] <- ah[mes >= 7] + 1
  ah
}

# Detrend a CP-frequency dataset in long format:
# expects columns: año, season, CP, N
# returns a numeric matrix (rows = years, cols = CPs) detrended per CP column
matriz_detrend <- function(x, año1, año2) {
  # long -> wide (columns = CPs)
  CPs_wide <- x %>%
    pivot_wider(names_from = CP, values_from = N)
  
  # replace missing CP counts with 0
  CPs_wide[is.na(CPs_wide)] <- 0
  
  # keep requested time window
  CPs_wide <- CPs_wide[CPs_wide$año %in% año1:año2, ]
  
  # drop identifiers (year/season) and detrend each CP column
  CPs_mat <- CPs_wide %>% select(-año, -season)
  
  CPs_mat_detr <- apply(CPs_mat, 2, function(v) residuals(lm(v ~ seq_along(v))))
  
  return(CPs_mat_detr)
}

# Read climate indices for a given season; detrend each index column
read_and_detrend_indices <- function(ruta_excel, season, año1) {
  sheet_name <- paste0("Todas_", season)
  df <- read_excel(ruta_excel, sheet = sheet_name)
  
  # keep years >= año1
  df <- df[df$año >= año1, ]
  
  # keep only index columns (drop year)
  idx <- df[, -1, drop = FALSE]
  
  # detrend each index series
  idx_detr <- apply(idx, 2, function(v) residuals(lm(v ~ seq_along(v))))
  
  colnames(idx_detr) <- colnames(idx)
  return(list(years = df$año, idx_detr = idx_detr))
}

# Spearman correlation matrix + p-value matrix (CPs vs indices)
spearman_cor_and_p <- function(mat1, mat2) {
  cor_mat <- matrix(NA_real_, nrow = ncol(mat1), ncol = ncol(mat2))
  p_mat   <- matrix(NA_real_, nrow = ncol(mat1), ncol = ncol(mat2))
  
  for (i in 1:ncol(mat1)) {
    for (j in 1:ncol(mat2)) {
      test <- cor.test(mat1[, i], mat2[, j], method = "spearman", use = "complete.obs")
      cor_mat[i, j] <- as.numeric(test$estimate)
      p_mat[i, j]   <- test$p.value
    }
  }
  
  rownames(cor_mat) <- colnames(mat1)
  colnames(cor_mat) <- colnames(mat2)
  rownames(p_mat)   <- colnames(mat1)
  colnames(p_mat)   <- colnames(mat2)
  
  return(list(cor = cor_mat, p = p_mat))
}

#####################################################################################################
# 1) Read CP daily classification and compute frequencies
#####################################################################################################

CPs <- read.csv(ruta_cps)

# Add season + hydro-year
CPs$season   <- month_to_season(CPs$mes)
CPs$año_hidro <- calc_año_hidro(CPs$año, CPs$mes)

# DJF: use hydro-year (to keep DJF together)
cps_frec_hidro <- CPs[CPs$season == "DJF", ] %>%
  group_by(season, año_hidro, CP) %>%
  summarise(N = n(), .groups = "drop")

# ANNUAL: use hydro-year
cps_frec_hidro_anual <- CPs %>%
  group_by(año_hidro, CP) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(season = "ANNUAL")

# MAM/JJA/SON: calendar-year frequencies (DJF excluded)
cps_frec <- CPs[CPs$season != "DJF", ] %>%
  group_by(season, año, CP) %>%
  summarise(N = n(), .groups = "drop")

# Rename hydro-year column to 'año' for compatibility downstream
names(cps_frec_hidro)[names(cps_frec_hidro) == "año_hidro"] <- "año"
names(cps_frec_hidro_anual)[names(cps_frec_hidro_anual) == "año_hidro"] <- "año"

#####################################################################################################
# 2) Restrict years to those available in the climate indices file
#####################################################################################################

ref_idx <- read_excel(ruta_excel, sheet = "Todas_ANNUAL")
ref_years <- ref_idx$año[ref_idx$año >= año1]

cps_frec       <- cps_frec[cps_frec$año %in% ref_years, ]
cps_frec_hidro <- cps_frec_hidro[cps_frec_hidro$año %in% ref_years, ]
cps_frec_hidro_anual <- cps_frec_hidro_anual[cps_frec_hidro_anual$año %in% ref_years, ]

#####################################################################################################
# 3) Detrend CP-frequency matrices (one per season + annual)
#####################################################################################################

CPs_DJF    <- matriz_detrend(cps_frec_hidro, año1, año2)
CPs_MAM    <- matriz_detrend(cps_frec[cps_frec$season == "MAM", ], año1, año2)
CPs_JJA    <- matriz_detrend(cps_frec[cps_frec$season == "JJA", ], año1, año2)
CPs_SON    <- matriz_detrend(cps_frec[cps_frec$season == "SON", ], año1, año2)
CPs_ANNUAL <- matriz_detrend(cps_frec_hidro_anual, año1, año2)

#####################################################################################################
# 4) Example: plot CP vs one index (raw + moving mean + differenced moving mean)
#####################################################################################################
# This section is optional and mainly for a “figure-ready” time-series comparison.

season <- "MAM"
CP_sel <- 3
index_name <- "SAM"                      # must be a column name in the corresponding sheet
title_name <- paste0("CP", CP_sel, " and ", index_name, " (", season, ")")

# CP series for selected season and CP (calendar year)
s <- cps_frec %>% filter(season == season, CP == CP_sel) %>% arrange(año)

idx_obj <- read_and_detrend_indices(ruta_excel, season, año1)
df_ind <- read_excel(ruta_excel, sheet = paste0("Todas_", season))
df_ind <- df_ind[df_ind$año >= año1, ]

# Build plotting dataframe (match years explicitly)
years_use <- año1:año2
df_plot <- data.frame(
  año    = years_use,
  CP     = s$N[s$año %in% years_use],
  indice = df_ind[[index_name]][df_ind$año %in% years_use]
)

# 9-yr moving mean
df_plot$indice_suav <- rollmean(df_plot$indice, k = 9, fill = NA)
df_plot$CP_suav     <- rollmean(df_plot$CP,     k = 9, fill = NA)

# "detrended" here is actually first-difference of the smoothed series
df_plot$indice_suav_detrended <- c(NA, diff(df_plot$indice_suav))
df_plot$CP_suav_detrended     <- c(NA, diff(df_plot$CP_suav))

# Dual-axis scaling parameters (purely for visualization)
signo <- 1
factor_transform <- signo * 10
offset_transform <- 25

RHO <- cor.test(df_plot$CP, df_plot$indice, method = "spearman")$estimate * signo
RHO <- round(RHO, 2)

g1 <- ggplot(df_plot) +
  geom_line(aes(x = año, y = CP), linewidth = 1) +
  geom_line(aes(x = año, y = indice * factor_transform + offset_transform), color = "red", linewidth = 1) +
  geom_smooth(method = "lm", aes(x = año, y = indice * factor_transform + offset_transform),
              color = "red", se = FALSE, linewidth = 0.5, linetype = "dashed") +
  geom_smooth(method = "lm", aes(x = año, y = CP),
              color = "black", se = FALSE, linewidth = 0.5, linetype = "dashed") +
  scale_y_continuous(
    name = "N (days)",
    sec.axis = sec_axis(~ (. - offset_transform) / factor_transform, name = "Index")
  ) +
  theme_bw() +
  xlab("year") +
  ggtitle(title_name)
# + annotate("text", x = max(df_plot$año, na.rm=TRUE), y = offset_transform, label = paste("RHO =", RHO))

plot(g1)

setwd(out_dir)
write.csv(df_plot, "series.csv", row.names = FALSE)
ggsave(plot = g1, filename = "series.png", height = 2.5, width = 6)

#####################################################################################################
# 5) Correlation heatmap (ALT): stack seasons and facet by climate index
#####################################################################################################
# It returns a tidy dataframe with: Var1 (CP), Var2 (Index), value (rho), p, season

figura_tidy <- function(season, CPs_mat) {
  
  idx_obj <- read_and_detrend_indices(ruta_excel, season, año1)
  detrended_series <- idx_obj$idx_detr
  
  # Compute Spearman rho + p-values
  res <- spearman_cor_and_p(CPs_mat, detrended_series)
  
  # Convert to long tidy format
  cor_long <- melt(res$cor, varnames = c("Var1", "Var2"), value.name = "value")
  p_long   <- melt(res$p,   varnames = c("Var1", "Var2"), value.name = "p")
  
  df <- merge(cor_long, p_long, by = c("Var1", "Var2"))
  df$season <- season
  
  return(df)
}

df_son    <- figura_tidy("SON",    CPs_SON)
df_djf    <- figura_tidy("DJF",    CPs_DJF)
df_mam    <- figura_tidy("MAM",    CPs_MAM)
df_jja    <- figura_tidy("JJA",    CPs_JJA)
df_annual <- figura_tidy("ANNUAL", CPs_ANNUAL)

df_plot <- bind_rows(df_annual, df_djf, df_mam, df_jja, df_son)

# Formatting
df_plot$season <- factor(df_plot$season, levels = c("ANNUAL", "DJF", "MAM", "JJA", "SON"))
df_plot$value  <- round(df_plot$value, 1)

# Plot option A: y = CPs (Var1), facet by index (Var2), x = season
# Labels only where p <= 0.05
g <- ggplot(df_plot, aes(x = season, y = Var1, fill = value)) +
  geom_tile(color = "grey", linewidth = 0.3) +
  geom_text(
    data = df_plot %>% filter(p <= 0.05),
    aes(label = value),
    color = "black", size = 2
  ) +
  scale_fill_gradient2(
    low = "deepskyblue", high = "brown1", mid = "white", midpoint = 0,
    na.value = "white", limits = c(-0.6, 0.6),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    oob = squish, name = "cor"
  ) +
  facet_wrap(~ Var2, scales = "free", nrow = 2) +
  theme_minimal(base_size = 7) +
  labs(title = "", x = "", y = "CP") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right",
    legend.key.width  = unit(0.2, "cm"),
    legend.key.height = unit(0.4, "cm")
  )

plot(g)

setwd(out_dir)
ggsave(plot = g, filename = "cor_CPs_SIN_LAG_ALT.png", height = 3, width = 5, dpi = 300, bg = "white")
