# =============================================================================
# SACZ definition (Vera et al. 2010) and association with Circulation Patterns
# Reference: Vera et al. (2010), Climate Dynamics, https://doi.org/10.1007/s00382-010-0812-4
# =============================================================================

# ----------------------------
# Libraries (keep only what is used)
# ----------------------------
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ----------------------------
# User parameters (EDIT)
# ----------------------------
k <- 7
season_name <- "DJF"   # "ONDJFMA", "DJF", "MAM", "JJA", "SON"

base_dir <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st"
fig_dir  <- file.path(base_dir, "Figures")

# Period blocks available as separate NetCDF files
periods <- c("1979-1989", "1990-1999", "2000-2009", "2010-2023")

# OLR file template
olr_file <- function(periodo) file.path(base_dir, paste0("OLR_NCEP1_Vera_", periodo, ".nc"))

# CP classification file
cp_file <- file.path(fig_dir, paste0("dias_clasificados_UV_PCA_k-means", k, ".csv"))

# Monte Carlo
n_mc <- 1000
set.seed(1234)

# ----------------------------
# 0) Season definitions
# ----------------------------
season_matrix <- list(
  ONDJFMA = c(10, 11, 12, 1, 2, 3),
  DJF     = c(12, 1, 2),
  MAM     = c(3, 4, 5),
  JJA     = c(6, 7, 8),
  SON     = c(9, 10, 11)
)

stopifnot(season_name %in% names(season_matrix))
season_months <- season_matrix[[season_name]]

# ----------------------------
# 1) OLR extraction + SACZ definition
# ----------------------------
olr_extraction <- function(periodo, varname = "ulwrf") {
  archivo_nc <- olr_file(periodo)
  
  nc <- nc_open(archivo_nc)
  on.exit(nc_close(nc), add = TRUE)
  
  # Required variables
  stopifnot(varname %in% names(nc$var))
  stopifnot(all(c("lon", "lat", "time") %in% names(nc$dim)) || all(c("lon","lat","time") %in% names(nc$var)))
  
  olr  <- ncvar_get(nc, varname)
  lon  <- ncvar_get(nc, "lon")
  lat  <- ncvar_get(nc, "lat")
  time <- ncvar_get(nc, "time")
  
  # Bring time to first dimension (t, lon, lat) assuming original (lon, lat, t)
  # If your file is (lon,lat,time), aperm(c(3,1,2)) is correct.
  olr <- aperm(olr, c(3, 1, 2))
  
  # Flatten spatial dims -> daily mean over the (already subset) SACZ domain
  serie_olr <- matrix(olr, nrow = dim(olr)[1], ncol = prod(dim(olr)[2:3]))
  serie_olr <- data.frame(media = rowMeans(serie_olr, na.rm = TRUE))
  
  # Convert hours since 1800-01-01 to Date
  fechas <- as.POSIXct("1800-01-01 00:00:00", tz = "UTC") + time * 3600
  fechas <- as.Date(fechas)
  
  serie_olr <- serie_olr %>%
    mutate(
      año = as.integer(format(fechas, "%Y")),
      mes = as.integer(format(fechas, "%m")),
      dia = as.integer(format(fechas, "%d")),
      mes_dia = paste(mes, dia, sep = "-")
    )
  
  # Daily climatology by month-day and anomaly
  media_mes_dia <- aggregate(media ~ mes_dia, data = serie_olr, FUN = mean)
  names(media_mes_dia) <- c("mes_dia", "media_mes_dia")
  
  serie_olr <- merge(serie_olr, media_mes_dia, by = "mes_dia")
  serie_olr$indice <- serie_olr$media - serie_olr$media_mes_dia
  
  # SACZ condition (your rule)
  serie_olr$condicion <- "no SACZ"
  serie_olr$condicion[serie_olr$indice < 0] <- "SACZ event"
  
  # Return only needed columns
  serie_olr[, c("año", "mes", "dia", "media", "indice", "condicion")]
}

# Read all periods and bind
serie_olr <- do.call(rbind, lapply(periods, olr_extraction))

# ----------------------------
# 2) Merge with CP labels
# ----------------------------
cp <- read.csv(cp_file)
stopifnot(all(c("año", "mes", "dia", "CP") %in% names(cp)))

x <- merge(serie_olr, cp, by = c("año", "mes", "dia"))
x <- x[, c("dia", "mes", "año", "CP", "condicion", "indice", "media")]

# Save full daily table (useful for reproducibility)
write.csv(x, file.path(fig_dir, "SACZ_CP.csv"), row.names = FALSE)

# Seasonal subset
x_season <- x[x$mes %in% season_months, ]

# ----------------------------
# 3A) p(SACZ | CP): N(SACZ,CP) / N(CP) with Monte Carlo envelope
# ----------------------------
mc_p_SACZ_given_CP <- function(x_all, x_season, season_months, k = 7, n_mc = 1000) {
  # Evaluate both conditions to mimic your original (SACZ vs no SACZ)
  conds <- c("SACZ event", "no SACZ")
  
  out <- list()
  
  for (cond in conds) {
    tmp <- lapply(1:k, function(cp_i) {
      datos_cp <- x_all[x_all$CP == cp_i & x_all$mes %in% season_months, ]
      N <- nrow(datos_cp)
      if (N == 0) {
        return(data.frame(CP = cp_i, p = NA, p97.5 = NA, p2.5 = NA, sign = NA, condicion = cond))
      }
      
      p_obs <- mean(datos_cp$condicion == cond) * 100
      
      mc_vals <- replicate(n_mc, {
        samp <- sample(x_season$condicion, size = N, replace = FALSE)
        mean(samp == cond) * 100
      })
      
      p97.5 <- as.numeric(quantile(mc_vals, 0.975, na.rm = TRUE))
      p2.5  <- as.numeric(quantile(mc_vals, 0.025, na.rm = TRUE))
      
      sign <- "no"
      if (!is.na(p_obs) && p_obs > p97.5) sign <- "superior"
      if (!is.na(p_obs) && p_obs < p2.5)  sign <- "inferior"
      
      data.frame(CP = cp_i,
                 p = round(p_obs, 0),
                 p97.5 = round(p97.5, 0),
                 p2.5  = round(p2.5, 0),
                 sign = sign,
                 condicion = cond)
    })
    
    out[[cond]] <- do.call(rbind, tmp)
  }
  
  do.call(rbind, out)
}

df_A <- mc_p_SACZ_given_CP(x_all = x, x_season = x_season,
                           season_months = season_months, k = k, n_mc = n_mc)

# Optional: hide non-significant cells like your first figure (only for plotting)
df_A_plot <- df_A
df_A_plot$p[df_A_plot$sign == "no"] <- NA
df_A_plot$condicion[df_A_plot$condicion == "SACZ event"] <- "SACZ"

g_A <- ggplot(df_A_plot, aes(x = condicion, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = na.omit(df_A_plot), aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("superior" = "red", "inferior" = "cyan", "no" = NA), name = "") +
  geom_text(aes(label = p), size = 4) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "", y = "CP", title = "p(SACZ | CP) = N(SACZ,CP)/N(CP)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_A)

ggsave(
  filename = file.path(fig_dir, paste0("SACZ_given_CP_", season_name, "_def_Vera.png")),
  plot = g_A, dpi = 300, width = 1.5, height = 3, bg = "white"
)

# ----------------------------
# 3B) p(CP | SACZ): N(SACZ,CP) / N(SACZ) with Monte Carlo envelope
# ----------------------------
mc_p_CP_given_SACZ <- function(x_all, x_season, season_months, k = 7, n_mc = 1000) {
  conds <- c("SACZ event", "no SACZ")
  
  out <- list()
  
  for (cond in conds) {
    datos_cond <- x_all[x_all$condicion == cond & x_all$mes %in% season_months, ]
    N <- nrow(datos_cond)
    
    tmp <- lapply(1:k, function(cp_i) {
      if (N == 0) {
        return(data.frame(CP = cp_i, p = NA, p97.5 = NA, p2.5 = NA, sign = NA, condicion = cond))
      }
      
      p_obs <- mean(datos_cond$CP == cp_i) * 100
      
      mc_vals <- replicate(n_mc, {
        samp <- sample(x_season$CP, size = N, replace = FALSE)
        mean(samp == cp_i) * 100
      })
      
      p97.5 <- as.numeric(quantile(mc_vals, 0.975, na.rm = TRUE))
      p2.5  <- as.numeric(quantile(mc_vals, 0.025, na.rm = TRUE))
      
      sign <- "no"
      if (!is.na(p_obs) && p_obs > p97.5) sign <- "superior"
      if (!is.na(p_obs) && p_obs < p2.5)  sign <- "inferior"
      
      data.frame(CP = cp_i,
                 p = round(p_obs, 0),
                 p97.5 = round(p97.5, 0),
                 p2.5  = round(p2.5, 0),
                 sign = sign,
                 condicion = cond)
    })
    
    out[[cond]] <- do.call(rbind, tmp)
  }
  
  do.call(rbind, out)
}

df_B <- mc_p_CP_given_SACZ(x_all = x, x_season = x_season,
                           season_months = season_months, k = k, n_mc = n_mc)

df_B$condicion[df_B$condicion == "SACZ event"] <- "SACZ"
df_B$season <- season_name

# Save per-season CSV for later rbind panels
write.csv(df_B, file.path(fig_dir, paste0("sacz_", season_name, ".csv")), row.names = FALSE)

g_B <- ggplot(df_B, aes(x = condicion, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_B[!is.na(df_B$p), ], aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("superior" = "red", "inferior" = "cyan", "no" = NA), name = "") +
  geom_text(aes(label = p), size = 4) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "", y = "CP", title = season_name) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_B)

ggsave(
  filename = file.path(fig_dir, paste0("CP_given_SACZ_", season_name, "_def_Vera.png")),
  plot = g_B, dpi = 300, width = 1.3, height = 3, bg = "white"
)

# ----------------------------
# 4) Final multi-season panel (requires existing season CSVs)
# ----------------------------
# NOTE: update the list below to include all seasons you generated (e.g., DJF/MAM/SON/JJA)
df_plot <- rbind(
  read.csv(file.path(fig_dir, "sacz_SON.csv")),
  read.csv(file.path(fig_dir, "sacz_DJF.csv")),
  read.csv(file.path(fig_dir, "sacz_MAM.csv"))
)

df_plot$sign[df_plot$sign == "inferior"] <- "lower than the clim"
df_plot$sign[df_plot$sign == "superior"] <- "higher than the clim"
df_plot$condicion[df_plot$condicion == "no SACZ"] <- "non-SACZ"

g_final <- ggplot(df_plot, aes(x = season, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_plot[df_plot$sign != "no", ], aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  facet_wrap(~ condicion, ncol = 2) +
  scale_color_manual(name = "",
                     values = c("higher than the clim" = "red",
                                "lower than the clim" = "cyan",
                                "no" = NA)) +
  geom_text(data = df_plot[df_plot$sign != "no", ], aes(label = round(p, 1)), size = 4) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "", y = "CP", title = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_final)

ggsave(
  filename = file.path(fig_dir, "SACZ_CP_association_def_Vera_allSeasons.png"),
  plot = g_final, dpi = 300, width = 4, height = 3, bg = "white"
)
