rm(list=ls()) ; graphics.off() ; gc()

##############################################
# CPC precipitation extremes (p95, moving window)
# + Monte Carlo significance test by CP (seasonal)
##############################################

##############################################
# Packages
##############################################
library(ncdf4)
library(ggplot2)
library(gridExtra)
library(fields)
library(mapdata)
library(matrixStats)
library(metR)
library(RColorBrewer)
library(miceadds)
library(tidyr)
library(dplyr)
library(lubridate)
library(abind)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(stringr)
library(scales)
library(zoo)

##############################################
# 1) Read CPC NetCDF files and build daily table
##############################################

# Folder with CPC NetCDF files
directorio <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st/cpc_regional/"

# List all NetCDF files
archivos <- list.files(directorio, pattern = "*.nc", full.names = TRUE)

# Read lon/lat from the first file and create a grid index
nc <- nc_open(archivos[1])
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

df_grid <- expand.grid(lon = lon, lat = lat)
df_grid$grid <- 1:nrow(df_grid)   # gridpoint id (used later)

# NetCDF reader: returns one row per day, columns = gridpoints
leer_datos_nc <- function(archivo) {
  nc <- nc_open(archivo)
  
  # Read time and precipitation
  tiempo <- ncvar_get(nc, "time")
  lluvia <- ncvar_get(nc, "precip")
  
  # Convert time (hours since 1900-01-01) to Date
  fechas <- as.Date(tiempo / 24, origin = "1900-01-01")
  
  # Year/month/day
  datos <- data.frame(
    año = as.numeric(format(fechas, "%Y")),
    mes = as.numeric(format(fechas, "%m")),
    dia = as.numeric(format(fechas, "%d"))
  )
  
  # Flatten precip: each time slice becomes a vector of gridpoints
  lluvia <- apply(lluvia, 3, as.vector)
  
  # Combine dates + precip columns (one column per gridpoint)
  datos_lluvia <- cbind(datos, t(lluvia))
  
  nc_close(nc)
  gc()
  
  return(datos_lluvia)
}

# Read all files and stack into one dataframe (daily)
df_lluvia <- bind_rows(lapply(archivos, leer_datos_nc))

# Quick checks
head(df_lluvia[,1:10])
head(df_grid)

gc()

##############################################
# 2) Compute p95 thresholds using a 29-day moving window (wet days only)
##############################################

# Convert to data.table for speed
DT <- as.data.table(df_lluvia)

# Build Date and day-of-year (DOY)
DT[, date := as.Date(paste(año, mes, dia, sep = "-"))]
DT[, doy := yday(date)]

# Precip columns = everything except date metadata
precip_cols <- setdiff(names(DT), c("año", "mes", "dia", "date", "doy"))

# Convert to long format: one row per (day, gridpoint)
DT_long <- melt(
  DT,
  id.vars = c("año", "mes", "dia", "date", "doy"),
  measure.vars = precip_cols,
  variable.name = "grid",
  value.name = "precip"
)

# Keep only wet days (precip > 1 mm) and remove NA
DT_long <- DT_long[!is.na(precip) & precip > 1]

# Helper: DOY window indices (circular wrapping)
get_window <- function(d, window = 29, total = 365) {
  half <- floor(window / 2)
  indices <- ((d - half - 1):(d + half - 1)) %% total + 1
  return(indices)
}

# For each gridpoint and DOY: compute p95 within the 29-day moving window
thresholds <- DT_long[, {
  thresh <- sapply(1:365, function(d) {
    window_days <- get_window(d, window = 29, total = 365)
    vals <- precip[doy %in% window_days]
    if (length(vals) == 0) NA_real_ else as.numeric(quantile(vals, 0.95, na.rm = TRUE))
  })
  .(doy = 1:365, threshold = thresh)
}, by = grid]

# Attach threshold to each record (grid + doy)
DT_long <- merge(DT_long, thresholds, by = c("grid", "doy"), all.x = TRUE)

# Flag exceedances (extreme days)
DT_long[, exceed := precip > threshold]
DT_long$grid <- as.numeric(DT_long$grid)

# Keep only extreme events (optional output table)
DT_extremos <- DT_long[exceed == TRUE]

# Quick check / diagnostic
head(DT_extremos)

# Sort (useful for debugging)
setorder(DT_long, date, grid)

# Example: inspect one gridpoint (threshold over time)
example <- DT_long[DT_long$grid == 30,]
ggplot(data = example[1:300,]) + geom_line(aes(x = 1:300, y = threshold))

# Export extremes and grid table
setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
write.csv(DT_extremos, "dias_p95_CPC.csv", row.names = FALSE)
write.csv(df_grid, "gridpoints_p95_CPC.csv", row.names = FALSE)

rm(list=ls()) ; graphics.off() ; gc()








#####################################################################
# 3) Monte Carlo test (Olmo & Bettolli 2021-style)
#    Goal: for each season and CP, test whether observed extreme-day
#    frequency is higher than expected by random sampling of days.
#####################################################################

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(parallel)
library(ggplot2)
library(mapdata)
library(scales)
library(gridExtra)

##############################################
# Inputs: CP daily classification + CPC extremes
##############################################

k <- 7
tecnica <- "PCA_K_means_UV" # PCA_solo, PCA_k-means, k-means, PCA_K_means_UV, UV_K_means

# Daily CP labels (from wind clustering script)
archivo_csv <- paste(
  "C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures/dias_clasificados_UV_PCA_k-means",
  k, ".csv", sep = ""
)
df_clasificados <- read.csv(archivo_csv)

# CPC extremes (p95 exceedances) and grid coordinates
df_lluvia <- read.csv("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures/dias_p95_CPC.csv")
df_grid  <- read.csv("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures/gridpoints_p95_CPC.csv")

# Merge extremes with CP labels by date
df_lluvia <- merge(df_lluvia, df_clasificados, by = c("año", "mes", "dia"))
gc()

##############################################
# Utility: one season workflow (copied for DJF/MAM/JJA/SON)
# - Build complete day x grid table with exceed=0/1
# - Compute observed % extremes by CP and grid
# - Monte Carlo CI by random day sampling (same N per CP)
# - Keep only significant "high" cells (observed > 97.5%)
##############################################

###############################  DJF

# 0) Seasonal subsets
df_cp_djf  <- df_clasificados[df_clasificados$mes %in% c(12,1,2), c("año","mes","dia","CP")]
df_evt_djf <- df_lluvia[df_lluvia$mes %in% c(12,1,2), c("año","mes","dia","grid","exceed")]

# Ensure types
df_cp_djf$año <- as.integer(df_cp_djf$año); df_cp_djf$mes <- as.integer(df_cp_djf$mes); df_cp_djf$dia <- as.integer(df_cp_djf$dia)
df_cp_djf$CP  <- as.integer(df_cp_djf$CP)

df_evt_djf$año <- as.integer(df_evt_djf$año); df_evt_djf$mes <- as.integer(df_evt_djf$mes); df_evt_djf$dia <- as.integer(df_evt_djf$dia)
df_evt_djf$grid <- as.integer(df_evt_djf$grid)
df_evt_djf$exceed <- as.integer(df_evt_djf$exceed)  # TRUE/FALSE -> 1/0

df_grid$grid <- as.integer(df_grid$grid)

# 1) N days per CP (from classification table)
cp_count <- df_cp_djf %>% group_by(CP) %>% summarise(N = n(), .groups="drop")

# 2) Build complete table (all DJF days x all gridpoints), exceed = 0/1
DT_days <- as.data.table(unique(df_cp_djf))
DT_days[, day_id := .I]

DT_grid <- as.data.table(df_grid)[, .(grid, lon, lat)]

DTs <- CJ(day_id = DT_days$day_id, grid = DT_grid$grid)
DTs <- merge(DTs, DT_days, by = "day_id", all.x = TRUE)

DT_evt <- as.data.table(df_evt_djf)
DT_evt <- DT_evt[, .(exceed = 1L), by = .(año, mes, dia, grid)]  # collapse duplicates if any

setkey(DTs, año, mes, dia, grid)
setkey(DT_evt, año, mes, dia, grid)

DTs <- DT_evt[DTs]                 # left join: bring exceed where present
DTs[is.na(exceed), exceed := 0L]   # missing = non-extreme

# 3) Observed % extremes per CP and grid
df_lluvia_CP <- DTs[, .(exceed = sum(exceed)), by = .(CP, grid)]
df_mapa <- merge(df_lluvia_CP, DT_grid, by="grid", all.x=TRUE)
df_mapa <- merge(df_mapa, as.data.table(cp_count), by="CP", all.x=TRUE)
df_mapa[, exceed := exceed / N * 100]

gc()

# 4) Monte Carlo: random sampling of days (same N per CP) to build CI
n_iter <- 5000
set.seed(1234)

n_days_total <- nrow(DT_days)
cp_levels <- sort(unique(DT_days$CP))

N_byCP <- as.data.table(cp_count); setkey(N_byCP, CP)

# Precompute samples: for each CP, sample Ncp day_ids repeatedly
samples_byCP <- lapply(cp_levels, function(cp_i){
  Ncp <- N_byCP[CP == cp_i, N]
  replicate(n_iter, sample.int(n_days_total, size = Ncp, replace = FALSE), simplify = FALSE)
})
names(samples_byCP) <- as.character(cp_levels)

# Parallel loop over gridpoints
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

res_ci <- foreach(g = unique(DTs$grid), .combine = rbind,
                  .packages = c("data.table")) %dopar% {
                    
                    # Exceed vector (0/1) for all days at this grid
                    tmp <- DTs[grid == g, .(day_id, exceed)]
                    tmp <- tmp[, .(exceed = max(exceed)), by = day_id]
                    
                    ex_vec <- integer(n_days_total)
                    ex_vec[tmp$day_id] <- tmp$exceed
                    
                    out_list <- vector("list", length(cp_levels))
                    
                    for (j in seq_along(cp_levels)) {
                      cp_i <- cp_levels[j]
                      samp_list <- samples_byCP[[as.character(cp_i)]]
                      vals <- vapply(samp_list, function(idx) mean(ex_vec[idx]) * 100, numeric(1))
                      
                      out_list[[j]] <- data.table(
                        grid = g,
                        CP = cp_i,
                        p2.5  = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
                        p97.5 = as.numeric(quantile(vals, 0.975, na.rm = TRUE))
                      )
                    }
                    
                    rbindlist(out_list)
                  }

stopCluster(cl)

# Compare observed vs CI
obs <- DTs[, .(exceed_obs = mean(exceed) * 100), by = .(CP, grid)]

sig_table <- merge(obs, res_ci, by = c("grid","CP"), all.x = TRUE)
sig_table[, sig := fifelse(exceed_obs > p97.5, "high",
                           fifelse(exceed_obs < p2.5,  "low", "ns"))]

# Mask map: keep only significant "high" cells
df_mapa <- merge(as.data.table(df_mapa),
                 sig_table[, .(grid, CP, sig)],
                 by = c("grid","CP"), all.x = TRUE)

df_mapa[, exceed := fifelse(sig == "high", exceed, 0)]
df_mapa[, sig := NULL]

# Plot DJF
world <- map_data("world")
limits_sel <- c(0,10)

g1 <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = exceed), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "red4", high = "blue4", mid="white",
                       na.value = "white", midpoint = 0, name = "",
                       limits = limits_sel, oob = squish) +
  facet_wrap(~ CP, nrow = 1) +
  labs(title = "DJF", x = "", y = "") +
  theme_bw() +
  coord_fixed(xlim = c(min(df_mapa$lon - 360), max(df_mapa$lon - 360)),
              ylim = c(min(df_mapa$lat), max(df_mapa$lat)), expand = 0) +
  theme(legend.position = "none")

plot(g1)

###############################  MAM

df_cp_mam  <- df_clasificados[df_clasificados$mes %in% c(3,4,5), c("año","mes","dia","CP")]
df_evt_mam <- df_lluvia[df_lluvia$mes %in% c(3,4,5), c("año","mes","dia","grid","exceed")]

df_cp_mam$año <- as.integer(df_cp_mam$año); df_cp_mam$mes <- as.integer(df_cp_mam$mes); df_cp_mam$dia <- as.integer(df_cp_mam$dia)
df_cp_mam$CP  <- as.integer(df_cp_mam$CP)

df_evt_mam$año <- as.integer(df_evt_mam$año); df_evt_mam$mes <- as.integer(df_evt_mam$mes); df_evt_mam$dia <- as.integer(df_evt_mam$dia)
df_evt_mam$grid <- as.integer(df_evt_mam$grid)
df_evt_mam$exceed <- as.integer(df_evt_mam$exceed)

df_grid$grid <- as.integer(df_grid$grid)

cp_count <- df_cp_mam %>% group_by(CP) %>% summarise(N = n(), .groups="drop")

DT_days <- as.data.table(unique(df_cp_mam))
DT_days[, day_id := .I]

DT_grid <- as.data.table(df_grid)[, .(grid, lon, lat)]

DTs <- CJ(day_id = DT_days$day_id, grid = DT_grid$grid)
DTs <- merge(DTs, DT_days, by = "day_id", all.x = TRUE)

DT_evt <- as.data.table(df_evt_mam)
DT_evt <- DT_evt[, .(exceed = 1L), by = .(año, mes, dia, grid)]

setkey(DTs, año, mes, dia, grid)
setkey(DT_evt, año, mes, dia, grid)

DTs <- DT_evt[DTs]
DTs[is.na(exceed), exceed := 0L]

df_lluvia_CP <- DTs[, .(exceed = sum(exceed)), by = .(CP, grid)]
df_mapa <- merge(df_lluvia_CP, DT_grid, by="grid", all.x=TRUE)
df_mapa <- merge(df_mapa, as.data.table(cp_count), by="CP", all.x=TRUE)
df_mapa[, exceed := exceed / N * 100]

gc()

n_iter <- 5000
set.seed(1234)

n_days_total <- nrow(DT_days)
cp_levels <- sort(unique(DT_days$CP))

N_byCP <- as.data.table(cp_count); setkey(N_byCP, CP)

samples_byCP <- lapply(cp_levels, function(cp_i){
  Ncp <- N_byCP[CP == cp_i, N]
  replicate(n_iter, sample.int(n_days_total, size = Ncp, replace = FALSE), simplify = FALSE)
})
names(samples_byCP) <- as.character(cp_levels)

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

res_ci <- foreach(g = unique(DTs$grid), .combine = rbind,
                  .packages = c("data.table")) %dopar% {
                    
                    tmp <- DTs[grid == g, .(day_id, exceed)]
                    tmp <- tmp[, .(exceed = max(exceed)), by = day_id]
                    
                    ex_vec <- integer(n_days_total)
                    ex_vec[tmp$day_id] <- tmp$exceed
                    
                    out_list <- vector("list", length(cp_levels))
                    
                    for (j in seq_along(cp_levels)) {
                      cp_i <- cp_levels[j]
                      samp_list <- samples_byCP[[as.character(cp_i)]]
                      vals <- vapply(samp_list, function(idx) mean(ex_vec[idx]) * 100, numeric(1))
                      
                      out_list[[j]] <- data.table(
                        grid = g,
                        CP = cp_i,
                        p2.5  = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
                        p97.5 = as.numeric(quantile(vals, 0.975, na.rm = TRUE))
                      )
                    }
                    
                    rbindlist(out_list)
                  }

stopCluster(cl)

obs <- DTs[, .(exceed_obs = mean(exceed) * 100), by = .(CP, grid)]

sig_table <- merge(obs, res_ci, by = c("grid","CP"), all.x = TRUE)
sig_table[, sig := fifelse(exceed_obs > p97.5, "high",
                           fifelse(exceed_obs < p2.5,  "low", "ns"))]

df_mapa <- merge(as.data.table(df_mapa),
                 sig_table[, .(grid, CP, sig)],
                 by = c("grid","CP"), all.x = TRUE)

df_mapa[, exceed := fifelse(sig == "high", exceed, 0)]
df_mapa[, sig := NULL]

world <- map_data("world")
limits_sel <- c(0,10)

g2 <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = exceed), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "red4", high = "blue4", mid="white",
                       na.value = "white", midpoint = 0, name = "",
                       limits = limits_sel, oob = squish) +
  facet_wrap(~ CP, nrow = 1) +
  labs(title = "MAM", x = "", y = "") +
  theme_bw() +
  coord_fixed(xlim = c(min(df_mapa$lon - 360), max(df_mapa$lon - 360)),
              ylim = c(min(df_mapa$lat), max(df_mapa$lat)), expand = 0) +
  theme(legend.position = "none")

plot(g2)

###############################  JJA

df_cp_jja  <- df_clasificados[df_clasificados$mes %in% c(6,7,8), c("año","mes","dia","CP")]
df_evt_jja <- df_lluvia[df_lluvia$mes %in% c(6,7,8), c("año","mes","dia","grid","exceed")]

df_cp_jja$año <- as.integer(df_cp_jja$año); df_cp_jja$mes <- as.integer(df_cp_jja$mes); df_cp_jja$dia <- as.integer(df_cp_jja$dia)
df_cp_jja$CP  <- as.integer(df_cp_jja$CP)

df_evt_jja$año <- as.integer(df_evt_jja$año); df_evt_jja$mes <- as.integer(df_evt_jja$mes); df_evt_jja$dia <- as.integer(df_evt_jja$dia)
df_evt_jja$grid <- as.integer(df_evt_jja$grid)
df_evt_jja$exceed <- as.integer(df_evt_jja$exceed)

df_grid$grid <- as.integer(df_grid$grid)

cp_count <- df_cp_jja %>% group_by(CP) %>% summarise(N = n(), .groups="drop")

DT_days <- as.data.table(unique(df_cp_jja))
DT_days[, day_id := .I]

DT_grid <- as.data.table(df_grid)[, .(grid, lon, lat)]

DTs <- CJ(day_id = DT_days$day_id, grid = DT_grid$grid)
DTs <- merge(DTs, DT_days, by = "day_id", all.x = TRUE)

DT_evt <- as.data.table(df_evt_jja)
DT_evt <- DT_evt[, .(exceed = 1L), by = .(año, mes, dia, grid)]

setkey(DTs, año, mes, dia, grid)
setkey(DT_evt, año, mes, dia, grid)

DTs <- DT_evt[DTs]
DTs[is.na(exceed), exceed := 0L]

df_lluvia_CP <- DTs[, .(exceed = sum(exceed)), by = .(CP, grid)]
df_mapa <- merge(df_lluvia_CP, DT_grid, by="grid", all.x=TRUE)
df_mapa <- merge(df_mapa, as.data.table(cp_count), by="CP", all.x=TRUE)
df_mapa[, exceed := exceed / N * 100]

gc()

n_iter <- 5000
set.seed(1234)

n_days_total <- nrow(DT_days)
cp_levels <- sort(unique(DT_days$CP))

N_byCP <- as.data.table(cp_count); setkey(N_byCP, CP)

samples_byCP <- lapply(cp_levels, function(cp_i){
  Ncp <- N_byCP[CP == cp_i, N]
  replicate(n_iter, sample.int(n_days_total, size = Ncp, replace = FALSE), simplify = FALSE)
})
names(samples_byCP) <- as.character(cp_levels)

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

res_ci <- foreach(g = unique(DTs$grid), .combine = rbind,
                  .packages = c("data.table")) %dopar% {
                    
                    tmp <- DTs[grid == g, .(day_id, exceed)]
                    tmp <- tmp[, .(exceed = max(exceed)), by = day_id]
                    
                    ex_vec <- integer(n_days_total)
                    ex_vec[tmp$day_id] <- tmp$exceed
                    
                    out_list <- vector("list", length(cp_levels))
                    
                    for (j in seq_along(cp_levels)) {
                      cp_i <- cp_levels[j]
                      samp_list <- samples_byCP[[as.character(cp_i)]]
                      vals <- vapply(samp_list, function(idx) mean(ex_vec[idx]) * 100, numeric(1))
                      
                      out_list[[j]] <- data.table(
                        grid = g,
                        CP = cp_i,
                        p2.5  = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
                        p97.5 = as.numeric(quantile(vals, 0.975, na.rm = TRUE))
                      )
                    }
                    
                    rbindlist(out_list)
                  }

stopCluster(cl)

obs <- DTs[, .(exceed_obs = mean(exceed) * 100), by = .(CP, grid)]

sig_table <- merge(obs, res_ci, by = c("grid","CP"), all.x = TRUE)
sig_table[, sig := fifelse(exceed_obs > p97.5, "high",
                           fifelse(exceed_obs < p2.5,  "low", "ns"))]

df_mapa <- merge(as.data.table(df_mapa),
                 sig_table[, .(grid, CP, sig)],
                 by = c("grid","CP"), all.x = TRUE)

df_mapa[, exceed := fifelse(sig == "high", exceed, 0)]
df_mapa[, sig := NULL]

world <- map_data("world")
limits_sel <- c(0,10)

g3 <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = exceed), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "red4", high = "blue4", mid="white",
                       na.value = "white", midpoint = 0, name = "",
                       limits = limits_sel, oob = squish) +
  facet_wrap(~ CP, nrow = 1) +
  labs(title = "JJA", x = "", y = "") +
  theme_bw() +
  coord_fixed(xlim = c(min(df_mapa$lon - 360), max(df_mapa$lon - 360)),
              ylim = c(min(df_mapa$lat), max(df_mapa$lat)), expand = 0) +
  theme(legend.position = "none")

plot(g3)

###############################  SON

df_cp_son  <- df_clasificados[df_clasificados$mes %in% c(9,10,11), c("año","mes","dia","CP")]
df_evt_son <- df_lluvia[df_lluvia$mes %in% c(9,10,11), c("año","mes","dia","grid","exceed")]

df_cp_son$año <- as.integer(df_cp_son$año); df_cp_son$mes <- as.integer(df_cp_son$mes); df_cp_son$dia <- as.integer(df_cp_son$dia)
df_cp_son$CP  <- as.integer(df_cp_son$CP)

df_evt_son$año <- as.integer(df_evt_son$año); df_evt_son$mes <- as.integer(df_evt_son$mes); df_evt_son$dia <- as.integer(df_evt_son$dia)
df_evt_son$grid <- as.integer(df_evt_son$grid)
df_evt_son$exceed <- as.integer(df_evt_son$exceed)

df_grid$grid <- as.integer(df_grid$grid)

cp_count <- df_cp_son %>% group_by(CP) %>% summarise(N = n(), .groups="drop")

DT_days <- as.data.table(unique(df_cp_son))
DT_days[, day_id := .I]

DT_grid <- as.data.table(df_grid)[, .(grid, lon, lat)]

DTs <- CJ(day_id = DT_days$day_id, grid = DT_grid$grid)
DTs <- merge(DTs, DT_days, by = "day_id", all.x = TRUE)

DT_evt <- as.data.table(df_evt_son)
DT_evt <- DT_evt[, .(exceed = 1L), by = .(año, mes, dia, grid)]

setkey(DTs, año, mes, dia, grid)
setkey(DT_evt, año, mes, dia, grid)

DTs <- DT_evt[DTs]
DTs[is.na(exceed), exceed := 0L]

df_lluvia_CP <- DTs[, .(exceed = sum(exceed)), by = .(CP, grid)]
df_mapa <- merge(df_lluvia_CP, DT_grid, by="grid", all.x=TRUE)
df_mapa <- merge(df_mapa, as.data.table(cp_count), by="CP", all.x=TRUE)
df_mapa[, exceed := exceed / N * 100]

gc()

n_iter <- 5000
set.seed(1234)

n_days_total <- nrow(DT_days)
cp_levels <- sort(unique(DT_days$CP))

N_byCP <- as.data.table(cp_count); setkey(N_byCP, CP)

samples_byCP <- lapply(cp_levels, function(cp_i){
  Ncp <- N_byCP[CP == cp_i, N]
  replicate(n_iter, sample.int(n_days_total, size = Ncp, replace = FALSE), simplify = FALSE)
})
names(samples_byCP) <- as.character(cp_levels)

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

res_ci <- foreach(g = unique(DTs$grid), .combine = rbind,
                  .packages = c("data.table")) %dopar% {
                    
                    tmp <- DTs[grid == g, .(day_id, exceed)]
                    tmp <- tmp[, .(exceed = max(exceed)), by = day_id]
                    
                    ex_vec <- integer(n_days_total)
                    ex_vec[tmp$day_id] <- tmp$exceed
                    
                    out_list <- vector("list", length(cp_levels))
                    
                    for (j in seq_along(cp_levels)) {
                      cp_i <- cp_levels[j]
                      samp_list <- samples_byCP[[as.character(cp_i)]]
                      vals <- vapply(samp_list, function(idx) mean(ex_vec[idx]) * 100, numeric(1))
                      
                      out_list[[j]] <- data.table(
                        grid = g,
                        CP = cp_i,
                        p2.5  = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
                        p97.5 = as.numeric(quantile(vals, 0.975, na.rm = TRUE))
                      )
                    }
                    
                    rbindlist(out_list)
                  }

stopCluster(cl)

obs <- DTs[, .(exceed_obs = mean(exceed) * 100), by = .(CP, grid)]

sig_table <- merge(obs, res_ci, by = c("grid","CP"), all.x = TRUE)
sig_table[, sig := fifelse(exceed_obs > p97.5, "high",
                           fifelse(exceed_obs < p2.5,  "low", "ns"))]

df_mapa <- merge(as.data.table(df_mapa),
                 sig_table[, .(grid, CP, sig)],
                 by = c("grid","CP"), all.x = TRUE)

df_mapa[, exceed := fifelse(sig == "high", exceed, 0)]
df_mapa[, sig := NULL]

world <- map_data("world")
limits_sel <- c(0,10)

g4 <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = exceed), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "red4", high = "blue4", mid="white",
                       na.value = "white", midpoint = 0, name = "",
                       limits = limits_sel, oob = squish) +
  facet_wrap(~ CP, nrow = 1) +
  labs(title = "SON", x = "", y = "") +
  theme_bw() +
  coord_fixed(xlim = c(min(df_mapa$lon - 360), max(df_mapa$lon - 360)),
              ylim = c(min(df_mapa$lat), max(df_mapa$lat)), expand = 0) +
  theme(legend.position = "none")

plot(g4)

##############################################
# 4) Combine seasonal panels and export figure
##############################################

graficos <- list(g1, g2, g3, g4)

g_f <- grid.arrange(
  grobs = graficos,
  nrow = 4,
  heights = c(1,1,1,1),
  padding = unit(0, "line")
)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
ggsave(plot = g_f, filename = paste("p95_seasons_", tecnica, k, ".png"), width = 8, height = 8.5, dpi = 200)
