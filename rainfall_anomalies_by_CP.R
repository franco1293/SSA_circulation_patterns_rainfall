rm(list=ls()) ; graphics.off() ; gc()

################################
# Packages
################################
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
library(scales)
library(dplyr)
library(lubridate)
library(abind)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(stringr)
library(zoo)

################################
# Input NetCDF files (CPC precipitation)
################################

# Folder containing CPC NetCDF files
directorio <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st/cpc_regional/"

# List all NetCDF files in the folder
archivos <- list.files(directorio, pattern = "*.nc", full.names = TRUE)

# Read lon/lat from the first file (used to build the grid index)
nc <- nc_open(archivos[1])
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
nc_close(nc)

# Grid table (lon/lat for each gridpoint id = row_number())
df_grid <- expand.grid(lon = lon, lat = lat)

################################
# NetCDF reader: extract daily precip and return a data.frame
################################
leer_datos_nc <- function(archivo) {
  nc <- nc_open(archivo)
  
  # Read time and precipitation variables
  tiempo <- ncvar_get(nc, "time")
  lluvia <- ncvar_get(nc, "precip")
  
  # Convert time (hours since 1900-01-01) to Date
  fechas <- as.Date(tiempo / 24, origin = "1900-01-01")
  
  # Date columns (year/month/day)
  datos <- data.frame(
    año = as.numeric(format(fechas, "%Y")),
    mes = as.numeric(format(fechas, "%m")),
    dia = as.numeric(format(fechas, "%d"))
  )
  
  # Flatten precip: each time slice becomes a vector (gridpoints)
  lluvia <- apply(lluvia, 3, as.vector)
  
  # Combine dates with precip-by-gridpoint columns
  datos_lluvia <- cbind(datos, t(lluvia))
  
  nc_close(nc)
  gc()
  
  return(datos_lluvia)
}

# Read all files and combine into one table (rows = days, columns = gridpoints)
df_lluvia <- bind_rows(lapply(archivos, leer_datos_nc))

# Quick checks
head(df_lluvia[, 1:10])
head(df_grid)

gc()

################################
# Remove daily climatology (for composites)
################################

# Daily date sequence for the full period (must match df_lluvia rows)
fechas <- seq(as.Date("1979-01-01"), as.Date("2023-12-31"), by = "day")

# Add Date and a calendar-day key (mm-dd)
df_z <- df_lluvia
df_z$fecha <- fechas
df_z[is.na(df_z)] <- 0

df_z <- df_z %>%
  mutate(mes_dia = format(fecha, "%m-%d"))

# Daily climatology using only wet days (precip >= 1 mm)
promedios_diarios <- df_z %>%
  group_by(mes_dia) %>%
  summarise(across(matches("^[1-9]"), ~ mean(if_else(. >= 1, ., NA_real_), na.rm = TRUE)))

# Order mes_dia chronologically (for smoothing)
promedios_diarios <- promedios_diarios %>%
  mutate(mes_dia = as.Date(mes_dia, format = "%m-%d")) %>%
  arrange(mes_dia)

# Moving-average smoothing (circular wrap at year boundaries)
n <- nrow(promedios_diarios)
extend <- floor(29 / 2)

promedios_extendidos <- bind_rows(
  tail(promedios_diarios, extend),
  promedios_diarios,
  head(promedios_diarios, extend)
)

promedios_suavizados <- promedios_extendidos %>%
  mutate(across(where(is.numeric),
                ~ rollapply(.x, width = 29, FUN = mean, align = "center", fill = NA, na.rm = TRUE)))

# Keep only the smoothed central part
promedios_moviles <- promedios_suavizados[(extend + 1):(extend + n), ]

# Replace NA with 0 only where the gridpoint exists (avoid changing pure-ocean columns)
promedios_diarios <- promedios_moviles %>%
  mutate(across(where(~ is.numeric(.) && any(!is.na(.))), ~ replace_na(., 0)))

# Quick diagnostic plot for one gridpoint column (example)
plot(y = promedios_diarios$'10', x = 1:length(promedios_diarios$'10'), type = "l")

# Restore mes_dia format as mm-dd for the join
promedios_diarios <- promedios_diarios %>%
  mutate(mes_dia = format(mes_dia, "%m-%d"))

# Build anomalies: set precip < 1 mm to NA before subtraction
df_z <- df_z %>%
  mutate(across(-c(año, mes, dia, fecha, mes_dia), ~ if_else(. < 1, NA_real_, .)))

anomalías <- df_z %>%
  left_join(promedios_diarios, by = "mes_dia", suffix = c("_original", "_promedio"))

promedio <- anomalías %>%
  select(contains("_promedio"))

anomalías <- anomalías %>%
  select(contains(c("año", "mes", "dia", "_original")))

anomalías$mes_dia <- NULL
gc()

# Final anomalies table (precip - daily climatology)
df_lluvia <- anomalías[, 4:ncol(anomalías)] - promedio
df_lluvia$año <- anomalías$año
df_lluvia$mes <- anomalías$mes
df_lluvia$dia <- anomalías$dia

dim(df_lluvia)

rm(anomalías, promedios_diarios, promedio)
gc()



################################
# Merge precipitation with CP classification (from wind clustering)
################################

# Choose classification type (k-means or PCA+k-means) by selecting the correct CSV
k <- 7
archivo_csv <- paste("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures/dias_clasificados_UV_PCA_k-means", k, ".csv", sep="")

# CP labels per day (year/month/day + CP)
df_clasificados <- read.csv(archivo_csv)

# Merge precipitation anomalies with CP labels
df_lluvia <- merge(df_lluvia, df_clasificados, by = c("año", "mes", "dia"))

head(df_lluvia[, c(1,2,3,4,5,ncol(df_lluvia))])
gc()




################################
# Mean precipitation anomaly per CP (all days)
################################

df_lluvia_CP <- df_lluvia[, 4:ncol(df_lluvia)] %>%
  group_by(CP) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

head(df_lluvia_CP[, 1:10])

# Number of valid values per CP (non-NA)
count_lluvia_CP <- df_lluvia[, 4:ncol(df_lluvia)] %>%
  group_by(CP) %>%
  summarise(across(everything(), ~ sum(!is.na(.))))
count_lluvia_CP




###############################  ANNUAL

# Long format: one row per (CP, gridpoint)
promedio_lluvia_por_CP_long <- df_lluvia_CP %>%
  pivot_longer(cols = -CP, names_to = "gridpoint", values_to = "lluvia_promedio") %>%
  mutate(gridpoint = str_remove(gridpoint, "_original") %>% as.numeric())

# Attach lon/lat by gridpoint id
df_mapa <- promedio_lluvia_por_CP_long %>%
  left_join(df_grid %>% mutate(gridpoint = row_number()), by = "gridpoint")

df_mapa[is.na(df_mapa)] <- 0
gc()

# Country boundaries
world <- map_data("world")

# Map (free scale)
g <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = lluvia_promedio), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "darkgoldenrod", high = "green4", mid = "white",
                       midpoint = 0, name = "mm") +
  facet_wrap(~ CP, nrow = 3) +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  coord_fixed(xlim = c(min(df_mapa$lon - 360), max(df_mapa$lon - 360)),
              ylim = c(min(df_mapa$lat), max(df_mapa$lat)), expand = 0) +
  theme(legend.position = "bottom")

plot(g)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")

# Export CP-mean precipitation fields
write.csv(df_mapa, paste("lluvia_por_CP", as.character(k), ".csv", sep=""), row.names = F)

# Export map figure
ggsave(plot = g, filename = paste("lluvia_por_CP", tecnica, k, ".png"), width = 7, height = 5, dpi = 200)

################################
# ANNUAL map with fixed color limits
################################

limit_llu <- c(-4, 4)

g <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = lluvia_promedio), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "darkgoldenrod", high = "green4", mid = "white",
                       midpoint = 0, name = "mm",
                       limits = limit_llu, oob = squish) +
  facet_wrap(~ CP, nrow = k) +
  labs(title = "", x = "", y = "") +
  theme_bw() + ggtitle("ANNUAL") +
  coord_fixed(xlim = c(min(df_mapa$lon - 360), max(df_mapa$lon - 360)),
              ylim = c(min(df_mapa$lat), max(df_mapa$lat)), expand = 0) +
  theme(legend.position = "right",
        legend.margin = margin(c(-0.7, 0, 0, 0), unit = "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(2, "cm"))

plot(g)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
ggsave(plot = g, filename = paste("lluvia_por_CP", tecnica, k, "esc_fija.png"), width = 12, height = 5, dpi = 200)





###############################  DJF

df_lluvia_CP <- df_lluvia[df_lluvia$mes %in% c(12,1,2), 4:ncol(df_lluvia)] %>%
  group_by(CP) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

head(df_lluvia_CP[, 1:10])

promedio_lluvia_por_CP_long <- df_lluvia_CP %>%
  pivot_longer(cols = -CP, names_to = "gridpoint", values_to = "lluvia_promedio") %>%
  mutate(gridpoint = str_remove(gridpoint, "_original") %>% as.numeric())

df_mapa <- promedio_lluvia_por_CP_long %>%
  left_join(df_grid %>% mutate(gridpoint = row_number()), by = "gridpoint")

df_mapa[is.na(df_mapa)] <- 0
gc()

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
write.csv(df_mapa, "lluvia_CP_DJF.csv", row.names = F)

world <- map_data("world")

g1 <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = lluvia_promedio), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "darkgoldenrod", high = "green4", mid = "white",
                       midpoint = 0, name = "mm",
                       limits = limit_llu, oob = squish) +
  facet_wrap(~ CP, nrow = k, strip.position = "left") +
  labs(title = "", x = "", y = "") +
  theme_bw() + ggtitle("DJF") +
  coord_fixed(xlim = c(min(lon - 360), max(lon - 360)),
              ylim = c(min(lat), max(lat)), expand = 0) +
  theme(legend.position = "none",
        strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside")

plot(g1)

###############################  MAM

df_lluvia_CP <- df_lluvia[df_lluvia$mes %in% c(3,4,5), 4:ncol(df_lluvia)] %>%
  group_by(CP) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

head(df_lluvia_CP[, 1:10])

promedio_lluvia_por_CP_long <- df_lluvia_CP %>%
  pivot_longer(cols = -CP, names_to = "gridpoint", values_to = "lluvia_promedio") %>%
  mutate(gridpoint = str_remove(gridpoint, "_original") %>% as.numeric())

df_mapa <- promedio_lluvia_por_CP_long %>%
  left_join(df_grid %>% mutate(gridpoint = row_number()), by = "gridpoint")

df_mapa[is.na(df_mapa)] <- 0
gc()

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
write.csv(df_mapa, "lluvia_CP_MAM.csv", row.names = F)

world <- map_data("world")

g2 <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = lluvia_promedio), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "darkgoldenrod", high = "green4", mid = "white",
                       midpoint = 0, name = "mm",
                       limits = limit_llu, oob = squish) +
  facet_wrap(~ CP, nrow = k) +
  labs(title = "", x = "", y = "") +
  theme_bw() + ggtitle("MAM") +
  coord_fixed(xlim = c(min(lon - 360), max(lon - 360)),
              ylim = c(min(lat), max(lat)), expand = 0) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text.y = element_blank())

plot(g2)

###############################  JJA

df_lluvia_CP <- df_lluvia[df_lluvia$mes %in% c(6,7,8), 4:ncol(df_lluvia)] %>%
  group_by(CP) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

head(df_lluvia_CP[, 1:10])

count_lluvia_CP <- df_lluvia[df_lluvia$mes %in% c(6,7,8), 4:ncol(df_lluvia)] %>%
  group_by(CP) %>%
  summarise(across(everything(), ~ sum(!is.na(.))))
count_lluvia_CP

# Ensure all CPs exist (fill missing CPs with NA rows)
todos_los_CP <- data.frame(CP = 1:k)
df_lluvia_CP <- todos_los_CP %>%
  left_join(df_lluvia_CP, by = "CP")

promedio_lluvia_por_CP_long <- df_lluvia_CP %>%
  pivot_longer(cols = -CP, names_to = "gridpoint", values_to = "lluvia_promedio") %>%
  mutate(gridpoint = str_remove(gridpoint, "_original") %>% as.numeric())

df_mapa <- promedio_lluvia_por_CP_long %>%
  left_join(df_grid %>% mutate(gridpoint = row_number()), by = "gridpoint")

df_mapa[is.na(df_mapa)] <- 0
gc()

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
write.csv(df_mapa, "lluvia_CP_JJA.csv", row.names = F)

world <- map_data("world")

g3 <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = lluvia_promedio), alpha = 1, na.rm = FALSE) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "darkgoldenrod", high = "green4", mid = "white",
                       midpoint = 0, name = "mm",
                       limits = limit_llu, oob = squish) +
  facet_wrap(~ CP, nrow = k) +
  labs(title = "", x = "", y = "") +
  theme_bw() + ggtitle("JJA") +
  coord_fixed(xlim = c(min(lon - 360), max(lon - 360)),
              ylim = c(min(lat), max(lat)), expand = 0) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.text.y = element_blank())

plot(g3)

###############################  SON

df_lluvia_CP <- df_lluvia[df_lluvia$mes %in% c(9,10,11), 4:ncol(df_lluvia)] %>%
  group_by(CP) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

head(df_lluvia_CP[, 1:10])

promedio_lluvia_por_CP_long <- df_lluvia_CP %>%
  pivot_longer(cols = -CP, names_to = "gridpoint", values_to = "lluvia_promedio") %>%
  mutate(gridpoint = str_remove(gridpoint, "_original") %>% as.numeric())

df_mapa <- promedio_lluvia_por_CP_long %>%
  left_join(df_grid %>% mutate(gridpoint = row_number()), by = "gridpoint")

df_mapa[is.na(df_mapa)] <- 0
gc()

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
write.csv(df_mapa, "lluvia_CP_SON.csv", row.names = F)

world <- map_data("world")

g4 <- ggplot(data = df_mapa) +
  geom_tile(aes(x = lon - 360, y = lat, fill = lluvia_promedio), alpha = 1) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .3) +
  scale_fill_gradient2(low = "darkgoldenrod", high = "green4", mid = "white",
                       midpoint = 0, name = "mm",
                       limits = limit_llu, oob = squish) +
  facet_wrap(~ CP, nrow = k) +
  labs(title = "", x = "", y = "") +
  theme_bw() + ggtitle("SON") +
  coord_fixed(xlim = c(min(lon - 360), max(lon - 360)),
              ylim = c(min(lat), max(lat)), expand = 0) +
  theme(legend.position = "right",
        strip.text = element_blank(),
        axis.text.y = element_blank())

plot(g4)

################################
# Combine seasonal maps and export
################################

graficos <- list(g1, g2, g3, g4) %>%
  lapply(function(grafico) {
    grafico + theme(plot.margin = unit(c(-1, 0, -1, -0.5), "cm"))
  })

g_f <- grid.arrange(
  grobs = graficos,
  nrow = 1,
  widths = c(1.345, 1, 1, 1.6),
  padding = unit(0, "line")
)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
ggsave(plot = g_f, filename = paste("lluvia_por_CP_seasons", tecnica, k, ".png"), height = 7, width = 6, dpi = 200)
