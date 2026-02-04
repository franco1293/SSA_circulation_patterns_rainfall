####################
### K-means + PCA ###
####################

# Daily circulation-type classification using ERA5 u850 and v850 (1979–2023):
# - Read and concatenate NetCDF files
# - Subset to the study region
# - Reshape into [gridpoint × day] matrices
# - Build daily anomalies (remove daily climatology)
# - PCA on standardized fields
# - k-means clustering in retained PC space
# - Cluster validation (pseudo-F, silhouette) and output products (frequencies, composites)

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
#library(loadeR)
#library(visualizeR)
library(miceadds)
library(dplyr)
library(lubridate)
library(abind)
library(parallel)
library(foreach)
library(doParallel)
library(cluster)
library(scales)
library(tidyr)

# NOTE: Before running k-means, decide whether to cluster using:
# - raw standardized fields (z) or
# - filtered anomalies (z_filtrado)
########################################################################################################## u

################################
# Study region limits
################################
lat_min <- -40
lat_max <- -10
lon_min <- -83  # Adjust if your data uses 0–360 longitudes
lon_max <- -40

################################
# Working directory (inputs)
################################
setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/ERA5/")

################################
# Read u850 (1979–1999)
################################
file_path <- "daily u850_1x1_1979-1999.nc"
nc_data <- nc_open(file_path)

# Read coordinates and time axis
lon  <- ncvar_get(nc_data, "longitude")
lat  <- ncvar_get(nc_data, "latitude")
time <- ncvar_get(nc_data, "time")

# Read variable: u component at 850 hPa
z <- ncvar_get(nc_data, "u")

nc_close(nc_data)

################################
# Read u850 (2000–2023)
################################
file_path <- "daily u850_1x1_2000-2023.nc"
nc_data <- nc_open(file_path)

z1 <- ncvar_get(nc_data, "u")

nc_close(nc_data)

# Concatenate both periods along the time dimension
z <- abind::abind(z, z1, along = 3)

rm(z1)
gc()

# Keep data ordered as z[lon, lat, time]
z <- aperm(z, c(1, 2, 3))

# Check final dimensions
print(dim(z))  # expected: [nlon, nlat, ndays]

################################
# Subset region
################################
lat_idx <- which(lat >= lat_min & lat <= lat_max)
lon_idx <- which(lon >= lon_min & lon <= lon_max)

lat <- lat[lat_idx]
lon <- lon[lon_idx]

z <- z[lon_idx, lat_idx, ]
dim(z)

################################
# Reshape into matrix [gridpoint × day]
################################
num_days <- dim(z)[3]
num_gridpoints <- length(lon) * length(lat)

z <- matrix(aperm(z, c(1, 2, 3)),
            nrow = num_gridpoints,
            ncol = num_days)

print(dim(z))  # expected: [num_gridpoints, num_days]

################################
# Daily climatology removal (for composites)
################################
# Build daily date sequence (must match the concatenated data period)
fechas <- seq(as.Date("1979-01-01"), as.Date("2023-12-31"), by = "day")

# Convert matrix to data frame: rows = days, columns = gridpoints
df_z <- as.data.frame(t(z))
df_z$fecha <- fechas

# Calendar-day key used to compute daily means
df_z <- df_z %>%
  mutate(mes_dia = format(fecha, "%m-%d"))

# Daily climatology (mean per calendar day)
promedios_diarios <- df_z %>%
  group_by(mes_dia) %>%
  summarise(across(starts_with("V"), mean))

# Join climatology and build anomalies (original - daily mean)
anomalías <- df_z %>%
  left_join(promedios_diarios, by = "mes_dia", suffix = c("_original", "_promedio"))

promedio <- anomalías %>%
  select(contains("_promedio"))

anomalías <- anomalías %>%
  select(contains("_original"))

gc()

# Filtered field (daily anomalies)
z_filtrado <- anomalías - promedio

################################
# Quick diagnostic plot (single gridpoint time series)
################################
n <- 100 # gridpoint index to plot (example)
g1 <- ggplot() +
  geom_line(aes(x = 1:length(promedio[, n]),  y = promedio[, n] / 10)) +
  geom_line(aes(x = 1:length(anomalías[, n]), y = anomalías[, n] / 10), col = "red", alpha = .5) +
  xlab("days") + ylab("z500")

g2 <- ggplot() +
  geom_line(aes(x = 1:length(z_filtrado[, n]), y = z_filtrado[, n] / 10), col = "blue") +
  xlab("days") + ylab("z500 filtered")

g <- grid.arrange(g1, g2)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")

ggsave(
  filename = paste("pruebas_series_", as.character(n), "_z500.png"),
  plot = g,
  device = "png",
  dpi = 300,
  width = 8,
  height = 5
)

rm(df_z, promedios_diarios, anomalías, nc_data)
gc()

# Convert filtered anomalies back to [gridpoint × day] matrix
zu_filtrado <- t(as.matrix(z_filtrado))
row.names(zu_filtrado) <- NULL
colnames(zu_filtrado) <- NULL

dim(zu_filtrado)

# Store u850 matrix (unfiltered)
zu <- z

rm(z)
gc()

########################################################################################################## v

################################
# Working directory (inputs)
################################
setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/ERA5/")

################################
# Read v850 (1979–1999)
################################
file_path <- "daily v850_1x1_1979-1999.nc"
nc_data <- nc_open(file_path)

lon  <- ncvar_get(nc_data, "longitude")
lat  <- ncvar_get(nc_data, "latitude")
time <- ncvar_get(nc_data, "time")

# Read variable: v component at 850 hPa
z <- ncvar_get(nc_data, "v")

nc_close(nc_data)

################################
# Read v850 (2000–2023)
################################
file_path <- "daily v850_1x1_2000-2023.nc"
nc_data <- nc_open(file_path)

z1 <- ncvar_get(nc_data, "v")

nc_close(nc_data)

# Concatenate both periods along the time dimension
z <- abind::abind(z, z1, along = 3)

rm(z1)
gc()

# Keep data ordered as z[lon, lat, time]
z <- aperm(z, c(1, 2, 3))

print(dim(z))  # expected: [nlon, nlat, ndays]

################################
# Subset region
################################
lat_idx <- which(lat >= lat_min & lat <= lat_max)
lon_idx <- which(lon >= lon_min & lon <= lon_max)

lat <- lat[lat_idx]
lon <- lon[lon_idx]

z <- z[lon_idx, lat_idx, ]

################################
# Reshape into matrix [gridpoint × day]
################################
num_days <- dim(z)[3]
num_gridpoints <- length(lon) * length(lat)

z <- matrix(aperm(z, c(1, 2, 3)),
            nrow = num_gridpoints,
            ncol = num_days)

print(dim(z))  # expected: [num_gridpoints, num_days]

################################
# Daily climatology removal (for composites)
################################
fechas <- seq(as.Date("1979-01-01"), as.Date("2023-12-31"), by = "day")

df_z <- as.data.frame(t(z))
df_z$fecha <- fechas

df_z <- df_z %>%
  mutate(mes_dia = format(fecha, "%m-%d"))

promedios_diarios <- df_z %>%
  group_by(mes_dia) %>%
  summarise(across(starts_with("V"), mean))

anomalías <- df_z %>%
  left_join(promedios_diarios, by = "mes_dia", suffix = c("_original", "_promedio"))

promedio <- anomalías %>%
  select(contains("_promedio"))

anomalías <- anomalías %>%
  select(contains("_original"))

gc()

z_filtrado <- anomalías - promedio

################################
# Quick diagnostic plot (single gridpoint time series)
################################
n <- 100
g1 <- ggplot() +
  geom_line(aes(x = 1:length(promedio[, n]),  y = promedio[, n] / 10)) +
  geom_line(aes(x = 1:length(anomalías[, n]), y = anomalías[, n] / 10), col = "red", alpha = .5) +
  xlab("days") + ylab("z500")

g2 <- ggplot() +
  geom_line(aes(x = 1:length(z_filtrado[, n]), y = z_filtrado[, n] / 10), col = "blue") +
  xlab("days") + ylab("z500 filtered")

g <- grid.arrange(g1, g2)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")

ggsave(
  filename = paste("pruebas_series_", as.character(n), "_z500.png"),
  plot = g,
  device = "png",
  dpi = 300,
  width = 8,
  height = 5
)

rm(df_z, promedios_diarios, anomalías, nc_data)
gc()

# Convert filtered anomalies back to [gridpoint × day] matrix
zv_filtrado <- t(as.matrix(z_filtrado))
row.names(zv_filtrado) <- NULL
colnames(zv_filtrado) <- NULL

dim(zv_filtrado)

# Store v850 matrix (unfiltered)
zv <- z
rm(z)
gc()

################################
# Combine u and v matrices
################################
z <- rbind(zu, zv)
z_filtrado <- rbind(zu_filtrado, zv_filtrado)

################################
# Standardize input for PCA + k-means
################################
k <- 7  # number of clusters (final classification)

# Standardize features (columns) after transposing to [day × features]
Xs <- scale(t(z), center = TRUE, scale = TRUE)

##################################################### PCA + k-means

##################################################
# PCA (keep leading components, then cluster in PC space)
##################################################

# Correlation matrix and SVD
RR <- cor(Xs)
DD <- svd(RR)

U  <- DD$u
V  <- DD$v
dv <- DD$d
dd <- diag(DD$d)

# PC scores (standardized)
A  <- dv^(-0.5)
AA <- diag(A)
Zs <- Xs %*% U %*% AA  # standardized PC scores

# Explained variance (%)
SumaDiag <- sum(dv)
VAR <- (dv / SumaDiag) * 100

# Number of PCs retained (cumulative variance > 90%)
m <- which(cumsum(VAR) > 90)[1]
m

rm(U, V, RR, AA, DD, dd, A, dv)
gc()


###########################

# Silhouette index
silhouette_optimized_serial <- function(X, k, ns, iter.max) {
  nk <- length(k)
  silhouette_avg <- numeric(nk)
  
  silhouette_wrapper <- function(k_val) {
    set.seed(12345)
    kmeans_result <- kmeans(X, centers = k_val, nstart = ns, iter.max = iter.max, algorithm = "MacQueen")
    cluster_labels <- kmeans_result$cluster
    sil <- silhouette(cluster_labels, dist(X))
    return(mean(sil[, 3]))
  }
  
  silhouette_avg <- sapply(k, silhouette_wrapper)
  return(list(k = k, silhouette_avg = silhouette_avg))
}

# Compute silhouette (example: 2–20)
resultados_silhouette_serial <- silhouette_optimized_serial(Zs[, 1:m], k = 2:20, ns = 25, iter.max = 10000)

# Plot silhouette curve
plot(resultados_silhouette_serial$k, resultados_silhouette_serial$silhouette_avg, type = 'b',
     xlab = "Number of clusters (k)", ylab = "Mean silhouette",
     main = "")

# Save silhouette plot
dev.copy(png, filename = "Silhouette_PCA_k-means_UV.png", width = 500, height = 400)
dev.off()



############################# Apply k-means (final)

k <- 7       # select the number of clusters
set.seed(12345)

# k-means on retained PC scores
fit <- kmeans(Zs[, 1:m], centers = k, iter.max = 100000, algorithm = "MacQueen")

# Cluster label per day
vector_cluster <- as.numeric(fit$cluster)

# Histogram: number of days per CP
g <- ggplot() +
  geom_histogram(aes(x = vector_cluster), binwidth = 1, fill = "cyan", col = "black") +
  ylab("Frequency") + xlab("CP") + ggtitle("Days in each CP")

plot(g)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
ggsave(paste("CPs_PCA_Kmeans_", k, "UV_days_CPs.png"), plot = g, width = 4, height = 3)

################################
# Daily CP table + annual frequencies
################################

fechas <- seq(as.Date("1979-01-01"), as.Date("2023-12-31"), by = "day")

fechas <- data.frame(
  dia = as.numeric(format(fechas, "%d")),
  mes = as.numeric(format(fechas, "%m")),
  año = as.numeric(format(fechas, "%Y"))
)

fechas$CP <- vector_cluster
head(fechas)

# Export daily CP assignment
write.csv(fechas, paste("dias_clasificados_UV_PCA_k-means", k, ".csv", sep = ""), row.names = F)

# Annual counts per CP
fechas <- fechas %>%
  group_by(año, CP) %>%
  summarise(conteo = n()) %>%
  ungroup()

head(fechas)

# Plot annual frequency per CP
g <- ggplot(fechas, aes(x = año, y = conteo)) +
  geom_line() +
  facet_wrap(~ CP, scales = "free_y", ncol = 2) +
  labs(title = "", x = "year", y = "Frequency") +
  theme_bw()

plot(g)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
ggsave(paste("CPs_PCA_Kmeans", k, "_UV_frecuencia.png", sep = ""), plot = g, width = 8, height = 5)

# Export annual frequency table
write.csv(fechas, paste("CPs_PCA_Kmeans", k, "_UV_frecuencia.csv", sep = ""), row.names = F)

gc()

################################
# Monthly climatological frequency (%)
################################

fechas <- seq(as.Date("1979-01-01"), as.Date("2023-12-31"), by = "day")

fechas <- data.frame(
  dia = as.numeric(format(fechas, "%d")),
  mes = as.numeric(format(fechas, "%m")),
  año = as.numeric(format(fechas, "%Y"))
)

fechas$CP <- vector_cluster
head(fechas)

n_dias <- fechas %>%
  group_by(mes) %>%
  summarise(tot_dias = n()) %>%
  ungroup()

fechas <- fechas %>%
  group_by(mes, CP) %>%
  summarise(conteo = n()) %>%
  ungroup()

fechas <- merge(n_dias, fechas, by = "mes")

fechas <- fechas %>%
  group_by(mes, CP) %>%
  summarise(porcentaje = conteo / tot_dias * 100) %>%
  ungroup()

head(fechas)

g <- ggplot(fechas, aes(x = mes, y = porcentaje, fill = as.factor(CP))) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 1) +
  scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "", x = "Month", y = "Frequency (%)", fill = "CP") +
  theme_bw(base_size = 15) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.text = element_text(face = "bold"))

plot(g)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
ggsave(paste("CPs_PCA_Kmeans_", k, "UV_frecuencia_mensual.png", sep = ""), plot = g, width = 8, height = 5)

gc()

################################
# Daily climatological frequency (%)
################################

fechas <- seq(as.Date("1979-01-01"), as.Date("2023-12-31"), by = "day")

fechas <- data.frame(
  dia = as.numeric(format(fechas, "%d")),
  mes = as.numeric(format(fechas, "%m")),
  año = as.numeric(format(fechas, "%Y"))
)

fechas$CP <- vector_cluster
head(fechas)

n_dias <- fechas %>%
  group_by(mes, dia) %>%
  summarise(tot_dias = n()) %>%
  ungroup()

fechas <- fechas %>%
  group_by(mes, dia, CP) %>%
  summarise(conteo = n()) %>%
  ungroup()

fechas <- merge(n_dias, fechas, by = c("mes", "dia"))

fechas <- fechas %>%
  group_by(mes, dia, CP) %>%
  summarise(porcentaje = conteo / tot_dias * 100) %>%
  ungroup()

fechas <- fechas %>%
  mutate(fecha_ficticia = as.Date(paste("2000", mes, dia, sep = "-"), format = "%Y-%m-%d"))

head(fechas)

# Color palette for CPs (up to 12)
colores <- brewer.pal(min(k, 12), "Paired")

g <- ggplot(fechas, aes(x = fecha_ficticia, y = porcentaje, fill = as.factor(CP))) +
  geom_bar(stat = "identity", position = "stack", color = NA, width = 1) +
  scale_x_date(date_breaks = "1 month", date_labels = "%m", expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = colores) +
  labs(title = "", x = "Month", y = "Frequency (%)", fill = "CP") +
  theme_bw(base_size = 15) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.text = element_text(face = "bold"))

plot(g)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
write.csv(fechas, "frec_clim_CP_1979-2023.csv", row.names = F)
ggsave(paste("CPs_PCA_Kmeans_", k, "UV_frecuencia_diaria.png", sep = ""), plot = g, width = 6, height = 3.5)

gc()

################################
# Seasonal climatological frequency (%)
################################

fechas <- seq(as.Date("1979-01-01"), as.Date("2023-12-31"), by = "day")

fechas <- data.frame(
  dia = as.numeric(format(fechas, "%d")),
  mes = as.numeric(format(fechas, "%m")),
  año = as.numeric(format(fechas, "%Y")),
  fecha = fechas
)

# Season tag per day
fechas <- fechas %>%
  mutate(season = case_when(
    mes %in% c(12, 1, 2) ~ "DJF",
    mes %in% c(3, 4, 5) ~ "MAM",
    mes %in% c(6, 7, 8) ~ "JJA",
    mes %in% c(9, 10, 11) ~ "SON"
  ))

fechas$CP <- vector_cluster

n_dias_season <- fechas %>%
  group_by(season) %>%
  summarise(tot_dias = n()) %>%
  ungroup()

fechas <- fechas %>%
  group_by(season, CP) %>%
  summarise(conteo = n()) %>%
  ungroup()

# Save seasonal counts for later use in seasonal composites
cantidad <- fechas

fechas <- merge(n_dias_season, fechas, by = "season")

fechas <- fechas %>%
  group_by(season, CP) %>%
  summarise(porcentaje = conteo / tot_dias * 100) %>%
  ungroup()

fechas$season <- factor(fechas$season, levels = c("DJF", "MAM", "JJA", "SON"))

g <- ggplot(fechas, aes(x = season, y = porcentaje, fill = as.factor(CP))) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 1) +
  labs(title = "", x = "Season", y = "Frequency (%)", fill = "CP") +
  theme_bw(base_size = 15) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "right",
        strip.background = element_rect(fill = "lightgray", color = "black"),
        strip.text = element_text(face = "bold"))

plot(g)

setwd("C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures")
ggsave(paste("CPs_PCA_Kmeans_", k, "UV_frecuencia_seasonal.png", sep = ""), plot = g, width = 6, height = 6)

gc()

################################
# CP centroid patterns (annual, unfiltered fields)
################################

xmn <- min(lon); xmx <- max(lon)
ymn <- min(lat); ymx <- max(lat)

# Build lon/lat grid vectors (matching gridpoint indexing)
longitude <- rep(lon, times = length(lat))
latitude <- rep(lat, each = length(lon))

# Compute CP mean u and v (annual)
Zu <- list()
for (cp in 1:k) {
  Zu[[cp]] <- rowMeans(zu[, which(vector_cluster == cp)], na.rm = TRUE)
}

Zv <- list()
for (cp in 1:k) {
  Zv[[cp]] <- rowMeans(zv[, which(vector_cluster == cp)], na.rm = TRUE)
}

datitos <- data.frame(
  Longitude = rep(longitude, k),
  Latitude = rep(latitude, k),
  Zu = unlist(Zu),
  Zv = unlist(Zv),
  CP = rep(1:k, each = length(longitude))
)

# Vector magnitude and direction
datitos$m <- Mag(datitos$Zu, datitos$Zv)
datitos$ang <- Angle(datitos$Zu, datitos$Zv)

print(head(datitos))

world <- map_data("world")

g <- ggplot() +
  geom_arrow(data = datitos, aes(x = Longitude, y = Latitude, mag = m, angle = ang),
             skip = 2, show.legend = F, size = 0.4) +
  facet_wrap(~CP, nrow = 3) + xlab("") + ylab("") +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = .6) +
  coord_fixed(xlim = c(min(lon), max(lon)), ylim = c(min(lat), max(lat)), expand = 0) +
  theme(panel.background = element_blank())

plot(g)

################################
# CP centroid patterns (annual, filtered anomalies)
################################

xmn <- min(lon); xmx <- max(lon)
ymn <- min(lat); ymx <- max(lat)

longitude <- rep(lon, times = length(lat))
latitude <- rep(lat, each = length(lon))

Zu <- list()
for (cp in 1:k) {
  Zu[[cp]] <- rowMeans(zu_filtrado[, which(vector_cluster == cp)], na.rm = TRUE)
}

Zv <- list()
for (cp in 1:k) {
  Zv[[cp]] <- rowMeans(zv_filtrado[, which(vector_cluster == cp)], na.rm = TRUE)
}

datitos <- data.frame(
  Longitude = rep(longitude, k),
  Latitude = rep(latitude, k),
  Zu = unlist(Zu),
  Zv = unlist(Zv),
  CP = rep(1:k, each = length(longitude))
)

datitos$m <- Mag(datitos$Zu, datitos$Zv)
datitos$ang <- Angle(datitos$Zu, datitos$Zv)

print(head(datitos))

world <- map_data("world")

# Annual counts per CP (for facet labels)
cantidad_ann <- aggregate(conteo ~ CP, FUN = sum, cantidad)
etiquetas <- setNames(cantidad_ann$conteo, cantidad_ann$CP)

g <- ggplot() +
  geom_arrow(data = datitos,
             aes(x = Longitude, y = Latitude, mag = m, angle = ang),
             skip = 2, show.legend = TRUE, size = 0.3, col = "blue") +
  facet_wrap(~CP, nrow = 1,
             labeller = labeller(CP = function(x) paste0("CP ", x, " (", etiquetas[x], ")"))) +
  xlab("") + ylab("") + ggtitle("ANNUAL") +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = 0.3) +
  coord_fixed(xlim = c(min(lon), max(lon)), ylim = c(min(lat), max(lat)), expand = 0) +
  theme(panel.background = element_blank())

plot(g)

ggsave(
  filename = paste("CPs_PCA_Kmeans_UV_", as.character(k), "_CPs.png", sep = ""),
  plot = g,
  device = "png",
  dpi = 300,
  width = 8,
  height = 5
)

# Export CP mean patterns (annual, filtered)
write.csv(datitos, paste("Cps_patterns_PCA_kmeans_", k, ".csv"), row.names = F)



#################### Composites by season

fechas <- seq(as.Date("1979-01-01"), as.Date("2023-12-31"), by = "day")

fechas <- data.frame(
  dia = as.numeric(format(fechas, "%d")),
  mes = as.numeric(format(fechas, "%m")),
  año = as.numeric(format(fechas, "%Y")),
  fecha = fechas
)

### DJF

Zu <- list()
for (cp in 1:k) {
  Zu[[cp]] <- rowMeans(zu_filtrado[, which(vector_cluster == cp & fechas$mes %in% c(12,1,2))], na.rm = TRUE)
}

Zv <- list()
for (cp in 1:k) {
  Zv[[cp]] <- rowMeans(zv_filtrado[, which(vector_cluster == cp & fechas$mes %in% c(12,1,2))], na.rm = TRUE)
}

datitos <- data.frame(
  Longitude = rep(longitude, k),
  Latitude = rep(latitude, k),
  Zu = unlist(Zu),
  Zv = unlist(Zv),
  CP = rep(1:k, each = length(longitude))
)

datitos$m <- Mag(datitos$Zu, datitos$Zv)
datitos$ang <- Angle(datitos$Zu, datitos$Zv)

print(head(datitos))

world <- map_data("world")

# CP counts within DJF (for facet labels)
etiquetas1 <- setNames(cantidad$conteo[cantidad$season == "DJF"], cantidad$CP[cantidad$season == "DJF"])

g1 <- ggplot() +
  geom_arrow(data = datitos,
             aes(x = Longitude, y = Latitude, mag = m, angle = ang),
             skip = 3, show.legend = TRUE, size = 0.3, col = "blue") +
  facet_wrap(~CP, nrow = 1,
             labeller = labeller(CP = function(x) paste0("CP ", x, " (", etiquetas1[x], ")"))) +
  xlab("") + ylab("") + ggtitle("DJF") +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = 0.3) +
  coord_fixed(xlim = c(min(lon), max(lon)), ylim = c(min(lat), max(lat)), expand = 0) +
  theme(panel.background = element_blank())

plot(g1)

### MAM

Zu <- list()
for (cp in 1:k) {
  Zu[[cp]] <- rowMeans(zu_filtrado[, which(vector_cluster == cp & fechas$mes %in% c(3,4,5))], na.rm = TRUE)
}

Zv <- list()
for (cp in 1:k) {
  Zv[[cp]] <- rowMeans(zv_filtrado[, which(vector_cluster == cp & fechas$mes %in% c(3,4,5))], na.rm = TRUE)
}

datitos <- data.frame(
  Longitude = rep(longitude, k),
  Latitude = rep(latitude, k),
  Zu = unlist(Zu),
  Zv = unlist(Zv),
  CP = rep(1:k, each = length(longitude))
)

datitos$m <- Mag(datitos$Zu, datitos$Zv)
datitos$ang <- Angle(datitos$Zu, datitos$Zv)

print(head(datitos))

world <- map_data("world")

# CP counts within MAM (for facet labels)
etiquetas2 <- setNames(cantidad$conteo[cantidad$season == "MAM"], cantidad$CP[cantidad$season == "MAM"])

g2 <- ggplot() +
  geom_arrow(data = datitos,
             aes(x = Longitude, y = Latitude, mag = m, angle = ang),
             skip = 3, show.legend = TRUE, size = 0.3, col = "blue") +
  facet_wrap(~CP, nrow = 1,
             labeller = labeller(CP = function(x) paste0("CP ", x, " (", etiquetas2[x], ")"))) +
  xlab("") + ylab("") + ggtitle("MAM") +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = 0.3) +
  coord_fixed(xlim = c(min(lon), max(lon)), ylim = c(min(lat), max(lat)), expand = 0) +
  theme(panel.background = element_blank())

plot(g2)

### JJA

Zu <- list()
for (cp in 1:k) {
  cols <- which(vector_cluster == cp & fechas$mes %in% c(6,7,8))
  
  if (length(cols) > 1) {
    Zu[[cp]] <- rowMeans(zu_filtrado[, cols, drop = FALSE], na.rm = TRUE)
  } else {
    if (length(cols) == 0) {
      Zu[[cp]] <- rep(NA, nrow(zu_filtrado))
    } else {
      Zu[[cp]] <- zu_filtrado[, cols]
    }
  }
}

Zv <- list()
for (cp in 1:k) {
  cols <- which(vector_cluster == cp & fechas$mes %in% c(6,7,8))
  
  if (length(cols) > 1) {
    Zv[[cp]] <- rowMeans(zv_filtrado[, cols, drop = FALSE], na.rm = TRUE)
  } else {
    if (length(cols) == 0) {
      Zv[[cp]] <- rep(NA, nrow(zv_filtrado))
    } else {
      Zv[[cp]] <- zv_filtrado[, cols]
    }
  }
}

datitos <- data.frame(
  Longitude = rep(longitude, k),
  Latitude = rep(latitude, k),
  Zu = unlist(Zu),
  Zv = unlist(Zv),
  CP = rep(1:k, each = length(longitude))
)

datitos$m <- Mag(datitos$Zu, datitos$Zv)
datitos$ang <- Angle(datitos$Zu, datitos$Zv)

print(head(datitos))

world <- map_data("world")

# CP counts within JJA (for facet labels)
etiquetas3 <- setNames(cantidad$conteo[cantidad$season == "JJA"], cantidad$CP[cantidad$season == "JJA"])

g3 <- ggplot() +
  geom_arrow(data = datitos,
             aes(x = Longitude, y = Latitude, mag = m, angle = ang),
             skip = 3, show.legend = TRUE, size = 0.3, col = "blue") +
  facet_wrap(~CP, nrow = 1,
             labeller = labeller(CP = function(x) paste0("CP ", x, " (", etiquetas3[x], ")"))) +
  xlab("") + ylab("") + ggtitle("JJA") +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = 0.3) +
  coord_fixed(xlim = c(min(lon), max(lon)), ylim = c(min(lat), max(lat)), expand = 0) +
  theme(panel.background = element_blank())

plot(g3)

### SON

Zu <- list()
for (cp in 1:k) {
  Zu[[cp]] <- rowMeans(zu_filtrado[, which(vector_cluster == cp & fechas$mes %in% c(9,10,11))], na.rm = TRUE)
}

Zv <- list()
for (cp in 1:k) {
  Zv[[cp]] <- rowMeans(zv_filtrado[, which(vector_cluster == cp & fechas$mes %in% c(9,10,11))], na.rm = TRUE)
}

datitos <- data.frame(
  Longitude = rep(longitude, k),
  Latitude = rep(latitude, k),
  Zu = unlist(Zu),
  Zv = unlist(Zv),
  CP = rep(1:k, each = length(longitude))
)

datitos$m <- Mag(datitos$Zu, datitos$Zv)
datitos$ang <- Angle(datitos$Zu, datitos$Zv)

print(head(datitos))

world <- map_data("world")

# CP counts within SON (for facet labels)
etiquetas4 <- setNames(cantidad$conteo[cantidad$season == "SON"], cantidad$CP[cantidad$season == "SON"])

g4 <- ggplot() +
  geom_arrow(data = datitos,
             aes(x = Longitude, y = Latitude, mag = m, angle = ang),
             skip = 3, show.legend = TRUE, size = 0.3, col = "blue") +
  facet_wrap(~CP, nrow = 1,
             labeller = labeller(CP = function(x) paste0("CP ", x, " (", etiquetas4[x], ")"))) +
  xlab("") + ylab("") + ggtitle("SON") +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "black", size = 0.3) +
  coord_fixed(xlim = c(min(lon), max(lon)), ylim = c(min(lat), max(lat)), expand = 0) +
  theme(panel.background = element_blank())

plot(g4)

################################
# Combine seasonal plots + annual plot and export
################################
graficos <- list(g1, g2, g3, g4, g) %>%
  lapply(function(grafico) {
    grafico + theme(plot.margin = unit(c(-1, 0, -1, -0.5), "cm"))
  })

g_f <- grid.arrange(grobs = graficos, nrow = 5, padding = unit(0, "line"))

ggsave(
  filename = paste("CPs_PCA_Kmeans_UV_", as.character(k), "_CPs_seasons.png", sep = ""),
  plot = g_f,
  device = "png",
  dpi = 300,
  width = 10/7*k,
  height = 9
)





################################
# Optional: within-CP similarity diagnostics (RHO + distance) using filtered data
################################

# CP mean pattern in filtered space
Z <- list()
for (cp in 1:k) {
  Z[[cp]] <- rowMeans(z_filtrado[, which(vector_cluster == cp)], na.rm = TRUE)
}

# Per-day correlation and Euclidean distance to its CP mean
RHO_df <- data.frame(
  CP   = as.character(vector_cluster),
  RHO  = vector(length = length(vector_cluster)),
  dist = vector(length = length(vector_cluster))
)

for (i in 1:nrow(RHO_df)) {
  RHO_df$RHO[i] <- cor(z_filtrado[, i], Z[[vector_cluster[i]]])
}

for (i in 1:nrow(RHO_df)) {
  RHO_df$dist[i] <- sqrt(sum((z_filtrado[, i] - Z[[vector_cluster[i]]])^2))
}

head(RHO_df)

# Normalize distance to [0,1] for plotting
escala <- max(RHO_df$dist)
RHO_df$dist <- RHO_df$dist / escala

# Long format for boxplots
RHO_df <- gather(RHO_df, metric, valor, c(2, 3))
RHO_df$metric <- factor(RHO_df$metric, c("RHO", "dist"))

g <- ggplot(RHO_df, aes(x = CP, y = valor)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, outlier.size = 2, aes(fill = metric)) +
  scale_fill_manual(values = c("white", "cyan")) +
  geom_hline(yintercept = 0, size = 2) +
  labs(
    title = "",
    x = "CP",
    y = "RHO",
    y.sec = "Distance"
  ) +
  scale_y_continuous(
    breaks = seq(-1, 1, 0.2),
    limits = c(-1, 1),
    sec.axis = sec_axis(~ . * escala, name = "Distance")
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

plot(g)

ggsave(
  filename = paste("RHO_dist_filtered_PCA_Kmeans_UV_", as.character(k), "_CPs.png", sep = ""),
  plot = g,
  device = "png",
  dpi = 300,
  width = 8,
  height = 4
)
