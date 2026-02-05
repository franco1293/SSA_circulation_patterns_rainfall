# =============================================================================
# SACZ (Vera et al. definition using OLR anomalies) + ENSO phases vs CPs
#
# =============================================================================

# ----------------------------
# Libraries (keep only what is used)
# ----------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(zoo)
library(lubridate)
library(stringr)
library(ggplot2)
library(ncdf4)   

# ----------------------------
# User parameters (EDIT)
# ----------------------------
k <- 7
season_name <- "DJF"          # "ONDJFMA", "DJF", "MAM", "JJA", "SON"
n_mc <- 10000
set.seed(1234)

base_dir <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st"
fig_dir  <- file.path(base_dir, "Figures")

nino_file <- file.path(base_dir, "nino34.xlsx")
cp_file   <- file.path(fig_dir, paste0("dias_clasificados_UV_PCA_k-means", k, ".csv"))

# OLR file template: OLR_NCEP1_Vera_<periodo>.nc
olr_template <- file.path(base_dir, "OLR_NCEP1_Vera_%s.nc")
olr_periods  <- c("1979-1989", "1990-1999", "2000-2009", "2010-2023")

# ----------------------------
# Season definitions
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
# 1) ENSO phases (Niño3.4 + NOAA operational criterion)
# ----------------------------
nino34 <- read_excel(nino_file) %>%
  pivot_longer(cols = `1`:`12`, names_to = "mes", values_to = "nino34") %>%
  mutate(
    mes = as.integer(mes),
    fecha = as.Date(paste(año, mes, "01", sep = "-"))
  ) %>%
  arrange(fecha) %>%
  mutate(
    nino34_3mes = rollmean(nino34, 3, fill = NA, align = "right"),
    fase_raw = case_when(
      nino34_3mes >= 0.5  ~ "El Niño",
      nino34_3mes <= -0.5 ~ "La Niña",
      TRUE                ~ "Neutral"
    )
  )

rle_fase <- rle(nino34$fase_raw)
fases_validadas <- ifelse(
  rle_fase$values %in% c("El Niño", "La Niña") & rle_fase$lengths >= 5,
  rle_fase$values,
  "Neutral"
)
nino34$fase_noaa <- inverse.rle(list(lengths = rle_fase$lengths, values = fases_validadas))
nino34 <- nino34[nino34$año %in% 1979:2023, c("año", "mes", "fase_noaa")]

# ----------------------------
# 2) OLR extraction + SACZ classification (Vera-style using OLR anomaly sign)
# ----------------------------
olr_extraction <- function(periodo) {
  archivo_nc <- sprintf(olr_template, periodo)
  
  nc <- nc_open(archivo_nc)
  
  # Variable name assumed "ulwrf" (edit if needed)
  olr  <- ncvar_get(nc, "ulwrf")
  lon  <- ncvar_get(nc, "lon")  # kept for completeness (not used below)
  lat  <- ncvar_get(nc, "lat")  # kept for completeness (not used below)
  time <- ncvar_get(nc, "time")
  
  nc_close(nc)
  
  # Reorder dims to [time, lon, lat] (original assumed [lon, lat, time])
  olr <- aperm(olr, c(3, 1, 2))
  
  # Domain-mean OLR per day (faster than flattening to a huge DF)
  media <- apply(olr, 1, mean, na.rm = TRUE)
  
  # Convert time to Date (hours since 1800-01-01 in your files)
  fechas <- as.POSIXct("1800-01-01 00:00:00", tz = "UTC") + time * 3600
  fechas <- as.Date(fechas)
  
  serie_olr <- data.frame(
    año = as.integer(format(fechas, "%Y")),
    mes = as.integer(format(fechas, "%m")),
    dia = as.integer(format(fechas, "%d")),
    media = media
  )
  
  # Calendar-day climatology (MM-DD) and daily anomaly
  serie_olr$mes_dia <- paste0(serie_olr$mes, "-", serie_olr$dia)
  clim <- aggregate(media ~ mes_dia, data = serie_olr, FUN = mean)
  names(clim) <- c("mes_dia", "media_mes_dia")
  
  serie_olr <- merge(serie_olr, clim, by = "mes_dia", all.x = TRUE)
  serie_olr$indice <- serie_olr$media - serie_olr$media_mes_dia
  
  # SACZ classification
  serie_olr$condicion <- "no SACZ"
  serie_olr$condicion[serie_olr$indice < 0] <- "SACZ event"
  
  serie_olr[, c("año", "mes", "dia", "media", "indice", "condicion")]
}

serie_olr <- do.call(rbind, lapply(olr_periods, olr_extraction))

# ----------------------------
# 3) Merge SACZ + CPs + ENSO
# ----------------------------
cp <- read.csv(cp_file)
stopifnot(all(c("año", "mes", "dia", "CP") %in% names(cp)))

x <- merge(serie_olr, cp, by = c("año", "mes", "dia"))
x <- x[, c("dia", "mes", "año", "CP", "condicion")]

x <- merge(x, nino34, by = c("año", "mes"))
names(x)[names(x) == "condicion"] <- "SACZ_sign"

# Combined state: SACZ_ENSO
x$SACZ_ENSO <- paste(x$SACZ_sign, x$fase_noaa, sep = "-")
SACZ_ENSO_comb <- sort(unique(x$SACZ_ENSO))

# Seasonal pool for Monte Carlo
x_season <- x[x$mes %in% season_months, ]

# ----------------------------
# Helper: Monte Carlo
# ----------------------------
mc_p_state_given_CP <- function(state_label) {
  # p(state | CP) = N(CP,state) / N(CP)
  do.call(rbind, lapply(1:k, function(cp_i) {
    datos_cp <- x[x$CP == cp_i & x$mes %in% season_months, ]
    N <- nrow(datos_cp)
    if (N == 0) return(data.frame(CP = cp_i, p = NA, p95 = NA, p025 = NA, sign = "no", condicion = state_label))
    
    p_obs <- mean(datos_cp$SACZ_ENSO == state_label) * 100
    
    mc_vals <- replicate(n_mc, {
      samp <- sample(x_season$SACZ_ENSO, size = N, replace = FALSE)
      mean(samp == state_label) * 100
    })
    
    p95  <- as.numeric(quantile(mc_vals, 0.975, na.rm = TRUE))
    p025 <- as.numeric(quantile(mc_vals, 0.025, na.rm = TRUE))
    
    sign <- "no"
    if (p_obs > p95)  sign <- "superior"
    if (p_obs < p025) sign <- "inferior"
    
    data.frame(CP = cp_i,
               p = round(p_obs, 1),
               p95 = round(p95, 1),
               p025 = round(p025, 1),
               sign = sign,
               condicion = state_label)
  }))
}

mc_p_CP_given_state <- function(state_label, x_pool) {
  # p(CP | state) = N(CP,state) / N(state)
  datos_state <- x_pool[x_pool$SACZ_ENSO == state_label, ]
  N <- nrow(datos_state)
  
  do.call(rbind, lapply(1:k, function(cp_i) {
    if (N == 0) return(data.frame(CP = cp_i, p = NA, p95 = NA, p025 = NA, sign = "no", condicion = state_label))
    
    p_obs <- mean(datos_state$CP == cp_i) * 100
    
    mc_vals <- replicate(n_mc, {
      samp <- sample(x_pool$CP, size = N, replace = FALSE)
      mean(samp == cp_i) * 100
    })
    
    p95  <- as.numeric(quantile(mc_vals, 0.975, na.rm = TRUE))
    p025 <- as.numeric(quantile(mc_vals, 0.025, na.rm = TRUE))
    
    sign <- "no"
    if (p_obs > p95)  sign <- "superior"
    if (p_obs < p025) sign <- "inferior"
    
    data.frame(CP = cp_i,
               p = round(p_obs, 1),
               p95 = round(p95, 1),
               p025 = round(p025, 1),
               sign = sign,
               condicion = state_label)
  }))
}

# =============================================================================
# A) p(SACZ_ENSO | CP) for selected season
# =============================================================================
df_A <- do.call(rbind, lapply(SACZ_ENSO_comb, mc_p_state_given_CP))

parts <- str_split_fixed(df_A$condicion, "-", 2)
df_A$SACZ <- parts[, 1]
df_A$nino <- parts[, 2]
df_A$SACZ[df_A$SACZ == "SACZ event"] <- "SACZ"

df_A$p <- round(df_A$p, 0)
df_A$p[df_A$sign == "no"] <- NA

g_A <- ggplot(df_A, aes(x = SACZ, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = na.omit(df_A), aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("superior" = "red", "inferior" = "cyan", "no" = NA), name = "") +
  geom_text(aes(label = p), size = 3) +
  facet_wrap(~ nino, ncol = 3) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "", y = "CP", title = "p(SACZ, ENSO | CP) = N(CP, SACZ, ENSO) / N(CP)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 9),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_A)

ggsave(
  filename = file.path(fig_dir, "SACZ-ENSO_porc_por_CP_def_Vera.png"),
  plot = g_A, dpi = 300, width = 4, height = 3, bg = "white"
)

# =============================================================================
# B) p(CP | SACZ_ENSO) for selected season
# =============================================================================
df_B <- do.call(rbind, lapply(SACZ_ENSO_comb, function(st) mc_p_CP_given_state(st, x_season)))

parts <- str_split_fixed(df_B$condicion, "-", 2)
df_B$SACZ <- parts[, 1]
df_B$nino <- parts[, 2]
df_B$SACZ[df_B$SACZ == "SACZ event"] <- "SACZ"

df_B$p <- round(df_B$p, 0)
df_B$p[df_B$sign == "no"] <- NA

g_B <- ggplot(df_B, aes(x = SACZ, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = na.omit(df_B), aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("superior" = "red", "inferior" = "cyan", "no" = NA), name = "") +
  geom_text(aes(label = p), size = 3) +
  facet_wrap(~ nino, ncol = 3) +
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
  filename = file.path(fig_dir, "al_rev_SACZ-ENSO_porc_por_CP_def_Vera.png"),
  plot = g_B, dpi = 300, width = 4, height = 3, bg = "white"
)

# =============================================================================
# 7) Alternative multi-season figure (example: DJF + MAM + SON, non-neutral only)
# =============================================================================
figura <- function(seas) {
  seas_months <- season_matrix[[seas]]
  x_s <- x[x$mes %in% seas_months, ]
  x_pool <- x_s  # MC pool restricted to season
  
  df <- do.call(rbind, lapply(SACZ_ENSO_comb, function(st) mc_p_CP_given_state(st, x_pool)))
  parts <- str_split_fixed(df$condicion, "-", 2)
  df$SACZ <- parts[, 1]
  df$nino <- parts[, 2]
  df$SACZ[df$SACZ == "SACZ event"] <- "SACZ"
  df$p <- round(df$p, 0)
  df$season <- seas
  df
}

df_alt <- rbind(figura("DJF"), figura("MAM"), figura("SON"))

df_alt <- df_alt[df_alt$nino != "Neutral", ]
df_alt$condicion <- gsub("SACZ event", "SACZ", df_alt$condicion)
df_alt$condicion <- gsub("no SACZ", "non-SACZ", df_alt$condicion)

g_alt <- ggplot(df_alt, aes(x = season, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_alt[df_alt$sign != "no", ], aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  facet_wrap(~ condicion, ncol = 2, scales = "free_x") +
  scale_color_manual(values = c("superior" = "red", "inferior" = "cyan", "no" = NA), name = "") +
  geom_text(data = df_alt[df_alt$sign != "no", ], aes(label = round(p, 1)), size = 3) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "", y = "CP", title = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_alt)

ggsave(
  filename = file.path(fig_dir, "al_rev_SACZ-ENSO_porc_por_CP_def_Vera.png"),
  plot = g_alt, dpi = 300, width = 2.8, height = 5, bg = "white"
)
