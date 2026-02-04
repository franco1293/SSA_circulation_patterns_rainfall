# =============================================================================
# Trends in CP frequency by season and annual
# =============================================================================

# ----------------------------
# Libraries
# ----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(trend)     
library(scales)    
library(grid)      

# ----------------------------
# User parameters (EDIT)
# ----------------------------
tecnica <- "UV_PCA_k-means7"   
año1 <- 1979
año2 <- 2019

in_dir  <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures"
out_dir <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures"

in_file <- file.path(in_dir, paste0("dias_clasificados_", tecnica, ".csv"))

# ----------------------------
# 1) Read CP classifications and define seasons + hydro-year
# ----------------------------
CPs <- read.csv(in_file)

# Basic checks (fail fast)
needed_cols <- c("año", "mes", "dia", "CP")
stopifnot(all(needed_cols %in% names(CPs)))

CPs <- CPs %>%
  mutate(
    season = case_when(
      mes %in% c(12, 1, 2)  ~ "DJF",
      mes %in% c(3, 4, 5)   ~ "MAM",
      mes %in% c(6, 7, 8)   ~ "JJA",
      mes %in% c(9, 10, 11) ~ "SON",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(season)) %>%
  filter(año %in% año1:año2)

# Hydro-year for DJF (Jul-Dec -> next year)
CPs <- CPs %>%
  mutate(
    año_hidro = if_else(mes >= 7, año + 1L, año)
  )

# ----------------------------
# 2) Seasonal frequencies
#    - DJF: by hydro-year
#    - MAM/JJA/SON: by calendar year
# ----------------------------

# DJF (hydro-year)
cps_frec_hidro <- CPs %>%
  filter(season == "DJF") %>%
  group_by(season, año_hidro, CP) %>%
  summarise(N = n(), .groups = "drop")

# Fill missing (año_hidro x CP) with 0
todos_los_años_hidro <- sort(unique(CPs$año_hidro))
todos_los_CPs <- sort(unique(CPs$CP))

cps_frec_hidro <- expand.grid(
  año_hidro = todos_los_años_hidro,
  CP = todos_los_CPs,
  season = "DJF"
) %>%
  as_tibble() %>%
  left_join(cps_frec_hidro, by = c("season", "año_hidro", "CP")) %>%
  mutate(N = ifelse(is.na(N), 0, N))

# Keep central years only (matches your original trimming)
# (Because DJF hydro-year needs buffers at ends)
cps_frec_hidro <- cps_frec_hidro %>%
  filter(año_hidro %in% (año1 + 1):(año2 - 1))

# MAM/JJA/SON (calendar year)
cps_frec <- CPs %>%
  filter(season != "DJF") %>%
  group_by(season, año, CP) %>%
  summarise(N = n(), .groups = "drop")

# Fill missing (año x CP x season) with 0 for MAM/JJA/SON
todos_los_años_cal <- sort(unique(CPs$año))

cps_frec <- expand.grid(
  año = todos_los_años_cal,
  CP = todos_los_CPs,
  season = c("MAM", "JJA", "SON")
) %>%
  as_tibble() %>%
  left_join(cps_frec, by = c("season", "año", "CP")) %>%
  mutate(N = ifelse(is.na(N), 0, N))

# ----------------------------
# 3) Plot seasonal time series by CP (DJF plotted with año_hidro)
# ----------------------------
g_ts <- ggplot() +
  geom_line(data = cps_frec, aes(x = año, y = N, col = season)) +
  geom_line(data = cps_frec_hidro, aes(x = año_hidro, y = N, col = season), show.legend = TRUE) +
  facet_wrap(~ CP, scales = "free") +
  scale_color_manual(
    name = "Season",
    values = c("DJF" = "black", "MAM" = "orange", "JJA" = "blue", "SON" = "green"),
    breaks = c("DJF", "MAM", "JJA", "SON")
  ) +
  labs(x = "Year", y = "Frequency (days)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

plot(g_ts)

# Optional save:
# out_ts <- file.path(out_dir, paste0("CPs_seasonal_frequency_", tecnica, ".png"))
# ggsave(filename = out_ts, plot = g_ts, height = 4, width = 8, bg = "white", dpi = 300)

# ----------------------------
# 4) Sen's slope trends by season
# ----------------------------
# DJF (hydro-year) trends
trend_djf <- cps_frec_hidro %>%
  group_by(season, CP) %>%
  summarise(
    sen_slope = round(trend::sens.slope(N)$estimates, 3),
    p_value   = round(trend::sens.slope(N)$p.value, 2),
    .groups = "drop"
  )

# MAM/JJA/SON trends
trend_seas <- cps_frec %>%
  group_by(season, CP) %>%
  summarise(
    sen_slope = round(trend::sens.slope(N)$estimates, 3),
    p_value   = round(trend::sens.slope(N)$p.value, 2),
    .groups = "drop"
  )

trend_seasonal <- bind_rows(trend_seas, trend_djf) %>%
  mutate(asterisk = ifelse(p_value <= 0.05, "*", ""))

# ----------------------------
# 5) Annual trends (calendar year)
# ----------------------------
cps_frec_anual <- CPs %>%
  group_by(año, CP) %>%
  summarise(N = n(), .groups = "drop")

trend_anual <- cps_frec_anual %>%
  group_by(CP) %>%
  summarise(
    sen_slope = round(trend::sens.slope(N)$estimates, 3),
    p_value   = round(trend::sens.slope(N)$p.value, 2),
    .groups = "drop"
  ) %>%
  mutate(
    season = "ANNUAL",
    asterisk = ifelse(p_value <= 0.05, "*", "")
  )

# ----------------------------
# 6) Heatmap of trends (days/10yr)
# ----------------------------
trend_all <- bind_rows(trend_seasonal, trend_anual) %>%
  mutate(
    season = factor(season, levels = c("ANNUAL", "DJF", "MAM", "JJA", "SON")),
    CP = as.factor(CP)
  )

g_hm <- ggplot(trend_all, aes(x = season, y = CP, fill = sen_slope * 10)) +
  geom_tile(color = "black", show.legend = TRUE) +
  scale_fill_gradient2(
    low = "deepskyblue", mid = "white", high = "brown1",
    midpoint = 0, na.value = "white",
    name = "[days/10yr]"
  ) +
  geom_text(
    data = trend_all %>% filter(p_value <= 0.05),
    aes(label = round(sen_slope * 10, 1)),
    color = "black", size = 2
  ) +
  labs(x = "", y = "CP") +
  theme_minimal(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(0.3, "cm")
  )

plot(g_hm)

out_hm <- file.path(out_dir, paste0("trend_CPs_seasonal_", tecnica, ".png"))
ggsave(filename = out_hm, plot = g_hm, height = 2, width = 2, bg = "white", dpi = 300)
