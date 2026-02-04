############################################################
# CP transitions and persistence analysis
############################################################

rm(list=ls()) ; graphics.off() ; gc()

library(ggpubr)
library(scales)
library(dplyr)
library(tidyr)
library(boot)
library(ncdf4)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(RColorBrewer)
library(lubridate)
library(abind)
library(data.table)
library(stringr)
library(readxl)
library(corrplot)
library(writexl)

#####################################################################################################
# 0) Paths / inputs
#####################################################################################################

# Daily CP classification file (CPs derived from 850-hPa winds via PCA + k-means)
file_cps <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures/dias_clasificados_UV_PCA_k-means7.csv"

# Output directory for figures
out_dir <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st/Figures"

#####################################################################################################
# 1) ANNUAL transition matrix + Monte Carlo (shuffle whole CP sequence)
#####################################################################################################

# Read CP classification (daily)
CPs <- read.csv(file_cps, sep="")

# Add season label (DJF/MAM/JJA/SON) based on month
CPs$season <- "DJF"
CPs$season[CPs$mes %in% c(3,4,5)] <- "MAM"
CPs$season[CPs$mes %in% c(6,7,8)] <- "JJA"
CPs$season[CPs$mes %in% c(9,10,11)] <- "SON"

# Ensure chronological order (needed for lead/lag transitions)
CPs <- CPs %>% arrange(año, mes, dia)

# Build "next-day CP" to count transitions CP(t) -> CP(t+1)
CPs <- CPs %>%
  mutate(CP_siguiente = lead(CP)) %>%
  filter(!is.na(CP_siguiente))  # drop last day (no next day)

# Observed transition counts
conteo_transiciones <- CPs %>%
  group_by(CP, CP_siguiente) %>%
  summarise(n = n(), .groups = "drop")

# Observed transition probabilities (row-normalized, %)
matrices_transicion <- conteo_transiciones %>%
  group_by(CP) %>%
  mutate(probabilidad = n / sum(n) * 100) %>%
  ungroup()

# All possible CP -> CP_siguiente combinations (used to fill missing transitions with 0)
combinaciones_posibles <- expand.grid(
  CP = unique(CPs$CP),
  season = unique(CPs$season),        # not used in annual case, but kept to match your original structure
  CP_siguiente = unique(CPs$CP)
)

# Monte Carlo: shuffle the entire CP sequence, recompute transitions, store probabilities
n_iteraciones <- 10000
set.seed(123)

resultados <- vector("list", n_iteraciones)

for (i in 1:n_iteraciones) {
  # Random CP order (destroys temporal dependence; preserves marginal CP frequencies)
  CPs_aleat <- data.frame(CP = sample(CPs$CP))
  
  # Transition counts in randomized sequence
  CPs_aleat <- CPs_aleat %>%
    mutate(CP_siguiente = lead(CP)) %>%
    filter(!is.na(CP_siguiente))
  
  conteo_tran_alea <- CPs_aleat %>%
    group_by(CP, CP_siguiente) %>%
    summarise(n = n(), .groups = "drop")
  
  matrices_tran_alea <- conteo_tran_alea %>%
    group_by(CP) %>%
    mutate(probabilidad = n / sum(n) * 100) %>%
    ungroup()
  
  # Fill missing transitions with 0 to keep a full matrix each iteration
  x <- merge(matrices_tran_alea, combinaciones_posibles, by = c("CP", "CP_siguiente"), all = TRUE)
  x[is.na(x)] <- 0
  
  resultados[[i]] <- x %>% mutate(iteracion = i)
}

resultados_df <- bind_rows(resultados)

# Monte Carlo threshold: 90th percentile of randomized probabilities, per transition
p90 <- resultados_df %>%
  group_by(CP, CP_siguiente) %>%
  summarise(percentil_90 = quantile(probabilidad, 0.90, na.rm = TRUE), .groups = "drop")

# Merge observed probabilities with p90 threshold
matrices_transicion <- merge(matrices_transicion, p90, by=c("CP", "CP_siguiente"), all=TRUE)
matrices_transicion[is.na(matrices_transicion)] <- 0

# Mark "significant" transitions: observed prob > p90, otherwise masked
matrices_transicion$sign <- "NO"
matrices_transicion$sign[matrices_transicion$probabilidad > matrices_transicion$percentil_90] <- "SI"

# Mask non-significant transitions in the plotted matrix
matrices_transicion$probabilidad[matrices_transicion$sign=="NO"] <- 0
matrices_transicion$probabilidad <- round(matrices_transicion$probabilidad, 0)
matrices_transicion$probabilidad[matrices_transicion$probabilidad == 0] <- NA

# Plot: annual transition matrix (only significant cells show value/color)
g <- ggplot(matrices_transicion, aes(x = CP_siguiente, y = CP, fill = probabilidad)) +
  geom_tile(col="black") +
  geom_text(aes(label = probabilidad), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "red", na.value="white") +
  scale_x_continuous(breaks = sort(unique(matrices_transicion$CP_siguiente)), expand = c(0,0)) +
  scale_y_continuous(breaks = sort(unique(matrices_transicion$CP)), expand = c(0,0)) +
  labs(title = "", x = "To", y = "From", fill = "%") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

plot(g)

setwd(out_dir)
ggsave("CP_Transitions_ANNUAL.png", plot = g, width = 4, height = 3)

#####################################################################################################
# 2) SEASONAL transition matrices + Monte Carlo (shuffle within each season)
#####################################################################################################

CPs <- read.csv(file_cps, sep="")

# Add season label
CPs$season <- "DJF"
CPs$season[CPs$mes %in% c(3,4,5)] <- "MAM"
CPs$season[CPs$mes %in% c(6,7,8)] <- "JJA"
CPs$season[CPs$mes %in% c(9,10,11)] <- "SON"

# Chronological order for transitions
CPs <- CPs %>% arrange(año, mes, dia)

# Next-day CP (global lead; season is used later for grouping/transitions)
CPs <- CPs %>%
  mutate(CP_siguiente = lead(CP)) %>%
  filter(!is.na(CP_siguiente))

# Observed seasonal transition counts
conteo_transiciones <- CPs %>%
  group_by(CP, CP_siguiente, season) %>%
  summarise(n = n(), .groups = "drop")

# Observed seasonal probabilities (row-normalized within each season)
matrices_transicion <- conteo_transiciones %>%
  group_by(CP, season) %>%
  mutate(probabilidad = n / sum(n) * 100) %>%
  ungroup()

# Full grid of season x CP x CP_siguiente combinations (to fill missing transitions)
combinaciones_posibles <- expand.grid(
  CP = unique(CPs$CP),
  season = unique(CPs$season),
  CP_siguiente = unique(CPs$CP)
)

# Monte Carlo: shuffle CP labels within each season, recompute transitions
n_iteraciones <- 10000
set.seed(123)

resultados <- vector("list", n_iteraciones)

for (i in 1:n_iteraciones) {
  # Shuffle CP within each season (preserves CP frequencies per season)
  CPs_aleat <- CPs %>%
    group_by(season) %>%
    mutate(CP = sample(CP)) %>%
    ungroup()
  
  # Recompute next-day CP within each season (transitions constrained to season blocks)
  CPs_aleat <- CPs_aleat %>%
    group_by(season) %>%
    mutate(CP_siguiente = lead(CP)) %>%
    filter(!is.na(CP_siguiente)) %>%
    ungroup()
  
  conteo_tran_alea <- CPs_aleat %>%
    group_by(season, CP, CP_siguiente) %>%
    summarise(n = n(), .groups = "drop")
  
  matrices_tran_alea <- conteo_tran_alea %>%
    group_by(season, CP) %>%
    mutate(probabilidad = n / sum(n) * 100) %>%
    ungroup()
  
  x <- merge(matrices_tran_alea, combinaciones_posibles,
             by = c("season", "CP", "CP_siguiente"), all = TRUE)
  x[is.na(x)] <- 0
  
  resultados[[i]] <- x %>% mutate(iteracion = i)
}

resultados_df <- bind_rows(resultados)

# Seasonal Monte Carlo threshold: p90 per (season, CP, CP_siguiente)
p90 <- resultados_df %>%
  group_by(season, CP, CP_siguiente) %>%
  summarise(percentil_90 = quantile(probabilidad, 0.90, na.rm = TRUE), .groups = "drop")

# Merge observed with thresholds
matrices_transicion <- merge(matrices_transicion, p90, by = c("season", "CP", "CP_siguiente"), all=TRUE)
matrices_transicion[is.na(matrices_transicion)] <- 0

# Significance flag
matrices_transicion$sign <- "NO"
matrices_transicion$sign[matrices_transicion$probabilidad > matrices_transicion$percentil_90] <- "SI"

# Prepare for plotting (show only significant values as text; NA elsewhere)
matrices_transicion$probabilidad <- round(matrices_transicion$probabilidad, 0)
matrices_transicion$probabilidad[matrices_transicion$probabilidad == 0] <- NA
matrices_transicion$season <- factor(matrices_transicion$season, levels = c("DJF", "MAM", "JJA", "SON"))

# Plot seasonal transition matrices (labels only for significant cells)
g <- ggplot(matrices_transicion, aes(x = CP_siguiente, y = CP, fill = probabilidad)) +
  geom_tile(col = "black") +
  geom_text(
    data = matrices_transicion[matrices_transicion$sign == "SI", ],
    aes(label = probabilidad),
    color = "black", size = 5
  ) +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  scale_x_continuous(breaks = sort(unique(matrices_transicion$CP_siguiente)), expand = c(0,0)) +
  scale_y_continuous(breaks = sort(unique(matrices_transicion$CP)), expand = c(0,0)) +
  facet_wrap(~ season) +
  labs(title = "", x = "To", y = "From", fill = "%") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))

plot(g)

setwd(out_dir)
ggsave("CP_Transitions_SEASONAL.png", plot = g, width = 6, height = 5)

#####################################################################################################
# 3) CP persistence by season (run-length statistics)
#####################################################################################################

CPs <- read.csv(file_cps, sep="")

# Add season label
CPs$season <- "DJF"
CPs$season[CPs$mes %in% c(3,4,5)] <- "MAM"
CPs$season[CPs$mes %in% c(6,7,8)] <- "JJA"
CPs$season[CPs$mes %in% c(9,10,11)] <- "SON"

# Sort by season + time (persistence is computed within each season block)
CPs <- CPs %>% arrange(season, año, mes, dia)

# Identify CP runs within each season (consecutive days with same CP)
CPs <- CPs %>%
  group_by(season) %>%
  mutate(CP_anterior = lag(CP)) %>%
  drop_na(CP, CP_anterior) %>%
  mutate(nueva_secuencia = ifelse(CP != CP_anterior, 1, 0)) %>%
  mutate(secuencia_id = cumsum(nueva_secuencia)) %>%
  ungroup()

# Run-length (duration) per CP and season
persistencia <- CPs %>%
  group_by(season, secuencia_id, CP) %>%
  summarise(duracion = n(), .groups = "drop") %>%
  ungroup()

# Persistence summary statistics per CP and season
estadisticas_persistencia <- persistencia %>%
  group_by(season, CP) %>%
  summarise(
    promedio = mean(duracion),
    mediana  = median(duracion),
    maximo   = max(duracion),
    minimo   = min(duracion),
    .groups  = "drop"
  )

print(estadisticas_persistencia)

estadisticas_persistencia$season <- factor(estadisticas_persistencia$season, levels = c("DJF", "MAM", "JJA", "SON"))

# Heatmap: maximum persistence
g1 <- ggplot(estadisticas_persistencia, aes(y = CP, x = season, fill = maximo)) +
  geom_tile(col="black") +
  geom_text(aes(label = maximo), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_y_continuous(breaks = 1:7, expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(title = "Maximum persistence", y = "CP", x = "", fill = "days") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank())

# Heatmap: mean persistence
g2 <- ggplot(estadisticas_persistencia, aes(y = CP, x = season, fill = round(promedio,1))) +
  geom_tile(col="black") +
  geom_text(aes(label = round(promedio,1)), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_y_continuous(breaks = 1:7, expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(title = "Mean persistence", y = "CP", x = "", fill = "days") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank())

g <- grid.arrange(g2, g1, ncol=2)

setwd(out_dir)
ggsave("CP_persistence.png", plot = g, width = 6, height = 3)
