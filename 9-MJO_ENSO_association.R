# =============================================================================
# MJO + ENSO phases and their association with Circulation Patterns (CPs)
# Notes / assumptions
#   - ENSO criterion uses 3-month running mean, aligned to the RIGHT (month t is
#     mean(t-2..t)), and NOAA “>=5 consecutive months” applied to the classified
#     phases (El Niño / La Niña). Everything else -> Neutral.
# =============================================================================

# ----------------------------
# Libraries (keep only what is used)
# ----------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(zoo)
library(stringr)
library(ggplot2)

# ----------------------------
# User parameters (EDIT)
# ----------------------------
k <- 7
season_name <- "SON"          # "ONDJFMA", "DJF", "MAM", "JJA", "SON"
n_mc <- 10000
set.seed(1234)

base_dir <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st"
fig_dir  <- file.path(base_dir, "Figures")

nino_file <- file.path(base_dir, "nino34.xlsx")
mjo_int_file <- file.path(base_dir, "MJO_intensity.xlsx")
mjo_ph_file  <- file.path(base_dir, "MJO_phase.xlsx")
cp_file <- file.path(fig_dir, paste0("dias_clasificados_UV_PCA_k-means", k, ".csv"))

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
# 1) ENSO phases from Niño3.4 (NOAA operational criterion)
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

# NOAA operational definition: at least 5 consecutive months for El Niño/La Niña
rle_fase <- rle(nino34$fase_raw)
fases_validadas <- ifelse(
  rle_fase$values %in% c("El Niño", "La Niña") & rle_fase$lengths >= 5,
  rle_fase$values,
  "Neutral"
)
nino34$fase_noaa <- inverse.rle(list(lengths = rle_fase$lengths, values = fases_validadas))

# Keep analysis window (EDIT if needed)
nino34 <- nino34[nino34$año %in% 1979:2023, c("año", "mes", "fase_noaa")]

# ----------------------------
# 2) MJO daily: intensity + phase, inactive days set to phase 0
# ----------------------------
mjo_int <- read_excel(mjo_int_file)
mjo_int$fecha <- as.Date(mjo_int$dias_julianos - 2440587.5, origin = "1970-01-01")
mjo_int$dia <- as.integer(format(mjo_int$fecha, "%d"))
mjo_int$mes <- as.integer(format(mjo_int$fecha, "%m"))
mjo_int$año <- as.integer(format(mjo_int$fecha, "%Y"))

mjo_ph <- read_excel(mjo_ph_file)
mjo_ph$fecha <- as.Date(mjo_ph$dias_julianos - 2440587.5, origin = "1970-01-01")
mjo_ph$dia <- as.integer(format(mjo_ph$fecha, "%d"))
mjo_ph$mes <- as.integer(format(mjo_ph$fecha, "%m"))
mjo_ph$año <- as.integer(format(mjo_ph$fecha, "%Y"))

x <- merge(mjo_int, mjo_ph, by = c("año", "mes", "dia"))
x$MJO_phase[x$MJO_intensity < 1] <- 0

# ----------------------------
# 3) Merge CP daily labels
# ----------------------------
cp <- read.csv(cp_file)
stopifnot(all(c("año", "mes", "dia", "CP") %in% names(cp)))

x <- merge(x, cp, by = c("año", "mes", "dia"))

# ----------------------------
# 4) Merge monthly ENSO phase to daily table by (year, month)
# ----------------------------
x <- merge(x, nino34, by = c("año", "mes"))
x$MJO_ENSO <- paste0("MJO ", x$MJO_phase, "+", x$fase_noaa)

MJO_ENSO_comb <- sort(unique(x$MJO_ENSO))

# Seasonal subset pool for Monte Carlo
x_season <- x[x$mes %in% season_months, ]

# ----------------------------
# Helper: Monte Carlo envelopes
# ----------------------------
mc_p_state_given_CP <- function(state_label) {
  # p(state | CP) = N(CP,state) / N(CP)
  do.call(rbind, lapply(1:k, function(cp_i) {
    datos_cp <- x[x$CP == cp_i & x$mes %in% season_months, ]
    N <- nrow(datos_cp)
    if (N == 0) return(data.frame(CP = cp_i, p = NA, p95 = NA, p025 = NA, sign = "no", condicion = state_label))
    
    p_obs <- mean(datos_cp$MJO_ENSO == state_label) * 100
    
    mc_vals <- replicate(n_mc, {
      samp <- sample(x_season$MJO_ENSO, size = N, replace = FALSE)
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

mc_p_CP_given_state <- function(state_label) {
  # p(CP | state) = N(CP,state) / N(state)
  datos_state <- x[x$MJO_ENSO == state_label & x$mes %in% season_months, ]
  N <- nrow(datos_state)
  
  do.call(rbind, lapply(1:k, function(cp_i) {
    if (N == 0) return(data.frame(CP = cp_i, p = NA, p95 = NA, p025 = NA, sign = "no", condicion = state_label))
    
    p_obs <- mean(datos_state$CP == cp_i) * 100
    
    mc_vals <- replicate(n_mc, {
      samp <- sample(x_season$CP, size = N, replace = FALSE)
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

# ----------------------------
# 5A) p(MJO_ENSO | CP) for all combined states (optional diagnostic plot)
# ----------------------------
df_A <- do.call(rbind, lapply(MJO_ENSO_comb, mc_p_state_given_CP))

# Split condicion -> MJO and ENSO
parts <- str_split_fixed(df_A$condicion, "\\+", 2)
df_A$MJO  <- parts[, 1]
df_A$nino <- parts[, 2]
df_A$MJO[df_A$MJO == "MJO 0"] <- "no MJO"

df_A$p <- round(df_A$p, 0)

# (Optional) plot all, faceted by ENSO
g_A <- ggplot(df_A, aes(x = MJO, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_A[!is.na(df_A$p), ], aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("superior" = "red", "inferior" = "cyan", "no" = NA), name = "") +
  geom_text(aes(label = p), size = 3) +
  facet_wrap(~ nino, ncol = 3) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "", y = "CP", title = "p(MJO_ENSO | CP) = N(CP, MJO_ENSO) / N(CP)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 9),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_A)

# ----------------------------
# 5B) p(CP | MJO_ENSO): write per-season CSV + plot subset (non-neutral, active MJO)
# ----------------------------
df_B <- do.call(rbind, lapply(MJO_ENSO_comb, mc_p_CP_given_state))

parts <- str_split_fixed(df_B$condicion, "\\+", 2)
df_B$MJO  <- parts[, 1]
df_B$nino <- parts[, 2]
df_B$MJO[df_B$MJO == "MJO 0"] <- "no MJO"

df_B$p <- round(df_B$p, 0)
df_B$season <- season_name

# Save per-season table for later final figure
write.csv(df_B, file.path(fig_dir, paste0("MJO_ENSO_", season_name, ".csv")), row.names = FALSE)

# Plot only active MJO and non-neutral ENSO (like your original)
df_B_sub <- df_B[df_B$MJO != "no MJO" & df_B$nino != "Neutral", ]

g_B <- ggplot(df_B_sub, aes(x = MJO, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_B_sub[df_B_sub$sign != "no", ], aes(color = sign),
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

# ----------------------------
# 6) FINAL figure: combine DJF/MAM/JJA/SON per-season CSVs
# ----------------------------
df_plot <- rbind(
  read.csv(file.path(fig_dir, "MJO_ENSO_DJF.csv")),
  read.csv(file.path(fig_dir, "MJO_ENSO_MAM.csv")),
  read.csv(file.path(fig_dir, "MJO_ENSO_JJA.csv")),
  read.csv(file.path(fig_dir, "MJO_ENSO_SON.csv"))
)

df_plot$season <- factor(df_plot$season, levels = c("DJF", "MAM", "JJA", "SON"))

# Ensure consistent ordering (optional but makes facets stable)
df_plot$condicion <- factor(df_plot$condicion, levels = c(
  paste0("MJO ", 1:8, "+El Niño"), "MJO 0+El Niño",
  paste0("MJO ", 1:8, "+La Niña"), "MJO 0+La Niña",
  paste0("MJO ", 1:8, "+Neutral"), "MJO 0+Neutral"
))

df_plot$sign[df_plot$sign == "inferior"] <- "lower than the clim"
df_plot$sign[df_plot$sign == "superior"] <- "higher than the clim"

# Plot only active MJO + non-neutral ENSO, and annotate only significant cells
df_plot_sub <- df_plot[df_plot$MJO != "no MJO" & df_plot$nino != "Neutral", ]

g_final <- ggplot(df_plot_sub, aes(x = season, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_plot_sub[df_plot_sub$sign != "no", ], aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(name = "", values = c("higher than the clim" = "red",
                                           "lower than the clim" = "cyan",
                                           "no" = NA)) +
  geom_text(data = df_plot_sub[df_plot_sub$sign != "no", ],
            aes(label = round(p, 1)), size = 3) +
  facet_wrap(~ condicion, nrow = 2, scales = "free_x") +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "", y = "CP", title = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        axis.title.y = element_text(size = 12)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_final)

ggsave(
  filename = file.path(fig_dir, "al_rev_MJO-ENSO_porc_por_CP.png"),
  plot = g_final,
  dpi = 300,
  width = 9,
  height = 5,
  bg = "white"
)


