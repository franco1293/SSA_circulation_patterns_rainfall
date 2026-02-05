# =============================================================================
# MJO phase vs Circulation Patterns (CPs)
# =============================================================================

# ----------------------------
# Libraries
# ----------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ----------------------------
# User parameters (EDIT)
# ----------------------------
k <- 7
season_name <- "MAM"  # "ONDJFMA", "DJF", "MAM", "JJA", "SON"

base_dir <- "C:/Users/Franco/Documents/posdoctorado/WTs_1st"
fig_dir  <- file.path(base_dir, "Figures")

mjo_phase_file     <- file.path(base_dir, "MJO_phase.xlsx")
mjo_intensity_file <- file.path(base_dir, "MJO_intensity.xlsx")

cp_file <- file.path(
  fig_dir,
  paste0("dias_clasificados_UV_PCA_k-means", k, ".csv")
)

n_mc <- 10000
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
# 1) Read MJO phase + intensity and build daily table
# ----------------------------
x <- read_excel(mjo_phase_file)
z <- read_excel(mjo_intensity_file)

stopifnot(all(c("dias_julianos") %in% names(x)))
stopifnot(all(c("dias_julianos") %in% names(z)))

# Convert "dias_julianos" -> Date (keep your original conversion)
x <- x %>%
  mutate(
    fecha = as.Date(dias_julianos - 2440587.5, origin = "1970-01-01"),
    dia   = as.integer(format(fecha, "%d")),
    mes   = as.integer(format(fecha, "%m")),
    año   = as.integer(format(fecha, "%Y"))
  )

z <- z %>%
  mutate(
    fecha = as.Date(dias_julianos - 2440587.5, origin = "1970-01-01"),
    dia   = as.integer(format(fecha, "%d")),
    mes   = as.integer(format(fecha, "%m")),
    año   = as.integer(format(fecha, "%Y"))
  )

# Merge MJO phase + intensity by date parts
x <- merge(x, z, by = c("año", "mes", "dia"))

# Inactive MJO definition
stopifnot(all(c("MJO_phase", "MJO_intensity") %in% names(x)))
x$MJO_phase[x$MJO_intensity < 1] <- 0

# ----------------------------
# 2) Merge with CP classifications
# ----------------------------
y <- read.csv(cp_file)
stopifnot(all(c("año", "mes", "dia", "CP") %in% names(y)))

x <- merge(x, y, by = c("año", "mes", "dia"))

# Seasonal subset
x_season <- x[x$mes %in% season_months, ]
MJO_levels <- sort(unique(x$MJO_phase))

# ----------------------------
# 3A) p(MJO=i | CP=cp): % of CP days that have MJO phase i
# ----------------------------
mc_p_MJO_given_CP <- function(mjo_i, x_all, x_season, season_months, n_mc = 10000, k = 7) {
  res <- lapply(1:k, function(cp_i) {
    datos_cp <- x_all[x_all$CP == cp_i & x_all$mes %in% season_months, ]
    N <- nrow(datos_cp)
    if (N == 0) {
      return(data.frame(CP = cp_i, p = NA, p97.5 = NA, p2.5 = NA, sign = NA, MJO_phase = mjo_i))
    }
    
    p_obs <- mean(datos_cp$MJO_phase == mjo_i) * 100
    
    # Monte Carlo: sample MJO phases from seasonal pool (without replacement)
    mc_vals <- replicate(n_mc, {
      samp <- sample(x_season$MJO_phase, size = N, replace = FALSE)
      mean(samp == mjo_i) * 100
    })
    
    p97.5 <- as.numeric(quantile(mc_vals, 0.975, na.rm = TRUE))
    p2.5  <- as.numeric(quantile(mc_vals, 0.025, na.rm = TRUE))
    
    sign <- "no"
    if (!is.na(p_obs) && !is.na(p97.5) && p_obs > p97.5) sign <- "superior"
    if (!is.na(p_obs) && !is.na(p2.5)  && p_obs < p2.5)  sign <- "inferior"
    
    data.frame(CP = cp_i, p = round(p_obs, 1), p97.5 = round(p97.5, 1), p2.5 = round(p2.5, 1),
               sign = sign, MJO_phase = mjo_i)
  })
  
  do.call(rbind, res)
}

df_A <- do.call(rbind, lapply(MJO_levels, mc_p_MJO_given_CP,
                              x_all = x, x_season = x_season,
                              season_months = season_months, n_mc = n_mc, k = k))

df_A$MJO_phase <- ifelse(df_A$MJO_phase == 0, "no MJO", paste("MJO", df_A$MJO_phase))
df_A$p <- round(df_A$p, 0)

g_A <- ggplot(df_A, aes(x = MJO_phase, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_A[!is.na(df_A$p), ], aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("superior" = "red", "inferior" = "cyan", "no" = NA), name = "") +
  geom_text(aes(label = p), size = 4) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "MJO phase", y = "CP",
       title = "p(MJO=i | CP=cp) = N(CP=cp and MJO=i) / N(CP=cp)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_A)

# Optional save (avoid overwriting the final panel name)
ggsave(
  filename = file.path(fig_dir, paste0("MJO_given_CP_", season_name, ".png")),
  plot = g_A, dpi = 300, width = 4, height = 3, bg = "white"
)

# ----------------------------
# 3B) p(CP=cp | MJO=i): % of MJO phase i days that have CP=cp
# ----------------------------
mc_p_CP_given_MJO <- function(mjo_i, x_all, x_season, season_months, n_mc = 10000, k = 7) {
  datos_mjo <- x_all[x_all$MJO_phase == mjo_i & x_all$mes %in% season_months, ]
  N <- nrow(datos_mjo)
  
  res <- lapply(1:k, function(cp_i) {
    if (N == 0) {
      return(data.frame(CP = cp_i, p = NA, p97.5 = NA, p2.5 = NA, sign = NA, MJO_phase = mjo_i))
    }
    
    p_obs <- mean(datos_mjo$CP == cp_i) * 100
    
    # Monte Carlo: sample CP labels from seasonal pool (without replacement)
    mc_vals <- replicate(n_mc, {
      samp <- sample(x_season$CP, size = N, replace = FALSE)
      mean(samp == cp_i) * 100
    })
    
    p97.5 <- as.numeric(quantile(mc_vals, 0.975, na.rm = TRUE))
    p2.5  <- as.numeric(quantile(mc_vals, 0.025, na.rm = TRUE))
    
    sign <- "no"
    if (!is.na(p_obs) && !is.na(p97.5) && p_obs > p97.5) sign <- "superior"
    if (!is.na(p_obs) && !is.na(p2.5)  && p_obs < p2.5)  sign <- "inferior"
    
    data.frame(CP = cp_i, p = round(p_obs, 1), p97.5 = round(p97.5, 1), p2.5 = round(p2.5, 1),
               sign = sign, MJO_phase = mjo_i)
  })
  
  do.call(rbind, res)
}

df_B <- do.call(rbind, lapply(MJO_levels, mc_p_CP_given_MJO,
                              x_all = x, x_season = x_season,
                              season_months = season_months, n_mc = n_mc, k = k))

df_B$MJO_phase <- ifelse(df_B$MJO_phase == 0, "no MJO", paste("MJO", df_B$MJO_phase))
df_B$p <- round(df_B$p, 0)
df_B$season <- season_name

# Save CSV per season (this is what you later rbind)
write.csv(df_B, file.path(fig_dir, paste0("MJO_", season_name, ".csv")), row.names = FALSE)

g_B <- ggplot(df_B, aes(x = MJO_phase, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_B[!is.na(df_B$p), ], aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(values = c("superior" = "red", "inferior" = "cyan", "no" = NA), name = "") +
  geom_text(aes(label = p), size = 4) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "MJO phase", y = "CP",
       title = paste0(season_name, "  —  p(CP=cp | MJO=i)")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 10),
        legend.position = "none") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_B)

ggsave(
  filename = file.path(fig_dir, paste0("CP_given_MJO_", season_name, ".png")),
  plot = g_B, dpi = 300, width = 4, height = 3, bg = "white"
)

# ----------------------------
# 4) Final multi-season figure (requires MJO_DJF/MAM/JJA/SON csv already saved)
# ----------------------------
df_plot <- rbind(
  read.csv(file.path(fig_dir, "MJO_DJF.csv")),
  read.csv(file.path(fig_dir, "MJO_MAM.csv")),
  read.csv(file.path(fig_dir, "MJO_JJA.csv")),
  read.csv(file.path(fig_dir, "MJO_SON.csv"))
)

# Keep only active MJO phases (exclude phase 0)
df_plot <- df_plot[df_plot$MJO_phase != "no MJO", ]
df_plot$season <- factor(df_plot$season, levels = c("DJF", "MAM", "JJA", "SON"))

# Rename sign labels for legend clarity
df_plot$sign[df_plot$sign == "inferior"] <- "lower than the clim"
df_plot$sign[df_plot$sign == "superior"] <- "higher than the clim"

g_final <- ggplot(df_plot, aes(x = season, y = factor(CP), fill = p)) +
  geom_tile(color = "black") +
  geom_tile(data = df_plot[df_plot$sign != "no", ], aes(color = sign),
            size = 1, alpha = .5, width = 0.9, height = 0.9) +
  scale_color_manual(
    name = "",
    values = c("higher than the clim" = "red", "lower than the clim" = "cyan", "no" = NA)
  ) +
  facet_wrap(~ MJO_phase, nrow = 1) +
  geom_text(data = df_plot[df_plot$sign != "no", ], aes(label = round(p, 1)), size = 4) +
  scale_fill_gradient2(low = "red1", high = "red4", na.value = "white",
                       name = "p (%)", limits = c(0, 100)) +
  labs(x = "", y = "CP", title = "") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0))

plot(g_final)

ggsave(
  filename = file.path(fig_dir, "MJO_CP_association_allSeasons.png"),
  plot = g_final, dpi = 300, width = 9, height = 3.5, bg = "white"
)
