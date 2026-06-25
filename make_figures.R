# make_figures_clean.R
# Generates the manuscript figures from processed SDP-DICE output data.
# Input workbook: figure_data_learning.xlsx

rm(list = ls())

library(ggplot2)
library(readxl)
library(gridExtra)
library(ggridges)
library(patchwork)
library(dplyr)
library(scales)

# Optional: use BBC style if the bbplot package is available.
if (requireNamespace("bbplot", quietly = TRUE)) {
  library(bbplot)
  base_style <- function() bbc_style()
} else {
  base_style <- function() theme_minimal()
}

DATA_FILE <- "figure_data_learning.xlsx"

# -------------------------------------------------------------------------
# Common plotting settings
# -------------------------------------------------------------------------

cols <- c("With Learning" = "#1380A1", "Without Learning" = "#990000")
ltys <- c("With Learning" = "dashed", "Without Learning" = "solid")

# Common y-axis scales for trajectory figures.
y_cp_limits   <- c(100, 30000)
y_cp_breaks   <- c(100, 500, 1000, 5000, 10000, 30000)
y_miu_limits  <- c(0, 1)
y_miu_breaks  <- seq(0, 1, 0.25)
y_mat_limits  <- c(800, 1200)
y_mat_breaks  <- seq(800, 1200, 100)
y_temp_limits <- c(0, 7)
y_temp_breaks <- seq(0, 7, 1)
y_i_limits    <- c(0, 1500)
y_i_breaks    <- seq(0, 1500, 250)
y_gwp_limits  <- c(0, 6000)
y_gwp_breaks  <- seq(0, 6000, 1000)

read_sheet <- function(sheet) {
  dat <- read_excel(DATA_FILE, sheet = sheet)
  dat$scenario <- factor(dat$scenario, levels = c("With Learning", "Without Learning"))
  dat
}

make_line_panel <- function(dat, title, y_scale) {
  p <- ggplot(dat, aes(x = year, y = val, group = scenario, color = scenario)) +
    base_style() +
    geom_hline(yintercept = 0, size = 1, colour = "#333333") +
    ggtitle(title) +
    geom_line(aes(linetype = scenario), size = 1.2) +
    scale_color_manual(values = cols, drop = FALSE) +
    scale_linetype_manual(values = ltys, drop = FALSE) +
    scale_x_continuous(name = element_blank(), limits = c(2000, 2300)) +
    theme(legend.position = "none", legend.title = element_blank())

  p + y_scale
}

make_trajectory_figure <- function(prefix, outfile) {
  p1 <- make_line_panel(
    read_sheet(paste0(prefix, "_cp")),
    "(a) Carbon Price (2005 USD/tCO2)",
    scale_y_continuous(name = element_blank(), trans = "log2", limits = y_cp_limits,
                       breaks = y_cp_breaks, labels = comma)
  )
  p2 <- make_line_panel(
    read_sheet(paste0(prefix, "_miu")),
    "(b) Carbon Control Rate",
    scale_y_continuous(name = element_blank(), limits = y_miu_limits, breaks = y_miu_breaks,
                       labels = percent_format(accuracy = 1))
  )
  p3 <- make_line_panel(
    read_sheet(paste0(prefix, "_mat")),
    "(c) Carbon Concentration (GtC)",
    scale_y_continuous(name = element_blank(), limits = y_mat_limits, breaks = y_mat_breaks)
  ) + geom_hline(yintercept = 800, size = 1, colour = "#333333")
  p4 <- make_line_panel(
    read_sheet(paste0(prefix, "_temp")),
    "(d) Surface Temperature (°C)",
    scale_y_continuous(name = element_blank(), limits = y_temp_limits, breaks = y_temp_breaks)
  )
  p5 <- make_line_panel(
    read_sheet(paste0(prefix, "_i")),
    "(e) Investment (Trillions USD)",
    scale_y_continuous(name = element_blank(), limits = y_i_limits, breaks = y_i_breaks,
                       labels = comma)
  ) + theme(legend.position = "bottom", legend.text = element_text(size = 13.5))
  p6 <- make_line_panel(
    read_sheet(paste0(prefix, "_gwp")),
    "(f) Gross World Output (Trillions USD)",
    scale_y_continuous(name = element_blank(), limits = y_gwp_limits, breaks = y_gwp_breaks,
                       labels = comma)
  ) + theme(legend.position = "bottom", legend.text = element_text(size = 13.5))

  fig <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2)
  ggsave(outfile, plot = fig, width = 16, height = 12, dpi = 300)
}

# -------------------------------------------------------------------------
# Figures 2-5 and Appendix Figure A10: trajectory figures
# -------------------------------------------------------------------------

make_trajectory_figure("8c",  "ECS8C.png")
make_trajectory_figure("6c",  "ECS6C.png")
make_trajectory_figure("1c",  "ECS1C.png")
make_trajectory_figure("32c", "ECS3C.png")
make_trajectory_figure("sens_shock42X", "ECS_sens_042.png")

# -------------------------------------------------------------------------
# Figure 1: evolution of ECS beliefs
# -------------------------------------------------------------------------

lambda <- 1.2

generate_ecs_data <- function(mu, sigma, year, x_range = seq(0, 10, 0.05)) {
  density_vals <- 1 / sqrt(2 * pi * sigma^2) * lambda / x_range^2 *
    exp(-0.5 * (1 - lambda / x_range - mu)^2 / sigma^2)
  data.frame(x = x_range, density = density_vals, year = as.factor(year))
}

create_ridgeline_plot <- function(params, true_ecs, plot_title) {
  ecs_data <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    generate_ecs_data(params$mu[i], params$sigma[i], params$year[i])
  }))

  year_colors <- c(
    "2020" = "#3B518A", "2050" = "#5F5399", "2100" = "#7B5190",
    "2150" = "#A04F7E", "2200" = "#C5515F", "2250" = "#E05745"
  )
  ecs_data$year <- factor(ecs_data$year, levels = sort(unique(ecs_data$year)))

  ggplot(ecs_data, aes(x = x, y = year, height = density, fill = year, color = year)) +
    geom_density_ridges(stat = "identity", scale = 3, alpha = 0.7, size = 0.75,
                        rel_min_height = 0.01, show.legend = FALSE) +
    geom_segment(aes(x = 0, xend = 10, y = year, yend = year, color = year), size = 0.5) +
    geom_vline(xintercept = true_ecs, linetype = "dashed", color = "#D62728", size = 0.75) +
    annotate("text", x = true_ecs + 0.15, y = 1.4, label = "True ECS Value",
             color = "#D62728", hjust = 0, fontface = "bold", size = 3.5) +
    scale_fill_manual(values = year_colors) +
    scale_color_manual(values = year_colors) +
    scale_x_continuous(name = "Equilibrium Climate Sensitivity (°C)",
                       breaks = seq(0, 10, 1), limits = c(0, 10), expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    geom_hline(yintercept = 0.4, size = 1, color = "#333333") +
    theme_minimal() +
    theme(panel.grid.major.x = element_line(color = "#EEEEEE"),
          panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 18, face = "bold", hjust = 0),
          axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 15)),
          axis.title.y = element_blank(), axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14, face = "bold", color = "#333333"),
          legend.position = "none", plot.margin = margin(t = 20, r = 30, b = 20, l = 20)) +
    labs(title = plot_title)
}

params_8C <- data.frame(
  year = c(2020, 2050, 2100, 2150, 2200, 2250),
  mu = c(0.6, 0.749401, 0.827556, 0.839808, 0.844876, 0.846021),
  sigma = c(0.1, sqrt(0.004013), sqrt(0.000913), sqrt(0.000359), sqrt(0.000196), sqrt(0.000127))
)
params_6C <- data.frame(
  year = c(2020, 2050, 2100, 2150, 2200, 2250),
  mu = c(0.6, 0.71243, 0.7776371, 0.79101543778, 0.7938874847, 0.79586583222),
  sigma = c(0.1, sqrt(0.0042073152), sqrt(0.00105607923), sqrt(0.000450877168),
            sqrt(0.0002585068), sqrt(0.00017615477))
)
params_3C <- data.frame(
  year = c(2020, 2050, 2100, 2150, 2200, 2250),
  mu = c(0.6, 0.609839916, 0.61744706, 0.620655277, 0.621400555, 0.622026583),
  sigma = c(0.1, sqrt(0.004817705), sqrt(0.001547363), sqrt(0.000775851),
            sqrt(0.000500941), sqrt(0.000367819))
)
params_1C <- data.frame(
  year = c(2020, 2050, 2100, 2150, 2200, 2250),
  mu = c(0.6, 0.366368, 0.138276, 0.024355, -0.032081, -0.065933),
  sigma = c(0.1, sqrt(0.007124), sqrt(0.004238), sqrt(0.002823),
            sqrt(0.002124), sqrt(0.001713))
)

plot_3C <- create_ridgeline_plot(params_3C, 3.2, "(a) True ECS = 3.2°C")
plot_8C <- create_ridgeline_plot(params_8C, 8, "(b) True ECS = 8°C")
plot_6C <- create_ridgeline_plot(params_6C, 6, "(c) True ECS = 6°C")
plot_1C <- create_ridgeline_plot(params_1C, 1, "(d) True ECS = 1°C")
fig1 <- (plot_3C + plot_8C) / (plot_6C + plot_1C)
ggsave("Learning_ECS_Dist.png", fig1, width = 16, height = 12, dpi = 600)

# Appendix Figure A9: ECS belief sensitivity to temperature variability
params_042X <- data.frame(
  year = c(2020, 2050, 2100, 2150, 2200, 2250),
  mu = c(0.6, 0.747681996970726, 0.825108918898649, 0.839644096452692,
         0.843388348221091, 0.845751680583133),
  sigma = c(0.1, 0.0632514045034756, 0.0300361439501199, 0.0188790919862679,
            0.0139428684918792, 0.0112345140618276)
)
plot_8C_default <- create_ridgeline_plot(params_8C, 8, "(a) True ECS = 8°C with the Default Temperature Shock")
plot_8C_042X <- create_ridgeline_plot(params_042X, 8, "(b) True ECS = 8°C with a Larger Temperature Shock")
figA9 <- plot_8C_default + plot_8C_042X
ggsave("Learning_ECS_Dist_SENS8C.png", figA9, width = 16, height = 6, dpi = 600)

# -------------------------------------------------------------------------
# Figure 7: initial ECS prior distribution
# -------------------------------------------------------------------------

mu <- 0.6
sigma <- 0.1

ecs_pdf <- function(x, lambda = 1.2, mu = 0.6, sigma = 0.1) {
  ifelse(x > 0, (1 / sqrt(2 * pi * sigma^2)) * lambda / x^2 *
           exp(-0.5 * (1 - lambda / x - mu)^2 / sigma^2), 0)
}

fig7 <- ggplot(data.frame(x = c(0, 10)), aes(x = x)) +
  theme_minimal() +
  stat_function(fun = ecs_pdf, geom = "line", color = "#3B518A", size = 1.2) +
  stat_function(fun = ecs_pdf, geom = "area", fill = "#3B518A", alpha = 0.15) +
  geom_ribbon(data = data.frame(x = c(2, 5)), aes(x = x, ymin = 0, ymax = 0.6),
              fill = "grey30", alpha = 0.15, inherit.aes = FALSE) +
  geom_vline(xintercept = c(2, 5), linetype = "dashed", color = "#D62728", size = 0.75) +
  geom_point(data = data.frame(x = 3.2, y = 0), aes(x = x, y = y),
             color = "#D62728", size = 4, inherit.aes = FALSE) +
  geom_point(data = data.frame(x = c(1, 6, 8), y = 0), aes(x = x, y = y),
             color = "#7A2180", size = 4, inherit.aes = FALSE) +
  annotate("text", x = 3.9, y = 0.59, label = "IPCC’s Very Likely Range",
           color = "#555555", fontface = "bold", size = 4) +
  scale_x_continuous(name = "Equilibrium Climate Sensitivity (°C)", breaks = seq(0, 10, 1),
                     limits = c(0, 10), expand = c(0, 0)) +
  scale_y_continuous(name = "Density", limits = c(0, 0.6), expand = c(0, 0.01)) +
  theme(panel.grid.major.x = element_line(color = "#EEEEEE"),
        panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(color = "#EEEEEE"),
        panel.grid.minor.y = element_blank(), legend.position = "none",
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 15)))
scenario_labels <- data.frame(
  x = c(3.2, 1, 6, 8), y = 0.05,
  label = c("Best Estimate (3.2°C)", "Low Scenario (1°C)", "High Scenario (6°C)", "Extreme Scenario (8°C)"),
  color = c("#D62728", "#7A2180", "#7A2180", "#7A2180")
)
fig7 <- fig7 +
  geom_text(data = scenario_labels, aes(x = x, y = y, label = label, color = I(color)),
            size = 3, vjust = -0.5, fontface = "bold", inherit.aes = FALSE) +
  geom_segment(data = scenario_labels, aes(x = x, xend = x, y = 0, yend = 0.04, color = I(color)),
               linetype = "dotted", size = 0.5, inherit.aes = FALSE)
ggsave("ECS_Initial_Distribution.png", fig7, width = 11, height = 7, dpi = 300)

# -------------------------------------------------------------------------
# Figure 8: damage functions
# -------------------------------------------------------------------------

dam_fn <- read_excel(DATA_FILE, sheet = "dam_fn")
dam_fn$scenario <- ifelse(dam_fn$scenario == "Modified Damage Function (Howard and Sterner, 2018)",
                          "Howard and Sterner (2017)", dam_fn$scenario)
fig8 <- ggplot(dam_fn, aes(x = temp, y = val, color = scenario, linetype = scenario)) +
  geom_hline(yintercept = seq(0, 1, by = 0.25), color = "#EEEEEE", size = 0.5) +
  annotate("rect", xmin = 0, xmax = 2, ymin = 0, ymax = 1, fill = "#7fbf7b", alpha = 0.1) +
  annotate("rect", xmin = 2, xmax = 4, ymin = 0, ymax = 1, fill = "#fcae60", alpha = 0.1) +
  annotate("rect", xmin = 4, xmax = 6, ymin = 0, ymax = 1, fill = "#f46d43", alpha = 0.1) +
  annotate("rect", xmin = 6, xmax = 11, ymin = 0, ymax = 1, fill = "#d73027", alpha = 0.1) +
  geom_line(size = 1.2) +
  scale_x_continuous(name = "Global Mean Temperature Increase (°C)", breaks = seq(0, 10, by = 5),
                     limits = c(0, 11), expand = c(0, 0)) +
  scale_y_continuous(name = "Damage (% of Gross World Product)", labels = percent_format(accuracy = 1),
                     limits = c(0, 1), expand = c(0, 0)) +
  scale_color_manual(values = c("#3B518A", "#D62728")) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "#EEEEEE"), panel.grid.minor.y = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 15)))
ggsave("Damage_Function.png", fig8, width = 10, height = 7, dpi = 300)

# -------------------------------------------------------------------------
# Figure 6: Value of learning across ECS cases
# -------------------------------------------------------------------------

value_data <- data.frame(
  ecs = c(1, 3.2, 6, 8),
  value_trillion = c(23, 2, 15, 37),
  scenario = c("True ECS = 1°C", "True ECS = 3.2°C", "True ECS = 6°C", "True ECS = 8°C"),
  color = c("#7A2180", "#D62728", "#7A2180", "#7A2180")
)
value_data$density <- sapply(value_data$ecs, ecs_pdf)

x_seq <- seq(0.5, 10, by = 0.01)
dist_data <- data.frame(x = x_seq, density = sapply(x_seq, ecs_pdf))
y_scale <- 0.012
fig6 <- ggplot() +
  geom_area(data = dist_data, aes(x = x, y = density), fill = "#E8E8E8",
            alpha = 0.5, color = "#CCCCCC", size = 0.5) +
  geom_ribbon(data = data.frame(x = c(2, 5)), aes(x = x, ymin = 0, ymax = 0.6),
              fill = "#DDDDDD", alpha = 0.3) +
  geom_vline(xintercept = c(2, 5), linetype = "dashed", color = "#999999", size = 0.5) +
  geom_col(data = value_data, aes(x = ecs, y = value_trillion * y_scale),
           width = 0.25, fill = "#2E86AB", alpha = 0.8, color = "#1F5F7A", size = 0.5) +
  geom_text(data = value_data, aes(x = ecs, y = value_trillion * y_scale + 0.03,
                                   label = paste0("$", value_trillion, "T")),
            size = 3.8, fontface = "bold", color = "#1F5F7A") +
  geom_text(data = value_data, aes(x = ecs, y = -0.05, label = scenario, color = I(color)),
            size = 3.2, fontface = "bold", vjust = 1) +
  geom_point(data = value_data, aes(x = ecs, y = density), size = 4, color = "#FF6B35", alpha = 0.8) +
  annotate("text", x = 3.5, y = 0.55, label = "IPCC Very Likely Range\n(2°C - 5°C)",
           color = "#666666", size = 3.5, fontface = "italic") +
  scale_x_continuous(name = "Equilibrium Climate Sensitivity (°C)",
                     breaks = c(1, 2, 3.2, 4, 5, 6, 7, 8, 9), limits = c(0.5, 9.5), expand = c(0, 0)) +
  scale_y_continuous(name = "Probability Density", limits = c(-0.12, 0.6), expand = c(0, 0),
                     sec.axis = sec_axis(~ ./y_scale, name = "Value of Learning ($ Trillions)")) +
  theme_minimal() +
  theme(panel.grid.major.x = element_line(color = "#F5F5F5", size = 0.5),
        panel.grid.minor.x = element_blank(), panel.grid.major.y = element_line(color = "#F5F5F5", size = 0.5),
        panel.grid.minor.y = element_blank(), axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 20)),
        axis.title.y.left = element_text(size = 12, face = "bold", margin = margin(r = 15), color = "#999999"),
        axis.title.y.right = element_text(size = 12, face = "bold", margin = margin(l = 15), color = "#2E86AB"))
ggsave("value_of_learning_plot.png", fig6, width = 16, height = 12, dpi = 600, bg = "white")

# -------------------------------------------------------------------------
# Appendix Figure A11: alternative ECS priors for EVL sensitivity
# -------------------------------------------------------------------------

roe_baker_pdf <- function(x, sigma_rb = 0.1) ecs_pdf(x, lambda = 1.2, mu = 0.6, sigma = sigma_rb)

# AR6-consistent lognormal prior: 66% likely interval matched to 2.5-4.0°C.
q_lo <- 2.5; p_lo <- 0.17
q_hi <- 4.0; p_hi <- 0.83
sd_ln_thin <- (log(q_hi) - log(q_lo)) / (qnorm(p_hi) - qnorm(p_lo))
mu_ln_thin <- log(q_lo) - sd_ln_thin * qnorm(p_lo)

prior_data <- data.frame(
  ecs = rep(x_seq, 3),
  density = c(roe_baker_pdf(x_seq, 0.1),
              roe_baker_pdf(x_seq, 0.12),
              dlnorm(x_seq, meanlog = mu_ln_thin, sdlog = sd_ln_thin)),
  prior = rep(c("Baseline (σ=0.1)", "Fatter-tail (σ=0.12)", "Thin-tail (AR6)"), each = length(x_seq))
)
figA11 <- ggplot(prior_data, aes(x = ecs, y = density, color = prior)) +
  geom_line(size = 1.1) +
  geom_vline(xintercept = c(2, 5), linetype = "dashed", color = "#999999") +
  annotate("text", x = 3.5, y = 0.62, label = "IPCC Very Likely Range (2°C-5°C)",
           size = 3.5, color = "#555555") +
  scale_x_continuous(name = "Equilibrium Climate Sensitivity (°C)", limits = c(0.5, 10)) +
  scale_y_continuous(name = "Density", limits = c(0, 0.7)) +
  scale_color_manual(values = c("#3B518A", "#D62728", "#7B5190")) +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank())
ggsave("Alternative_ECS_Priors.png", figA11, width = 10, height = 6, dpi = 300)

