library(posterior)
draws <- as_draws_rvars(fit)

# --- phi --- #
phi_rv <- draws$phi
phi_matrix <- posterior::as_draws_matrix(phi_rv)

phi_mean_site_year <- matrix(
  colMeans(phi_matrix),
  nrow = dim(phi_rv)[1],
  byrow = FALSE
)

phi_summary <- apply(phi_mean_site_year, 2, function(x) {
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
})

phi_df <- data.frame(
  year = 2019:2025,
  mean = phi_summary["mean", ],
  se   = phi_summary["se", ],
  param = "Persistence"
)


# --- gamma --- #
gamma_rv <- draws$gamma
gamma_matrix <- posterior::as_draws_matrix(gamma_rv)

gamma_mean_site_year <- matrix(
  colMeans(gamma_matrix),
  nrow = dim(gamma_rv)[1],
  byrow = FALSE
)

gamma_summary <- apply(gamma_mean_site_year, 2, function(x) {
  c(mean = mean(x), se = sd(x) / sqrt(length(x)))
})

gamma_df <- data.frame(
  year = 2019:2025,
  mean = gamma_summary["mean", ],
  se   = gamma_summary["se", ],
  param = "Colonization"
)

pg_df <- rbind(phi_df, gamma_df)

library(ggplot2)

ggplot(pg_df, aes(x = year, y = mean, color = param, fill = param)) +
  geom_ribbon(
    aes(ymin = mean - se, ymax = mean + se),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(size = 1.3) +
  geom_point(size = 3.5) +
  theme_classic(base_size = 18) +
  labs(
    x = "Year",
    y = "Probability",
    color = "Parameter",
    fill = "Parameter"
  ) +
  
  # ----- Axis increments -----
scale_x_continuous(
  breaks = seq(min(pg_df$year), max(pg_df$year), by = 1)
) +
  scale_y_continuous(
    limits = c(0,1),
    breaks = seq(0, 1, by = 0.1),
    expand = c(0,0.01)
  ) +
  
  # ----- Colors -----
scale_color_manual(values = c(
  Persistence = "steelblue4",
  Colonization = "darkgreen"
)) +
  scale_fill_manual(values = c(
    Persistence = "steelblue",
    Colonization = "darkgreen"
  ))
