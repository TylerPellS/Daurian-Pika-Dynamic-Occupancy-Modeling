library(posterior)
post <- rstan::extract(fit)

# phi
phi0_vec <- post$phi_0
beta_phi_mat <- post$beta_phi  # draws × Kphi

k <- 1# index of the factor covariate in beta_phi_mat

# Factor levels
levels <- c(0, 1)  # Unburned = 0, Burned = 1

# Number of posterior draws
nd <- length(phi0_vec)

# Initialize matrix to store probabilities
phi_factor <- matrix(NA, nrow = nd, ncol = length(levels))

# Compute linear predictor for each factor level
for(i in seq_along(levels)){
  eta <- phi0_vec + beta_phi_mat[, k] * levels[i]
  phi_factor[, i] <- plogis(eta)  # transform to probability
}

# Summarize for plotting
phi_factor_df <- data.frame(
  level = c("Normal", "Flood"),
  mean = apply(phi_factor, 2, mean),
  low95 = apply(phi_factor, 2, quantile, 0.025),
  high95 = apply(phi_factor, 2, quantile, 0.975)
)

library(ggplot2)

flood <- ggplot(phi_factor_df, aes(x = level, y = mean)) +
  
  # Bars
  geom_col(fill = "steelblue", alpha = 0.8, width = 0.6) +
  
  # Error bars
  geom_errorbar(aes(ymin = low95, ymax = high95),
                width = 0.15,
                color = "black",
                size = 1.0) +
  
  # Optional: add points on top of bars
  geom_point(aes(y = mean),
             size = 4,
             color = "steelblue4") +
  
  # Labels
  labs(
    x = "Flood Status",
    y = "Persistence Probability (φ)",
    title = "Effect of Flood on Persistence"
  ) +
  
  # Axes limits and breaks
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), expand = c(0,0.01)) +
  
  # Theme
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text  = element_text(size = 16),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray40"),
    plot.margin = margin(15, 20, 15, 20)
  )
library(patchwork)
# Combine plots side-by-side and remove y-axis elements from the second plot
combined_plot <- flood + one +
  plot_layout(guides = "collect") & # Collects legends if present
  theme(axis.title.y.right = element_blank(), # Remove title for right plot's y-axis
        axis.text.y.right = element_blank(),  # Remove text for right plot's y-axis
        axis.ticks.y.right = element_blank()) # Remove ticks for right plot's y-axis

# Display the combined plot
combined_plot
