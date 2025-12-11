plot_fitted_gamma <- function(fit, stan_data, ci = 0.95) {
  library(posterior)
  library(ggplot2)
  
  draws <- as_draws_df(fit)
  
  # Extract posterior draws
  gamma0_samples   <- draws$gamma_0
  beta_gam_draws   <- as.matrix(draws[, grep("^beta_gam", names(draws))])
  
  n_iter <- nrow(beta_gam_draws)
  nsite  <- stan_data$nsite
  nyear  <- stan_data$nyear
  Xgam   <- stan_data$Xgam  # [site, year, Kgam]
  
  # output matrix: γ(iter, year)
  gamma_year_iter <- matrix(NA, n_iter, nyear)
  
  for (iter in 1:n_iter) {
    beta_i  <- beta_gam_draws[iter, ]
    gam0_i  <- gamma0_samples[iter]
    
    # compute γ_hat[site, year]
    gamma_hat <- array(NA, c(nsite, nyear))
    
    for (s in 1:nsite) {
      for (t in 1:nyear) {
        lp <- gam0_i + sum(Xgam[s,t,] * beta_i)
        gamma_hat[s, t] <- plogis(lp)
      }
    }
    
    # average γ across sites for this iteration → γ(t)
    gamma_year_iter[iter, ] <- colMeans(gamma_hat)
  }
  
  # posterior summaries
  gamma_mean  <- colMeans(gamma_year_iter)
  gamma_lower <- apply(gamma_year_iter, 2, quantile, probs = (1 - ci)/2)
  gamma_upper <- apply(gamma_year_iter, 2, quantile, probs = 1 - (1 - ci)/2)
  
  df <- data.frame(
    year  = 1:nyear,
    mean  = gamma_mean,
    lower = gamma_lower,
    upper = gamma_upper
  )
  
  # plot
  ggplot(df, aes(x = year)) +
    geom_line(aes(y = mean), color = "darkgreen", size = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "darkgreen", alpha = 0.25) +
    ylab("Colonization (gamma)") +
    xlab("Year") +
    theme_classic() +
    ggtitle("Fitted γ(t)")
}
plot_fitted_gamma(fit, stan_data)
