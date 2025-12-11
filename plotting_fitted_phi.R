plot_fitted_phi <- function(fit, stan_data, ci = 0.95) {
  library(posterior)
  library(ggplot2)
  
  draws <- as_draws_df(fit)
  
  # Extract posterior draws
  phi0_samples   <- draws$phi_0
  beta_phi_draws <- as.matrix(draws[, grep("^beta_phi", names(draws))])
  
  n_iter <- nrow(beta_phi_draws)
  nsite  <- stan_data$nsite
  nyear  <- stan_data$nyear
  Xphi   <- stan_data$Xphi  # [site, year, Kphi]
  
  # output matrix: φ(iter, year)
  phi_year_iter <- matrix(NA, n_iter, nyear)
  
  for (iter in 1:n_iter) {
    beta_i <- beta_phi_draws[iter, ]
    phi0_i <- phi0_samples[iter]
    
    # compute φ_hat[site, year]
    phi_hat <- array(NA, c(nsite, nyear))
    
    for (s in 1:nsite) {
      for (t in 1:nyear) {
        lp <- phi0_i + sum(Xphi[s,t,] * beta_i)
        phi_hat[s, t] <- plogis(lp)
      }
    }
    
    # average φ across sites for this iteration → φ(t)
    phi_year_iter[iter, ] <- colMeans(phi_hat)
  }
  
  # posterior summaries
  phi_mean  <- colMeans(phi_year_iter)
  phi_lower <- apply(phi_year_iter, 2, quantile, probs = (1 - ci)/2)
  phi_upper <- apply(phi_year_iter, 2, quantile, probs = 1 - (1 - ci)/2)
  
  df <- data.frame(
    year  = 1:nyear,
    mean  = phi_mean,
    lower = phi_lower,
    upper = phi_upper
  )
  
  # plot
  ggplot(df, aes(x = year)) +
    geom_line(aes(y = mean), color = "blue", size = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.25) +
    ylab("Survival (phi)") +
    xlab("Year") +
    theme_classic() +
    ggtitle("Fitted φ(t)")
}

plot_fitted_phi(fit, stan_data)
