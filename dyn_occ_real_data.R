



library(dplyr)
library(tidyr)
library(purrr)
library(rstan)
library(parallel)


# read in data
pika <- read.csv("dpdata.csv")
head(pika)
pika$site <- rep(c(1:87))
pika <- pika[,c(14,1:12)] # looks good


pika$yearsc <- pika$Year - 2022

# Global Data Settings
nsite  <- 87
nyear  <- 7
nsurv  <- 4
survey_years <- c(2019, 2022, 2023, 2025)

#  Xpsi ---- matrix of site x covariates x year for initial occupancy
Xpsi <- pika %>%
  filter(Year == 2019) %>%   # Xpsi is only for initial occupancy
  arrange(site) %>%
  select(sGrassHt2019, sLichen2019, sReliefPlot) %>%
  as.matrix()

Kpsi <- ncol(Xpsi)

# array[nsite] of nyear × Kphi matrices
phi_covs <- c("PreFlood", "sOIS", "yearsc")

Xphi_list <- pika %>%
  arrange(site, Year) %>%
  group_by(site) %>%
  summarise(mat = list(as.matrix(select(cur_data(), all_of(phi_covs))))) %>%
  pull(mat)

# Each element must be nyear × Kphi
stopifnot(length(Xphi_list) == nsite)
stopifnot(all(sapply(Xphi_list, nrow) == nyear))

Kphi <- length(phi_covs)

#Xgam array[nsite] of nyear × Kgam matrices

gam_covs <- c("sRelief200L", "yearsc")

Xgam_list <- pika %>%
  arrange(site, Year) %>%
  group_by(site) %>%
  summarise(mat = list(as.matrix(select(cur_data(), all_of(gam_covs))))) %>%
  pull(mat)

Kgam <- 2

#Y matrix nsite × nsurv
Y <- pika %>%
  filter(Year %in% survey_years) %>%
  arrange(site, Year) %>%
  select(site, Year, Occu) %>%
  pivot_wider(names_from = Year, values_from = Occu) %>%
  arrange(site) %>%
  select(all_of(as.character(survey_years))) %>%
  as.matrix()

# intervals
interval <- c(3, 1, 2)

# --- FINAL STAN DATA ----
stan_data <- list(
  nsite = nsite,
  nyear = nyear,
  nsurv = nsurv,
  interval = interval,
  Kpsi = Kpsi,
  Kphi = Kphi,
  Kgam = Kgam,
  Xpsi = Xpsi,
  Xphi = Xphi_list,
  Xgam = Xgam_list,
  Y = Y
)

modd <- stan_model("dyn_occ_v3.stan")

fit <- sampling(modd, data = stan_data, warmup = 500, iter=1000, chains = 4)





# ---- Diagnostics ---- #
draws = as_draws_df(fit)

bayesplot::mcmc_trace(draws, "psi1_0")
bayesplot::mcmc_trace(draws,"phi_0")
bayesplot::mcmc_trace(draws,"gamma_0")
bayesplot::mcmc_trace(draws,"beta_psi[2]") #beta coefficient of each cov on psi
bayesplot::mcmc_trace(draws,"beta_phi[2]") #beta coefficient of each cov on phi
bayesplot::mcmc_trace(draws,"beta_gam[2]") #beta coefficient of each cov on gamma
mean(draws$`beta_psi[1]`) #0.8 grass height 2019
mean(draws$`beta_psi[2]`) #0.84 lichen 2019
mean(draws$`beta_psi[3]`) #1.03 relief plot level
mean(draws$`beta_phi[1]`) #-0.83 preflood
mean(draws$`beta_phi[2]`) #0.71 OIS
mean(draws$`beta_gam[1]`) #0.36 relief 200m

s <- summary(fit)$summary # I'm looking at rhat and n_eff for convergence and autocorrelation
s # rhats are all generally below 1.01

acf(draws$gamma_0)
acf(draws$psi1_0)
acf(draws$phi_0)
acf(draws$`beta_psi[1]`)
acf(draws$`beta_psi[2]`)
acf(draws$`beta_psi[3]`)
acf(draws$`beta_phi[1]`)
acf(draws$`beta_phi[2]`)
acf(draws$`beta_gam[1]`)
# there is no autocorrelation issues with our data

# ---- Trying to plot observed vs. model ---- #

library(dplyr)

occ_summary <- pika %>%
  group_by(Year) %>%
  summarise(
    mean_occu = mean(Occu),
    se_occu = sd(Occu) / sqrt(n())
  )

library(ggplot2)

ggplot(occ_summary, aes(x = Year, y = mean_occu)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_occu - se_occu,
                    ymax = mean_occu + se_occu),
                width = 0.2) +
  theme_classic() +
  ylab("Mean Occupancy") +
  xlab("Year")

# ---- Model ---- #
library(posterior)
library(dplyr)
library(tidyr)

draws <- as_draws_rvars(fit)   # nicer format

# Extract posterior rvars
psi1   <- draws$psi1           # length nsite
phi    <- draws$phi            # array [site, year]
gamma  <- draws$gamma          # array [site, year]

nsite  <- dim(phi)[1]
nyear  <- dim(phi)[2]


# compute psi dynamics for each site
psi <- array(NA, dim = c(nsite, nyear))

psi1_draws <- posterior::as_draws_matrix(psi1)
psi1_vec <- colMeans(psi1_draws)

phi_draws <- posterior::as_draws_matrix(phi)
phi_mat <- matrix(
  colMeans(phi_draws),
  nrow = nsite,
  byrow = FALSE
)

gamma_draws <- posterior::as_draws_matrix(gamma)
gamma_mat <- matrix(
  colMeans(gamma_draws),
  nrow = nsite,
  byrow = FALSE
)

for (s in 1:nsite) {
  psi[s,1] <- psi1_vec[s]
  for (t in 2:nyear) {
    psi[s,t] <- psi[s,t-1] * phi_mat[s,t-1] +
      (1 - psi[s,t-1]) * gamma_mat[s,t-1]
  }
}


psi_summary <- apply(psi, 2, function(x) {
  mean_val <- mean(x)
  se_val   <- sd(x) / sqrt(length(x))  # standard error
  c(mean = mean_val, se = se_val)
})

# Convert to data frame
psi_df <- data.frame(
  year = 2019:2025,
  mean = psi_summary["mean", ],
  se   = psi_summary["se", ]
)

library(ggplot2)

ggplot(psi_df, aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.3) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  theme_classic() +
  labs(y = "Occupancy Probability (ψ)", x = "Year") 


# ---- Attempting to plot together ---- #
# Observed occupancy summary (already made)
occ_summary2 <- occ_summary %>%
  mutate(type = "Observed") %>%
  rename(mean = mean_occu,
         se   = se_occu)

# Modeled occupancy summary (already made)
psi_df2 <- psi_df %>%
  mutate(type = "Modeled")

plot_df <- bind_rows(
  occ_summary2 %>% select(year = Year, mean, se, type),
  psi_df2
)

library(ggplot2)

ggplot(plot_df, aes(x = year, y = mean, color = type, fill = type)) +
  
  # Modeled SE ribbon
  geom_ribbon(
    data = subset(plot_df, type == "Modeled"),
    aes(ymin = mean - se, ymax = mean + se),
    alpha = 0.25,
    color = NA
  ) +
  
  # Observed SE bars
  geom_errorbar(
    data = subset(plot_df, type == "Observed"),
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.15,
    color = "black",
    linewidth = 1.0
  ) +
  
  # Lines
  geom_line(size = 1.3) +
  
  # Points
  geom_point(size = 3.5) +
  
  # Axis labels
  labs(
    x = "Year",
    y = "Occupancy (ψ)",
    color = "Series",
    fill   = "Series"
  ) +
  
  # Colors (same as before)
  scale_color_manual(values = c("Observed" = "black", "Modeled" = "blue")) +
  scale_fill_manual(values = c("Observed" = "black", "Modeled" = "blue")) +
  
  # Clean theme with *larger, readable text*
  theme_classic(base_size = 18) +   # ⬅️ main improvement
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.position = "right",
    plot.margin = margin(15, 20, 15, 20)
  )


