"dpdata.csv" is the raw dataset used in creating the model for Pika occupancy.

"dyn_occ_real_data.R" is the R code for setting up the data and executing the model through R. It also includes code for running diagnostics on the model after it has been run.

"dyn_occ_V3.stan" Is the stan file which contains the bayesian dynamic occupancy model that gets run through R.

The other R code scripts are supplemental code for plotting trend plots on phi, gamma and psi. It also includes code for running partial effects plots for each variable. I recommend using "Flood_Vs_Persistence.R", "Partial_Effects_Plots.R", and "Modeled_Phi_Vs_Gamma.R". These are the most up-to-date code for plotting aspects of the model.
