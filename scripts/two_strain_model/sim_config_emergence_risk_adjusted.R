#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')


# Uncontrolled epidemic simulation ----

# Controlled epidemic simulations set up ---- 

# Event times ====
# NPI ####
# Best case scenario ####
npi_start <- 1:max_time
npi_duration <- max_time - npi_start
npi_intensity <- seq(0.2, 1, 0.2) #Five levels of npi intensity (could correspond to the five stages in South Africa, for e.g)

# Vaccination ====
# Best case scenario ####
vax_start <- 1:max_time
campaign_duration <- max_time - vax_start
vax_cov <- seq(0.2, 1, 0.20)

# Complete set of scenarios ====
# All combinations of vaccination coverage and NPI intensity  
control_scenarios <- expand.grid(npi_intensity = npi_intensity, vax_coverage = vax_cov) 

control_scenarios_rep <- control_scenarios %>% 
    slice(rep(1:n(), times = length(variant_emergence_times))) # repeatedly bind the df to itself


# Repeat the variant emergence times for binding
variant_emergence_times_df <- data.frame(variant_emergence_day = rep(variant_emergence_times, each = nrow(control_scenarios))) 

control_and_emergence_scenarios <- cbind(control_scenarios_rep, variant_emergence_times_df)

control_and_emergence_scenarios_rep <- control_and_emergence_scenarios %>% 
    slice(rep(1:n(), times = length(vax_start))) # repeatedly bind the df to itself
    

# Full simulation table ====
simulation_table <- control_and_emergence_scenarios_rep %>% 
    mutate(campaign_duration = campaign_duration,
           npi_duration = npi_duration)

