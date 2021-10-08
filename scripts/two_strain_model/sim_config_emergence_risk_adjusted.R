#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)

# Helper scripts ----
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/two_strain_model.R')


# Uncontrolled epidemic simulation ----

# Controlled epidemic simulations set up ---- 

# Event times ====
# NPI ####
# Best case scenario ####
# npi_start <- 1:max_time
# npi_duration <- max_time - npi_start
# npi_intensity <- seq(0.2, 1, 0.2) #Five levels of npi intensity (could correspond to the five stages in South Africa, for e.g)

#Simple case: No NPI 
npi_start <- max_time
npi_duration <- 0
npi_intensity <- 0 #Five levels of npi intensity (could correspond to the five stages in South Africa, for e.g)


# Vaccination ====
# Best case scenario ####
# vax_start <- 1:max_time
# campaign_duration <- max_time/2 - vax_start
# vax_cov <- seq(0.2, 1, 0.20)

# Simple case: vaccination starts on day 1 and can achieve 100% coverage but at different rates
vax_start <- 1
vax_cov <- 1
vax_rate_vec <- seq(0.1, 1, 0.05)

#Find the time it will take to achieve the assumed coverate with the given rates
campaign_controls_df <- vax_rate_vec %>% 
    purrr::map_df(function(x) {
        calc_campaign_duration(vax_rate = x, vax_coverage = vax_cov, coverage_correction = 0.9999)
    })

# Full simulation table with variant emergence times appended ====
simulation_table <- campaign_controls_df %>%
    slice(rep(1:n(), times = length(variant_emergence_times))) %>% 
    mutate(variant_emergence_day = rep(variant_emergence_times, 
                                       each = nrow(campaign_controls_df)
                                       ),
           vax_start = 1,
           npi_start = npi_start,
           npi_intensity = npi_intensity,
           npi_duration = npi_duration
           )


#Baseline; no variant emerges and only wild type prevails
baseline_params <- campaign_controls_df %>%
    mutate(variant_emergence_day = max_time, 
           vax_start = 1,
           npi_start = npi_start,
           npi_intensity = npi_intensity,
           npi_duration = npi_duration
           )
    
