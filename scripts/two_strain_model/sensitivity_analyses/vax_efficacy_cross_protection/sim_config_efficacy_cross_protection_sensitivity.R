# Load packages
library(dplyr)

#Load helper scripts
source('./scripts/two_strain_model/two_strain_model.R')


# Run the model for this number of days
max_time <- 365

# Show the model results at these time points
eval_times <- seq(0, max_time, 1)

# When the variant will emerge
# variant_emergence_day <- 60
# variant_emergence_day_vec <- seq(60, max_time, 1)

# Population parameters
target_pop <- 50E6 # target population size
Iw_index_cases <- 50 # Initial/Index wild type cases
Im_index_cases <- 50 # Initial/Index variant cases


# Population initial conditions
pop_inits <- c(
  S = 1 - Iw_index_cases / target_pop,
  Iw = Iw_index_cases / target_pop,
  Im = 0,
  Iwm = 0,
  Imw = 0,
  RwSm = 0,
  RmSw = 0,
  R = 0,
  V = 0,
  VIw = 0,
  VIm = 0,
  RwSmV = 0,
  RmSwV = 0,
  K = 0
)


# ===============================
# Parameters for dynamics
# ===============================
R0_w <- 2.0 # R0 (wild type)
R0_m <- R0_w * 1.3 # R0 (variant) 30% more infectious
IP_w <- 14 # infectious period (wild type)
IP_m <- 14 # infectious period (variant)


# varying assumptions about vax efficacy and cross protection
efficacy_cross_protection_df <- data.frame(
  vax_efficacy_w = c(1, 0.95), # Vaccine efficacy against the wild-type
  vax_efficacy_m = c(0.80), # Vaccine efficacy against the variant)
  sigma_w = 0.85, # Cross-protection provided by wildtype
  sigma_m = 0.85 # Cross-protection provided by variant
)


dynamics_params <- efficacy_cross_protection_df %>%
  mutate(
    beta_w = R0_w / IP_w,
    beta_m = R0_m / IP_m,
    gamma_w = 1 / IP_w,
    gamma_m = 1 / IP_m
  )




# ===============================
# Parameters for control dynamics
# ===============================
coverage_correction <- 0.9991



#' Event data frame for introducing the variant into the model dynamics
#' (check ?deSolve::event for more on the structure of the event_df below)
#'

# Variant emergence times ===
variant_emergence_times <- c(1, 61, 121, 151, max_time)

event_df <- data.frame(
  var = c("S", "Im"), # Compartments to change at a set time
  value = c(-Im_index_cases / target_pop, Im_index_cases / target_pop), # introduce 10 variant cases
  method = c("add", "replace")
) # operation on state variables


#############################################################################
#' Intervention set up
#############################################################################


#Non-pharma interventions (NPIs)
npi_start <- 1 #Day of commencement in the model
npi_duration <- max_time - npi_start
npi_intensity <- seq(0, 0.3, by = 0.02) #Various levels of npi intensity 


# Vaccination ====
# Best case scenario ####
# Simple case: vaccination starts on day 1 and can achieve 100% coverage but at different rates
vax_start <- 1 #Day of commencement of vaccination in the model
vax_cov <- seq(0.1, 1, by = 0.025) #various levels of vaccination coverage
vax_speed_scenarios <- seq(1, 10, by = 0.25) #how many times faster than the daily vaccination rate
#daily_vax_rate <- as.vector(sapply(as.list(vax_cov/(max_time - vax_start)), function(x) {x*vax_speed_scenarios})) #daily rate of achieving the same coverage by the end of the period

#Find the time it will take to achieve the assumed coverage with the given rates
campaign_controls_df <- pmap_dfr(list(vax_cov, max_time, vax_start), 
                                 function(vax_cov, max_time, vax_start){
                                     calc_vax_rate(vax_coverage = vax_cov, 
                                                   coverage_correction = coverage_correction, 
                                                   max_time = max_time, 
                                                   campaign_start = vax_start
                                     )
                                 }) 

#Convert the calculate vaccination rates (hazards) to various speeds
campaign_controls_scenarios_df <- vax_speed_scenarios %>% 
    map_dfr(function(vax_speed) {
        campaign_controls_df %>% 
            group_by(vax_coverage) %>% 
            mutate(vax_rate = vax_rate*vax_speed, 
                   vax_speed = vax_speed,
                   vax_coverage = vax_coverage, 
                   campaign_duration = campaign_duration/vax_speed
            ) 
    }
    ) %>% 
    ungroup() %>% 
    arrange(vax_coverage)


campaign_controls_scenarios_df <- campaign_controls_scenarios_df %>% 
    slice(rep(1:n(), times = length(npi_intensity))) %>% 
    mutate(npi_intensity = rep(npi_intensity, each = nrow(campaign_controls_scenarios_df)))


# Full simulation table with variant emergence times appended ====
sim_config_table_full <- campaign_controls_scenarios_df %>%
    slice(rep(1:n(), times = length(variant_emergence_times))) %>% 
    mutate(variant_emergence_day = rep(variant_emergence_times, 
                                       each = nrow(campaign_controls_scenarios_df)
    ),
    vax_start = 1,
    npi_start = npi_start,
    npi_duration = npi_duration
    )


#Expand the simulation table to include the two scenarios of vaccine efficacy
perfect_efficacy_scenario_df <- sim_config_table_full %>% mutate(filter(dynamics_params, vax_efficacy_w == 1))
imperfect_efficacy_scenario_df <- sim_config_table_full %>% mutate(filter(dynamics_params, vax_efficacy_w != 1))

#merge the two data.frames
simulation_config_table <- bind_rows(perfect_efficacy_scenario_df, imperfect_efficacy_scenario_df)
