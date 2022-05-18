#Packages ----
library(dplyr)
library(purrr)

# Helper scripts ----
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/two_strain_model_functions.R')

# Event times ====
# NPI ####
# Best case scenario ####
# npi_start <- 1:max_time
# npi_duration <- max_time - npi_start
# npi_intensity <- seq(0.2, 1, 0.2) #Five levels of npi intensity (could correspond to the five stages in South Africa, for e.g)

#Simple case: No NPI 
npi_start <- 1
npi_duration <- max_time - npi_start
npi_intensity <- seq(0, 0.3, by = 0.02) #Five levels of npi intensity (could correspond to the five stages in South Africa, for e.g)


# Vaccination ====
# Best case scenario ####
# vax_start <- 1:max_time
# campaign_duration <- max_time/2 - vax_start
# vax_cov <- seq(0.2, 1, 0.20)

# Simple case: vaccination starts on day 1 and can achieve 100% coverage but at different rates
vax_start <- 1
vax_cov <- seq(0.1, 1, by = 0.025) #various levels of coverage
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
                                        }
                                        ) 

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
orv_npi_control_config_table <- campaign_controls_scenarios_df %>%
    slice(rep(1:n(), times = length(variant_emergence_times))) %>% 
    mutate(variant_emergence_day = rep(variant_emergence_times, 
                                       each = nrow(campaign_controls_scenarios_df)
                                       ),
           vax_start = 1,
           npi_start = npi_start,
           npi_duration = npi_duration
           )




