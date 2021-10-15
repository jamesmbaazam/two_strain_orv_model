#Packages ----
library(dplyr)

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
vax_cov <- c(0.5, 1) #various levels of coverage
vax_speed_scenarios <- seq(1, 10, length.out = 5) #how many times faster than the daily vaccination rate
daily_vax_rate <- as.vector(sapply(as.list(vax_cov/(max_time - vax_start)), function(x) {x*vax_speed_scenarios})) #daily rate of achieving the same coverage by the end of the period
#daily_vax_rate_vec <- daily_vaccinated_prop*target_pop

#Find the time it will take to achieve the assumed coverate with the given rates
campaign_controls_df <- data.frame(vax_rate = daily_vax_rate,
                                   vax_speed = rep(vax_speed_scenarios, 2),
                                   vax_coverage = rep(vax_cov, each = 5)
                                   )

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
    
