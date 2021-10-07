#Packages ----
library(deSolve)
library(dplyr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')

# control parameters (turned off for uncontrolled epidemic)
no_control_parms_df <- data.frame(npi_intensity = 0, vax_coverage = 0, 
                                  vax_start = 1, campaign_duration = 0, 
                                  npi_start = 1, npi_duration = 0
                                  )
# Simulations ====
no_control_epidemic <- variant_emergence_times %>% 
    purrr::map_df(function(x){simulate_model(pop_inits = pop_inits, 
                                             dynamics_parms = dynamics_params,
                                             control_parms = cbind(no_control_parms_df, data.frame(variant_emergence_day = x)),
                                             max_time = max_time, 
                                             dt = eval_times,
                                             events_table = event_df,
                                             return_dynamics = FALSE,
                                             browse = FALSE)}
                  )


#save the output
saveRDS(no_control_epidemic, file = './model_output/simulation_uncontrolled_epidemic.rds')

