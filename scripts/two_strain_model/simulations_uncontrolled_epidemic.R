#Packages ----
library(deSolve)
library(dplyr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_uncontrolled_epidemic.R')

# Simulations ====

variant_emergence_days <- seq(1, max_time, by = 1) #variant can emerge any day of the epidemic


no_control_epidemic_dynamics <- variant_emergence_days %>% 
    purrr::map_df(function(x){simulate_raw_dynamics(pop_inits = pop_inits, 
                                             dynamics_parms = dynamics_params,
                                             control_parms = cbind(no_control_parms_df, data.frame(variant_emergence_day = x)),
                                             max_time = max_time, 
                                             dt = eval_times,
                                             events_table = event_df,
                                             get_summaries = FALSE,
                                             browse = FALSE)}
                  )


#save the output
saveRDS(no_control_epidemic_dynamics, file = './model_output/uncontrolled_epidemic_dynamics.rds')


