#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sensitivity_analyses/sim_config_global_params_sensitivity.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')


# Controlled epidemic ----

#Simulations with npi and vaccination ====
#All scenarios
orv_npi_simulation_all_scenarios <- orv_npi_control_config_table %>% 
    rowwise() %>% 
    do({with(., 
             simulate_raw_dynamics(pop_inits = pop_inits, 
                                   dynamics_parms = dynamics_params,
                                   control_parms = .,
                                   max_time = max_time, 
                                   dt = eval_times,
                                   events_table = event_df,
                                   get_summaries = TRUE,
                                   browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()




#NPI = 0% ####
orv_npi_simulation_npi0 <- orv_npi_control_config_table %>% 
    filter(npi_intensity == 0) %>% #npi levels of 0%
    rowwise() %>% 
    do({with(., 
             simulate_raw_dynamics(pop_inits = pop_inits, 
                            dynamics_parms = dynamics_params,
                            control_parms = .,
                            max_time = max_time, 
                            dt = eval_times,
                            events_table = event_df,
                            get_summaries = TRUE,
                            browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()


#NPI = 20% ####
orv_npi_simulation_npi20 <- orv_npi_control_config_table %>% 
    filter(npi_intensity == 0.2) %>% #npi levels of 0%, 10%, and 50%
    rowwise() %>% 
    do({with(., 
             simulate_raw_dynamics(pop_inits = pop_inits, 
                                   dynamics_parms = dynamics_params,
                                   control_parms = .,
                                   max_time = max_time, 
                                   dt = eval_times,
                                   events_table = event_df,
                                   get_summaries = TRUE,
                                   browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()

#NPI = 50% ####
orv_npi_simulation_npi50 <- orv_npi_control_config_table %>% 
    filter(npi_intensity == 0.5) %>% #npi levels of 0%, 10%, and 50%
    rowwise() %>% 
    do({with(., 
             simulate_raw_dynamics(pop_inits = pop_inits, 
                                   dynamics_parms = dynamics_params,
                                   control_parms = .,
                                   max_time = max_time, 
                                   dt = eval_times,
                                   events_table = event_df,
                                   get_summaries = TRUE,
                                   browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()



#save the simulations
saveRDS(object = orv_full_simulation, file = './model_output/controlled_epidemic_dynamics_sequential_run_npi0.rds')
saveRDS(object = orv_full_simulation, file = './model_output/controlled_epidemic_dynamics_sequential_run_npi20.rds')
saveRDS(object = orv_full_simulation, file = './model_output/controlled_epidemic_dynamics_sequential_run_npi50.rds')

beep(sound = 3)


