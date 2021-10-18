#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')


# Controlled epidemic

#Baseline; no variant emerges and only wild type prevails
baseline_no_variant <- baseline_params %>% 
    rowwise() %>% 
    do({with(., 
             simulate_raw_dynamics(pop_inits = pop_inits, 
                            dynamics_parms = dynamics_params,
                            control_parms = .,
                            max_time = max_time, 
                            dt = eval_times,
                            events_table = event_df,
                            browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()  


#Simulations with variant emergence
orv_full_simulation <- simulation_table %>% 
    rowwise() %>% 
    do({with(., 
             simulate_raw_dynamics(pop_inits = pop_inits, 
                            dynamics_parms = dynamics_params,
                            control_parms = .,
                            max_time = max_time, 
                            dt = eval_times,
                            events_table = event_df,
                            browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()

#save the simulation
saveRDS(object = orv_full_simulation, file = './model_output/controlled_epidemic_dynamics_sequential_run.rds')

beep()


