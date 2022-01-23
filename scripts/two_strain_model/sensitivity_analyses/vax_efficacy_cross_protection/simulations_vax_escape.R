#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sensitivity_analyses/vax_efficacy_cross_protection/sim_config_efficacy_cross_protection_sensitivity.R')



# Controlled epidemic with vax escape ----

#Simulations with npi and vaccination ====
#All scenarios
vax_escape_sensitivity_analysis_summaries <- simulation_config_table %>% 
    slice(1:2) %>% 
    rowwise() %>% 
    do({with(., 
             simulate_raw_dynamics(model_func = ts_model_vax_escape, 
                                   get_summaries_func = extract_vax_escape_mod_summaries, 
                                   pop_inits = pop_inits, 
                                   dynamics_parms = NULL,
                                   control_parms = .,
                                   max_time = max_time, 
                                   dt = eval_times,
                                   events_table = event_df,
                                   get_summaries = F,
                                   browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()

#play a sound when the simulation is over
beep(sound = 3)