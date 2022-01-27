#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model_functions.R')
source('./scripts/two_strain_model/simulation_functions.R')
source('./scripts/two_strain_model/vax_escape_analyses/vax_escape_sim_global_params.R')
source('./scripts/two_strain_model/vax_escape_analyses/vax_escape_sim_table_setup.R')



# Controlled epidemic ----

#############################

#Simulate the vax escape model with perfect vax efficacy to mimic old model
#############################
vax_escape_perfect_vax_simulations <- vax_escape_mod_sim_table %>% 
    rowwise() %>% 
    do({with(., 
             simulate_dynamics_vax_escape_ts_model(pop_inits = pop_inits, 
                                   dynamics_parms = dynamics_params,
                                   control_parms = .,
                                   max_time = max_time, 
                                   dt = eval_times,
                                   events_table = event_df,
                                   get_summaries = TRUE, #only extract the summaries 
                                   browse = FALSE
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()


#save the results
saveRDS(object = vax_escape_perfect_vax_simulations, file = './model_output/vax_escape_perfect_efficacy_summaries.rds')

beep(sound = 3)
