#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/vax_escape_analyses/sim_config_vax_escape.R')
source('./scripts/two_strain_model/sim_config_global_params.R')


# Controlled epidemic ----

#Simulations with npi and vaccination ====
#All scenarios
vax_escape_model_dynamics <- vax_escape_sim_config_table %>% 
    rowwise() %>% 
    do({with(., 
             simulate_raw_dynamics(model_func = ts_model_vax_escape, 
                                   pop_inits = pop_inits, 
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


vax_escape_dynamics_summaries <- vax_escape_model_dynamics %>% 
    mutate(time, S, wildtype = Iw + Imw + VIw, variant = Im + Iwm + VIm, 
           rec = RmSw + RwSm + R, vaxed = V + RwSmV + RmSwV, K = K, .keep = 'unused')


vax_escape_dynamics_summaries_long <- vax_escape_dynamics_summaries %>% 
    pivot_longer(cols = c(S, wildtype:vaxed), 
                 names_to = 'health_state', 
                 values_to = 'pop') %>% 
    mutate(pop = pop*target_pop)


vax_escape_dynamics_summaries_long %>% 
    filter(health_state == 'wildtype') %>% 
    ggplot() + 
    geom_line(aes(x = time, y = pop, color = health_state)) + 
    facet_wrap('variant_emergence_day')
