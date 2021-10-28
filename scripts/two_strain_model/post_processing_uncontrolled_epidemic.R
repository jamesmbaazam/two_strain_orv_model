#Packages ----
library(tidyverse)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_uncontrolled_epidemic.R')

#Read the model output
no_control_epidemic_dynamics <- readRDS('./model_output/uncontrolled_epidemic_dynamics.rds')


#Rescale the population proportions to the actual sizes
no_control_epidemic_rescaled <- no_control_epidemic_dynamics %>%
    mutate(across(.cols = S:K, .fns = ~ .x*target_pop)) #rescale the population proportions to total sizes



#Extract the summaries
no_control_epidemic_summaries <- no_control_epidemic_rescaled %>% 
    group_split(variant_emergence_day) %>% 
    purrr::map_dfr(function(df){extract_model_summaries(df)}) %>% 
    relocate(variant_emergence_day, .before = npi_intensity) %>% 
    as_tibble()



outbreak_size_max <- no_control_epidemic_summaries %>% 
    slice_max(order_by = total_cases) %>% 
    select(-c(peak_cases, total_vaccinated))


outbreak_size_min <- no_control_epidemic_summaries %>% 
    slice_min(order_by = total_cases) %>% 
    select(-c(peak_cases, total_vaccinated))

peak_incidence_max <- no_control_epidemic_summaries %>% 
    slice_max(order_by = peak_cases) %>% 
    select(-c(total_cases, total_vaccinated))


peak_incidence_min <- no_control_epidemic_summaries %>% 
    slice_min(order_by = peak_cases) %>% 
    select(-c(total_cases, total_vaccinated))


#Range of outbreak size and the parameters producing them
range(no_control_epidemic_summaries$total_cases)
#df
outbreak_size_range <- bind_rows(outbreak_size_min, outbreak_size_max)
outbreak_size_range

#Range of peak incidence and the parameters producing them
range(no_control_epidemic_summaries$peak_cases)

#df
peak_incidence_range <- bind_rows(peak_incidence_min, peak_incidence_max)
peak_incidence_range





