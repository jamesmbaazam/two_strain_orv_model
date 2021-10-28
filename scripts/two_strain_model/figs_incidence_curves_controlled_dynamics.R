#Packages ----
library(scales)
library(patchwork)
library(mdthemes)
library(tidyverse)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')

#Load the model output
controlled_epidemic_dynamics <- readRDS('./model_output/control_dynamics_case_studies.rds')

#Rescale the population proportions to total sizes
controlled_epidemic_dynamics_rescaled <- controlled_epidemic_dynamics %>%
    select(!starts_with(c('S', 'R', 'V', 'K'), ignore.case = FALSE)) %>% 
    mutate(time = time,
           variant_emergence_day = variant_emergence_day,
              incidence = Iw + Iwm + Im + Imw,
              .keep = 'unused'
              ) %>% 
    mutate(across(.cols = incidence, .fns = ~ .x*target_pop))


#The incidence curves
incidence_curves <- ggplot(data = controlled_epidemic_dynamics_rescaled) + 
    geom_line(aes(x = time, 
                  y = incidence,
                  color = variant_emergence_day
                  ), 
              size = 0.1
              ) + 
    facet_wrap(~ as.factor(vax_speed))


print(incidence_curves)
