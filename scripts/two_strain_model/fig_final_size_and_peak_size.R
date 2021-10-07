#Packages ----
library(scales)
library(patchwork)
library(tidyverse)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')

#Read the model output
no_control_epidemic <- readRDS('./model_output/simulation_uncontrolled_epidemic.rds')


#Rescale the population proportions to the actual sizes
no_control_epidemic_rescaled <- no_control_epidemic %>%
    mutate(total_cases = total_cases*target_pop,
           peak_cases = peak_cases*target_pop) %>% 
    select(-c(npi_start, npi_intensity, npi_duration, 
              vax_start, campaign_duration, vax_coverage))


#Plot the final size
no_control_epidemic_final_size_plot <- ggplot(data = no_control_epidemic_rescaled) + 
    geom_line(aes(x = variant_emergence_day, 
                  y = total_cases
    ),
    size = 1.2
    ) +
    scale_y_continuous(labels = comma) +
    labs(x = 'Variant emergence day', 
         y = 'Total cases'
    ) +
    theme_minimal(base_size = 12)

plot(no_control_epidemic_final_size_plot)


#Plot the peak cases
no_control_epidemic_peak_cases_plot <- ggplot(data = no_control_epidemic_rescaled) + 
    geom_line(aes(x = variant_emergence_day, 
                  y = peak_cases
    ),
    size = 1.2
    ) +
    scale_y_continuous(labels = comma) +
    labs(x = 'Variant emergence day', 
         y = 'Peak daily incidence'
    ) +
    theme_minimal(base_size = 12)

plot(no_control_epidemic_peak_cases_plot)

#Plot the two size by size
plot(no_control_epidemic_final_size_plot | no_control_epidemic_peak_cases_plot)