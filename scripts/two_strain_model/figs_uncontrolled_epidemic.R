#Packages ----
library(scales)
library(patchwork)
library(tidyverse)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')

#Read the model output
no_control_epidemic_dynamics <- readRDS('./model_output/simulation_uncontrolled_epidemic_dynamics.rds')


#Rescale the population proportions to the actual sizes
no_control_epidemic_rescaled <- no_control_epidemic_dynamics %>%
    mutate(across(.cols = S:K, .fns = ~ .x*target_pop)) #rescale the population proportions to total sizes


ggplot(data = no_control_epidemic_dynamics) + 
    geom_line(aes(x = time, 
                  y = Iw + Im + Iwm + Imw
                  )
              ) +
    scale_y_log10(labels = comma)


#Plot the final size
outbreak_size_no_control_line_plot <- ggplot(data = no_control_epidemic_rescaled) + 
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

#plot(no_control_epidemic_final_size_plot)


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

#plot(no_control_epidemic_peak_cases_plot)

#Plot the two size by size
plot(no_control_epidemic_final_size_plot | no_control_epidemic_peak_cases_plot)