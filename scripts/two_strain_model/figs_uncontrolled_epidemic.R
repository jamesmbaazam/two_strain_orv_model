#Packages ----
library(scales)
library(patchwork)
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
    purrr::map_dfr(function(df){extract_model_summaries(df, no_control_parms_df)}) %>% 
    relocate(variant_emergence_day, .before = npi_intensity) %>% 
    as_tibble()



#Plot the incidence on log scale
incidence_curve_line_plot <- ggplot(data = no_control_epidemic_rescaled) + 
    geom_line(aes(x = time, 
                  y = Iw + Im + Iwm + Imw
                  )
              ) +
    scale_y_log10(labels = comma) +
    labs(title = 'Incidence curve for the uncontrolled epidemic', 
         x = 'Variant emergence day',
         y = 'Incidence (log-transformed)'
         )

#Plot in the viewer
print(incidence_curve_line_plot)

#Save the plot to file
ggsave(plot = incidence_curve_line_plot,
       filename = './figures/incidence_curve_line_plot.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )

#Plot the outbreak size
outbreak_size_no_control_line_plot <- ggplot(data = no_control_epidemic_summaries) + 
    geom_line(aes(x = variant_emergence_day, 
                  y = total_cases
    ),
    size = 1.2
    ) +
    scale_y_log10(labels = unit_format(unit = 'M', scale = 1E-6)) +
    labs(title = 'Outbreak sizes for the uncontrolled epidemic',
         x = 'Variant emergence day', 
         y = 'Outbreak size (log-transformed)'
    ) +
    theme_minimal(base_size = 12)

#Plot in the viewer
print(outbreak_size_no_control_line_plot)

#Save the plot to file
ggsave(plot = outbreak_size_no_control_line_plot,
       filename = './figures/outbreak_size_no_control_line_plot.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Plot the peak cases
peak_incidence_no_control_line_plot <- ggplot(data = no_control_epidemic_summaries) + 
    geom_line(aes(x = variant_emergence_day, 
                  y = peak_cases
    ),
    size = 1.2
    ) +
    scale_y_log10(labels = unit_format(unit = 'M', scale = 1E-6)) +
    labs(title = 'Peak incidence for the uncontrolled epidemic', 
         x = 'Variant emergence day', 
         y = 'Peak incidence (log-transformed)'
    ) +
    theme_minimal(base_size = 12)

print(peak_incidence_no_control_line_plot)


#Save the plot to file
ggsave(plot = peak_incidence_no_control_line_plot,
       filename = './figures/peak_incidence_no_control_line_plot.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
)