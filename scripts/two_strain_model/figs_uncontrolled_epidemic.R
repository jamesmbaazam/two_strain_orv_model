#Packages ----
library(dplyr)
library(purrr)
library(ggplot2)
library(scales)
library(patchwork)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_uncontrolled_epidemic.R')
source('./scripts/two_strain_model/two_strain_model_functions.R')
source('./scripts/two_strain_model/simulation_functions.R')

#Read the model output
no_control_epidemic_dynamics <- readRDS('./model_output/uncontrolled_epidemic_dynamics.rds')


#plot paths for thesis and git
git_plot_path <- './figures/'
thesis_plot_path <- "C:/Users/JAMESAZAM/Dropbox/My Academic Repository/_SACEMA/Academic/_PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/"


#Rescale the population proportions to the actual sizes
no_control_epidemic_rescaled <- no_control_epidemic_dynamics %>%
    mutate(across(.cols = S:K, .fns = ~ .x*target_pop)) #rescale the population proportions to total sizes



#Extract the summaries
no_control_epidemic_summaries <- no_control_epidemic_rescaled %>% 
    group_split(variant_emergence_day) %>% 
    purrr::map_dfr(function(df){extract_summaries_ts_model(dynamics_df = df, 
                                                           browse = FALSE
                                                           )
        }
        ) %>% 
    relocate(variant_emergence_day, .before = npi_intensity) %>% 
    as_tibble()



#Plot the incidence on log scale
# incidence_curve_line_plot <- ggplot(data = no_control_epidemic_rescaled) + 
#     geom_line(aes(x = time, 
#                   y = Iw + Im + Iwm + Imw
#                   )
#               ) +
#     scale_y_log10(labels = comma) +
#     labs(title = 'Incidence curve for the uncontrolled epidemic', 
#          x = 'Variant emergence day',
#          y = 'Incidence (log-transformed)'
#          )

#Plot in the viewer
# print(incidence_curve_line_plot)

#Save the plot to file
# ggsave(plot = incidence_curve_line_plot,
#        filename = './figures/incidence_curve_line_plot.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
#        )

#Plot the outbreak size
outbreak_size_no_control_line_plot <- ggplot(data = no_control_epidemic_summaries %>% 
                                                 filter(variant_emergence_day <= 365)
                                             ) + 
    geom_line(aes(x = variant_emergence_day, 
                  y = total_cases
    ),
    size = 1.2
    ) +
    scale_y_log10(labels = unit_format(unit = 'M', scale = 1E-6), 
                  breaks = seq(35E6, 50E6, 5E5)
                  ) +
    labs(#title = 'Outbreak sizes for the uncontrolled epidemic',
         x = 'Variant emergence day', 
         y = 'Outbreak size (log-transformed)'
    ) +
    theme_minimal(base_size = 18)

#Plot in the viewer
print(outbreak_size_no_control_line_plot)

#Save the plot to git folder
ggsave(plot = outbreak_size_no_control_line_plot,
       filename = 'outbreak_size_no_control_line_plot.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Save the plot to thesis folder
ggsave(plot = outbreak_size_no_control_line_plot,
       filename = 'outbreak_size_no_control_line_plot.png',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
)




#Plot the peak cases
peak_incidence_no_control_line_plot <- ggplot(data = no_control_epidemic_summaries %>% 
                                                  filter(variant_emergence_day <= 365)
                                              ) + 
    geom_line(aes(x = variant_emergence_day, 
                  y = peak_cases
    ),
    size = 1.2
    ) +
    scale_y_log10(labels = unit_format(unit = 'M', scale = 1E-6), 
                  breaks = seq(6.5E6, 13E6, 5E5)
                  ) +
    labs(#title = 'Peak incidence for the uncontrolled epidemic', 
         x = 'Variant emergence day', 
         y = 'Peak prevalence (log-transformed)'
    ) +
    theme_minimal(base_size = 18)

print(peak_incidence_no_control_line_plot)


#Save the plot to the git folder
ggsave(plot = peak_incidence_no_control_line_plot,
       filename = 'peak_incidence_no_control_line_plot.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
)

#Save the plot to the thesis folder
ggsave(plot = peak_incidence_no_control_line_plot,
       filename = 'peak_incidence_no_control_line_plot.png',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
)
