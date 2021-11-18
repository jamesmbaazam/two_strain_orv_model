#Packages ----
library(scales)
library(patchwork)
library(mdthemes)
library(tidyverse)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')



#plot paths for thesis and git
git_plot_path <- './figures/'
thesis_plot_path <- "C:/Users/JAMESAZAM/Dropbox/My Academic Repository/_SACEMA/Academic/_PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/"


#Load the model output
case_studies_dynamics_df <- readRDS('./model_output/dynamics_case_studies_controlled_vs_uncontrolled.rds')

#Rescale the population proportions to total sizes
case_studies_dynamics_rescaled <- case_studies_dynamics_df %>%
    group_by(variant_emergence_day, npi_intensity, vax_speed) %>% 
    mutate(incidence = Iw + Iwm + Im + Imw, 
           outbreak_size = max(K), #what is the outbreak size for each emergence day and when does it occur
           peak_time = which.max(incidence)
           ) %>% 
    mutate(across(.cols = c(incidence, outbreak_size), .fns = ~ .x*target_pop)) %>%  
    select(!starts_with(c('S', 'R', 'V', 'K'), ignore.case = FALSE)) %>% 
    ungroup()



#Add a column to label the control scenarios
case_studies_dynamics <- case_studies_dynamics_rescaled %>% 
    mutate(control_type = case_when(npi_intensity == 0 & vax_speed > 0 ~ 'vax_only', 
                                    npi_intensity > 0 ~ 'vax + NPIs',
                                    npi_intensity == 0 & vax_speed == 0 ~ 'no_control'
                                    ),
           control_type = as.factor(control_type)
           )


#Incidence curves ----
#Get the vax only and no control dynamics data
vax_vs_unmitigated_dynamics_df <- case_studies_dynamics %>% 
    filter(control_type %in% c('vax_only', 'no_control'),
           variant_emergence_day %in% c(seq(1, 151, by = 50), max_time)
           ) %>% 
    mutate(vax_speed = ifelse(control_type == 'no_control', NA, vax_speed))


#Vax only vs no control (log scaled) #### 
incidence_curves_vax_only <- ggplot(data = vax_vs_unmitigated_dynamics_df) + 
    geom_line(aes(x = time, 
                  y = incidence,
                  color = as.factor(variant_emergence_day)
                  ), 
              size = 1
              ) + 
    scale_y_log10(labels = comma) +
    facet_wrap(control_type ~ vax_speed, 
               labeller = labeller(control_type = label_value, 
                                   vax_speed = label_both
                                   ), 
               scales = 'fixed'
               ) +
    labs(color = 'Variant emergence',
         x = 'Days',
         y = 'Incidence',
         ) +
    theme_minimal(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(incidence_curves_vax_only)


#Save the plot to git folder
ggsave(plot = incidence_curves_vax_only,
       filename = 'incidence_curves_vax_only.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Save the plot to thesis folder
ggsave(plot = incidence_curves_vax_only,
       filename = 'incidence_curves_vax_only.png',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Get the vax + NPI control and no control dynamics data
vax_and_npi_dynamics <- case_studies_dynamics %>% 
    filter(control_type != 'vax_only',
           variant_emergence_day %in% c(seq(1, 151, by = 50), max_time)
           ) %>% 
    mutate(vax_speed = ifelse(control_type == 'no_control', NA, vax_speed),
           npi_intensity = ifelse(control_type == 'no_control', NA, npi_intensity)
           )


#Vax + NPI vs no control (log scaled) ####
incidence_curves_vax_and_npi <- ggplot(data = vax_and_npi_dynamics) + 
    geom_line(aes(x = time, 
                  y = incidence,
                  color = as.factor(variant_emergence_day)
    ), 
    size = 1
    ) + 
    scale_y_log10(labels = comma) +
    facet_wrap(control_type ~ vax_speed + npi_intensity, 
               labeller = labeller(control_type = label_value, 
                                   vax_speed = label_both,
                                   npi_intensity = label_both
                                   ), 
               scales = 'fixed'
    ) +
    labs(color = 'Variant emergence',
         x = 'Days',
         y = 'Incidence',
    ) +
    theme_minimal(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(incidence_curves_vax_and_npi)


#Save the plot to git folder
ggsave(plot = incidence_curves_vax_and_npi,
       filename = 'incidence_curves_vax_and_npi.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Save the plot to thesis folder
ggsave(plot = incidence_curves_vax_and_npi,
       filename = 'incidence_curves_vax_and_npi.png',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )




#Subset the unmitigated epidemic dynamics
dynamics_no_control <- case_studies_dynamics %>% 
    filter(control_type == 'no_control')


#Incidence curves (no control) ----
incidence_curve_unmitigated <- ggplot(data = dynamics_no_control %>% 
                                          filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time))
                                      ) +
    geom_line(aes(x = time, 
                  y = incidence,
                  color = as.factor(variant_emergence_day)
                  ),
              size = 1
    ) +
    scale_y_continuous(labels = comma) +
    labs(x = 'Time (days)', y = 'Incidence', color = 'Variant emergence day') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(incidence_curve_unmitigated)

#Save the plot to git folder
ggsave(plot = incidence_curve_unmitigated,
       filename = 'incidence_curve_unmitigated.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Save the plot to thesis folder
ggsave(plot = incidence_curve_unmitigated,
       filename = 'incidence_curve_unmitigated.png',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Peak timing versus variant emergence day (no control) ----
peak_timing_vs_emergence_unmitigated <- ggplot(data = dynamics_no_control) +
    geom_line(aes(x = variant_emergence_day, 
                  y = peak_time),
              size = 1
              ) +
    scale_y_continuous(breaks = seq(1, 200, 10), 
                       labels = seq(1, 200, 10)
                       ) +
    labs(x = 'Variant emergence day', y = 'Timing of peak') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 
   

print(peak_timing_vs_emergence_unmitigated)

#Save the plot to git folder
ggsave(plot = peak_timing_vs_emergence_unmitigated,
       filename = 'peak_timing_vs_emergence_unmitigated.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Save the plot to thesis folder
ggsave(plot = peak_timing_vs_emergence_unmitigated,
       filename = 'peak_timing_vs_emergence_unmitigated.png',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )
