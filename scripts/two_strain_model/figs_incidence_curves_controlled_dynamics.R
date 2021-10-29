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
controlled_epidemic_dynamics <- readRDS('./model_output/control_dynamics_case_studies.rds')

#Rescale the population proportions to total sizes
controlled_epidemic_dynamics_rescaled <- controlled_epidemic_dynamics %>%
    select(!starts_with(c('S', 'R', 'V', 'K'), ignore.case = FALSE)) %>% 
    mutate(time = time,
           variant_emergence_day = variant_emergence_day,
              incidence = Iw + Iwm + Im + Imw,
              .keep = 'unused'
              ) %>% 
    mutate(across(.cols = incidence, .fns = ~ .x*target_pop)) %>% 
    group_by(variant_emergence_day) %>% 
    mutate(outbreak_size = max(incidence), #what is the outbreak size for each emergence day and when does it occur
           peak_time = which.max(incidence)
           ) %>% 
    ungroup()


#The incidence curves
controlled_epidemic_inc_curves <- ggplot(data = controlled_epidemic_dynamics_rescaled %>% 
                               filter(variant_emergence_day %in% c(31, 365)) %>% 
                                   mutate(emergence_scenario = as.factor(ifelse(variant_emergence_day == 31, #Scenarios for emergence on day 31 versus no emergence
                                                                      'day_31', 
                                                                      'no_emergence'
                                                                      ))
                                          )
                           ) + 
    geom_line(aes(x = time, 
                  y = incidence,
                  color = emergence_scenario,
                  linetype = emergence_scenario
                  ), 
              size = 1
              ) + 
    scale_y_log10(labels = unit_format(unit = 'M', scale = 1E-6)) +
    facet_wrap(npi_intensity ~ vax_speed, labeller = 'label_both', scales = 'free') +
    labs(color = 'Variant emergence',
         linetype = 'Variant emergence',
         x = 'Days',
         y = 'Incidence (log-transformed)',
         ) +
    theme_minimal(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold')) 


print(controlled_epidemic_inc_curves)


#Save the plot to git folder
ggsave(plot = controlled_epidemic_inc_curves,
       filename = 'controlled_epidemic_inc_curves.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Save the plot to thesis folder
ggsave(plot = controlled_epidemic_inc_curves,
       filename = 'controlled_epidemic_inc_curves.png',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )
