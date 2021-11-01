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
    mutate(across(.cols = c(incidence, outbreak_size), .fns = ~ .x*target_pop),
           with_control = ifelse(vax_coverage == 0 & npi_intensity == 0, 'uncontrolled', 'controlled')) %>% 
    select(!starts_with(c('S', 'R', 'V', 'K'), ignore.case = FALSE)) %>% 
    ungroup()


#The incidence curves
incidence_curves <- ggplot(data = case_studies_dynamics_rescaled %>% 
                               filter(variant_emergence_day %in% c(31, 365)) %>% 
                                   mutate(emergence_scenario = as.factor(ifelse(variant_emergence_day == 31, #Scenarios for emergence on day 31 versus no emergence
                                                                      'day_31', 
                                                                      'no_emergence'
                                                                      ))
                                          )
                           ) + 
    geom_line(aes(x = time, 
                  y = incidence,
                  color = with_control,
                  linetype = emergence_scenario
                  ), 
              size = 1
              ) + 
    scale_y_continuous(labels = comma) +
    facet_wrap(npi_intensity ~ vax_speed, labeller = 'label_both', scales = 'free') +
    labs(color = 'Control scenario',
         linetype = 'Variant emergence',
         x = 'Days',
         y = 'Incidence',
         ) +
    theme_minimal(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold')) 


print(incidence_curves)


#Save the plot to git folder
ggsave(plot = incidence_curves,
       filename = 'incidence_curves.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Save the plot to thesis folder
ggsave(plot = incidence_curves,
       filename = 'incidence_curves.png',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )



#Summaries
case_studies_dynamics_rescaled %>% 
    group_by(variant_emergence_day, npi_intensity, vax_coverage, vax_speed) %>% 
    summarise(peak_incidence = max(incidence),
              outbreak_size = max(outbreak_size)
              ) %>% 
    ungroup() %>% 
    # filter(vax_coverage == 0.0) %>% 
    arrange(outbreak_size)
