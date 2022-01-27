#Packages ----
library(scales)
library(patchwork)
library(mdthemes)
library(tidyverse)
library(metR)

#helper scripts
source('./scripts/two_strain_model/simulation_functions.R')
source('./scripts/two_strain_model/vax_escape_analyses/vax_escape_sim_global_params.R')
source('./scripts/two_strain_model/two_strain_model_functions.R')
source('./scripts/two_strain_model/vax_escape_analyses/vax_escape_sim_table_setup.R')


#plot paths for thesis and git
git_plot_path <- './figures/'
thesis_plot_path <- "C:/Users/JAMESAZAM/Dropbox/My Academic Repository/_SACEMA/Academic/_PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/"

#Load the model output
vax_escape_model_perfect_efficacy <- readRDS('./model_output/vax_escape_analyses/vax_escape_perfect_efficacy_summaries.rds')

#Rescale the population proportions to total sizes
output_rescaled <- vax_escape_model_perfect_efficacy %>%
    mutate(across(.cols = c(total_cases, peak_cases, total_vaccinated), 
                  .fns = ~ .x*target_pop #rescale the population proportions to total sizes
    )
    ) %>% 
    select(-c(npi_duration, total_vaccinated)) %>% 
    rename(peak_prevalence = peak_cases)



#' Towards the outbreak size isoclines

outbreak_size_isocline_df <- output_rescaled %>% 
    filter(npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           total_cases <= 1000
           ) %>% 
    mutate(variant_emergence_day = as_factor(variant_emergence_day)) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity) %>% 
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

# Isoclines ----
# Outbreak size ====
outbreak_size_isocline <- ggplot(outbreak_size_isocline_df, 
                                 aes(x = vax_coverage, 
                                     y = min_speed, 
                                     color = variant_emergence_day
                                 )) + 
    geom_line(size = 1, show.legend = TRUE) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.30, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(title = 'Scenarios with cumulative cases <= 1000', 
         subtitle = 'Vaccine escape model with perfect efficacy',
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(outbreak_size_isocline)

#Save the files 
ggsave(outbreak_size_isocline,
       filename = 'outbreak_size_isocline_vax_escape_perfect_efficacy.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

#thesis directory (path on my laptop)
# ggsave(outbreak_size_isocline,
#        filename = 'outbreak_size_isocline_vax_escape_perfect_efficacy.eps',
#        path = thesis_plot_path,
#        width = 23.76,
#        height = 17.86,
#        units = 'cm')



# Peak prevalence ====

peak_prevalence_isocline_df <- output_rescaled %>% 
    filter(npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           peak_prevalence <= 300
           ) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity) %>% 
    mutate(variant_emergence_day = as_factor(variant_emergence_day)) %>% 
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

peak_prevalence_isocline <- ggplot(peak_prevalence_isocline_df, 
                                   aes(x = vax_coverage, 
                                       y = min_speed, 
                                       color = variant_emergence_day
                                   )) + 
    geom_line(size = 1, show.legend = TRUE) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.30, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(title = 'Scenarios with peak prevalence <= 300', 
         subtitle = 'Vaccine escape model with perfect efficacy',
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(peak_prevalence_isocline)


#Save the files 
ggsave(peak_prevalence_isocline,
       filename = 'peak_prevalence_isocline_vax_escape_perfect_efficacy.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

#thesis directory (path on my laptop)
# ggsave(peak_prevalence_isocline,
#        filename = 'peak_prevalence_isocline_vax_escape_perfect_efficacy.eps',
#        path = thesis_plot_path,
#        width = 23.76,
#        height = 17.86,
#        units = 'cm')