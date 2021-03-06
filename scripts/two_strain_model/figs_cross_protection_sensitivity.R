#Packages ----
library(tidyverse)
library(scales)
library(patchwork)
library(mdthemes)

#helper scripts
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')

#plot paths for thesis and git
git_plot_path <- './figures/'
thesis_plot_path <- "C:/Users/JAMESAZAM/Dropbox/My Academic Repository/_SACEMA/Academic/_PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/"

#Load the model output
full_cross_protection_model_summaries <- readRDS('./model_output/orv_npi_all_scenarios_dynamics_parallel.rds') %>% 
    mutate(R0m = R0_m, 
           sigma_w = dynamics_params$sigma_w,
           sigma_m = dynamics_params$sigma_m,
           vax_efficacy_w = 1,
           vax_efficacy_m = 1
    )

partial_cross_protection_model_summaries <- readRDS('./model_output/sensitivity_analyses/cross_protection/model_summaries_partial_cross_protection.rds')


#combine the two model outputs
cross_protection_sensitivity_df <- bind_rows(full_cross_protection_model_summaries, 
                                             partial_cross_protection_model_summaries
) %>% 
    select(-contains('efficacy'))



#Rescale the population proportions to total sizes
cross_protection_sensitivity_pop_rescale <- cross_protection_sensitivity_df %>%
    mutate(across(.cols = c(total_cases, peak_cases, total_vaccinated), 
                  .fns = ~ .x*target_pop #rescale the population proportions to total sizes
    )
    ) %>% 
    select(-c(npi_duration, total_vaccinated))





#' Towards the outbreak size isoclines

outbreak_size_isocline_df <- cross_protection_sensitivity_pop_rescale %>% 
    filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time), 
           npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           total_cases <= 1000
    ) %>% 
    mutate(variant_emergence_day = as_factor(variant_emergence_day)) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity) %>% 
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

# Isoclines ----
# Outbreak size ====
outbreak_size_isocline_cp_sensitivity <- ggplot(outbreak_size_isocline_df, 
                                                aes(x = vax_coverage, 
                                                    y = min_speed, 
                                                    color = variant_emergence_day
                                                )) + 
    geom_line(aes(linetype = factor(sigma_w)), 
              size = 1, 
              show.legend = TRUE
    ) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.30, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(#title = paste('Cumulative cases threshold <= 1000 at various NPI intensity levels'), 
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(outbreak_size_isocline_cp_sensitivity)

#Save the files 
ggsave(outbreak_size_isocline_cp_sensitivity,
       filename = 'outbreak_size_isocline_cp_sensitivity.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

ggsave(outbreak_size_isocline_cp_sensitivity,
       filename = 'outbreak_size_isocline_cp_sensitivity.eps',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')



# Peak prevalence ====

peak_prevalence_isocline_df <- controlled_epidemic_rescaled %>% 
    filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time), 
           npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           peak_cases <= 300
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
    labs(#title = paste('Peak prevalence <= 300 at various NPI intensity levels'), 
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
       filename = 'peak_prevalence_isocline.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

ggsave(peak_prevalence_isocline,
       filename = 'peak_prevalence_isocline.eps',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')





