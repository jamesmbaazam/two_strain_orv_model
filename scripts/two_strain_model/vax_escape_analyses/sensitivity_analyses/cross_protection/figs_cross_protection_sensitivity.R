#Packages ----
library(tidyverse)
library(scales)
library(patchwork)
library(mdthemes)

#helper scripts
source('./scripts/two_strain_model/two_strain_model_functions.R')
source('./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/global_params.R')
source('./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/intervention_params_table_setup_global.R')

#plot paths for thesis and git
git_plot_path <- './figures/'
thesis_plot_path <- "C:/Users/JAMESAZAM/Dropbox/My Academic Repository/_SACEMA/Academic/_PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/"

#Load the model output

# Original analysis
#' The baseline analysis assumes perfect vaccine efficacy against the wild-type 
#' and variant and perfect cross protection between the strains. The variant is
#' assumed to be 30% more transmissible than the wild-type (R0w = 1.3R0m)
baseline_analysis_results <- readRDS('./model_output/sensitivity_analyses/baseline_analysis_vR0_30_percent/vR0_30_percent_baseline_analysis_summaries.rds') %>% 
    rename(peak_prevalence = peak_cases)

#Cross protection sensitivity analysis results
cp_sensitivity_analysis_results <- readRDS('./model_output/sensitivity_analyses/cross_protection/cp_sensitivity_summaries.rds')


#combine the two model outputs
cp_sensitivity_analysis_all_results <- bind_rows(baseline_analysis_results, cp_sensitivity_analysis_results) 

#Remove redundant columns
cp_sensitivity_analysis_df <- cp_sensitivity_analysis_all_results %>% 
    select(-c(vax_rate, 
              npi_duration, 
              starts_with('vax_efficacy'), 
              cross_protection_m,
              R0m,
              total_vaccinated
              )
           )

#Rescale the population proportions to total sizes
cp_sensitivity_analysis_pop_rescaled <- cp_sensitivity_analysis_df %>%
    mutate(across(.cols = c(total_cases, peak_cases, total_vaccinated), 
                  .fns = ~ .x*target_pop)
           ) %>%  #rescale the population proportions to total sizes
    select(-c(npi_duration, total_vaccinated))





#' Towards the outbreak size isoclines

outbreak_size_cp_isocline_df <- cp_sensitivity_analysis_pop_rescaled %>% 
    filter(npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           total_cases <= 1000
           ) %>% 
    mutate(variant_emergence_day = as_factor(variant_emergence_day)) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity) %>% 
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

# Isoclines ----
# Outbreak size ====
outbreak_size_cp_isocline_sensitivity <- ggplot(outbreak_size_cp_isocline_df, 
                                 aes(x = vax_coverage, 
                                     y = min_speed, 
                                     color = variant_emergence_day
                                 )) + 
    geom_line(aes(linetype = factor(cross_protection_w)), 
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


print(outbreak_size_cp_isocline_sensitivity)

#Save the files 
ggsave(outbreak_size_cp_isocline_sensitivity,
       filename = 'outbreak_size_cp_isocline_sensitivity.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

ggsave(outbreak_size_cp_isocline_sensitivity,
       filename = 'outbreak_size_cp_isocline_sensitivity.eps',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')



# Peak prevalence ====

peak_prevalence_cp_isocline_df <- cp_sensitivity_analysis_pop_rescaled %>% 
    filter(npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           peak_cases <= 300) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity) %>% 
    mutate(variant_emergence_day = as_factor(variant_emergence_day)) %>% 
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

peak_prevalence_cp_isocline <- ggplot(peak_prevalence_cp_isocline_df, 
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


print(peak_prevalence_cp_isocline)


#Save the files 
ggsave(peak_prevalence_cp_isocline,
       filename = 'peak_prevalence_cp_isocline.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

ggsave(peak_prevalence_cp_isocline,
       filename = 'peak_prevalence_cp_isocline.eps',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')





