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
thesis_plot_path <- 'C:/Users/James Azam/Dropbox/My Academic Repository/_Degrees/PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/'

#Change how large numbers are printed
options(scipen = 10000)


#Load the model output

# Original analysis
#' The baseline analysis assumes perfect vaccine efficacy against the wild-type 
#' and variant and perfect cross protection between the strains. The variant is
#' assumed to be 30% more transmissible than the wild-type (R0w = 1.3R0m)
baseline_analysis_results <- readRDS('./model_output/sensitivity_analyses/baseline_analysis_vR0_30_percent/vR0_30_percent_baseline_analysis_summaries.rds') %>% 
    rename(peak_prevalence = peak_cases)

#Cross protection sensitivity analysis results
cp_05_085_sensitivity_analysis_results <- readRDS('./model_output/sensitivity_analyses/cross_protection/cp_05_to_085_sensitivity_analysis_summaries.rds') 
cp_0_sensitivity_analysis_results <- readRDS('./model_output/sensitivity_analyses/cross_protection/cp_0_sensitivity_analysis_summaries.rds')

#combine the two model outputs
cp_sensitivity_analysis_all_results <- bind_rows(cp_0_sensitivity_analysis_results, cp_05_085_sensitivity_analysis_results, baseline_analysis_results) 

#Remove redundant columns
cp_sensitivity_analysis_df <- cp_sensitivity_analysis_all_results %>% 
    select(-c(vax_rate, 
              npi_duration, 
              starts_with('vax_efficacy'), 
            #  cross_protection_m,
              R0m,
              total_vaccinated
              )
           )

#Rescale the population proportions to total sizes
cp_sensitivity_analysis_pop_rescaled <- cp_sensitivity_analysis_df %>%
    mutate(across(.cols = c(total_cases, peak_prevalence), 
                  .fns = ~ .x*target_pop
                  )
           ) 

#' Determine the subset of scenarios that meet the threshold outbreak size
#' and the minimum speed required 
outbreak_size_cp_isocline_df <- cp_sensitivity_analysis_pop_rescaled %>% 
    filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time), 
           npi_intensity %in% c(0.00, 0.10, 0.20, 0.30), 
           total_cases <= 1000
           ) %>% #only keep scenarios with 1000 total cases or less
    mutate(variant_emergence_day = as_factor(variant_emergence_day), 
           cross_protection_w = as_factor(cross_protection_w)
           ) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity,  cross_protection_m, cross_protection_w) %>% 
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

# Isoclines ----
# Outbreak size ====
outbreak_size_cp_isocline_sensitivity <- ggplot(outbreak_size_cp_isocline_df, 
                                 aes(x = vax_coverage, 
                                     y = min_speed, 
                                     color = variant_emergence_day
                                 )) + 
    geom_line(aes(linetype = cross_protection_w), 
              size = 1, 
              show.legend = TRUE
              ) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.30, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(title = paste('Sensitivity to cross protection assumptions'),
         subtitle = 'Strategies with cumulative cases <= 1000',
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day',
        linetype = 'Cross protection'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'right') 


print(outbreak_size_cp_isocline_sensitivity)

#Save the files 
ggsave(outbreak_size_cp_isocline_sensitivity,
       filename = 'sensitivity_analyses/outbreak_size_isocline_cp_sensitivity_analysis.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

# ggsave(outbreak_size_cp_isocline_sensitivity,
#        filename = 'outbreak_size_isocline_cp_sensitivity_analysis.eps',
#        path = thesis_plot_path,
#        width = 23.76,
#        height = 17.86,
#        units = 'cm')


# Outbreak size isocline for extreme assumptions of no or perfect cross protection

outbreak_size_cp_isocline_extremes_df <- outbreak_size_cp_isocline_df %>% 
filter(cross_protection_w %in% c(0, 1))


outbreak_size_cp_isocline_sensitivity_extremes <- ggplot(outbreak_size_cp_isocline_extremes_df, 
                                                aes(x = vax_coverage, 
                                                    y = min_speed, 
                                                    color = variant_emergence_day
                                                )) + 
    geom_line(aes(linetype = cross_protection_w), 
              size = 1, 
              show.legend = TRUE
    ) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.30, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(title = paste('Sensitivity to cross protection assumptions'),
         subtitle = 'Strategies with cumulative cases <= 1000',
         x = 'Vaccination coverage', 
         y = 'Vaccination speed', 
         color = 'Variant emergence day',
         linetype = 'Cross protection'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'right') 


print(outbreak_size_cp_isocline_sensitivity_extremes)



# Peak prevalence ====

peak_prevalence_cp_isocline_df <- cp_sensitivity_analysis_pop_rescaled %>% 
    filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time),
           npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           peak_prevalence <= 300) %>% 
     mutate(variant_emergence_day = as_factor(variant_emergence_day), 
           vax_coverage = as.factor(vax_coverage), 
           npi_intensity = as.factor(npi_intensity),
           cross_protection_w = as_factor(cross_protection_w)
           ) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity, cross_protection_w) %>% 
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

peak_prevalence_cp_isocline <- ggplot(peak_prevalence_cp_isocline_df, 
                                   aes(x = vax_coverage, 
                                       y = min_speed, 
                                       color = variant_emergence_day
                                   )) + 
    geom_line(aes(linetype = cross_protection_w), size = 1, show.legend = TRUE) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.30, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(title = 'Sensitivity to cross protection assumptions', 
         subtitle = 'Strategies with cumulative cases <= 1000',  
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day',
        linetype = 'Cross protection'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(peak_prevalence_cp_isocline)


#Save the files 
ggsave(peak_prevalence_cp_isocline,
       filename = '/sensitivity_analyses/peak_prevalence_cp_isocline.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

# ggsave(peak_prevalence_cp_isocline,
#        filename = 'peak_prevalence_cp_isocline.eps',
#        path = thesis_plot_path,
#        width = 23.76,
#        height = 17.86,
#        units = 'cm')





