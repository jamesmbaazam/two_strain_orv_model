#Packages ----
library(dplyr)
library(ggplot2)
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
cp_sensitivity_analysis_all_results <- readRDS('./model_output/sensitivity_analyses/cross_protection/cp_all_scenarios_sensitivity_analysis_summaries.rds')

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
# cp_sensitivity_analysis_pop_rescaled <- cp_sensitivity_analysis_df %>%
#     mutate(across(.cols = c(total_cases, peak_prevalence), 
#                   .fns = ~ .x*target_pop
#                   )
#            ) 


#' Determine the subset of scenarios that meet the threshold outbreak size
#' and the minimum speed required 
total_cases_thresholds <- c(0.01, 0.1, 0.15, 0.2)


cases_isocline_plot_list <- list()

for (cases_threhold in seq_along(total_cases_thresholds)) {
outbreak_size_cp_isocline_df <- cp_sensitivity_analysis_df %>% 
    filter(cross_protection_w %in% c(0.5, 1), 
           variant_emergence_day %in% c(1, 61, 121, 151, max_time), 
           npi_intensity %in% c(0.00, 0.10, 0.20, 0.30)
           ) %>% 
    filter(total_cases <= total_cases_thresholds[[cases_threhold]]) %>%  #keep scenarios with this level of total cases
    mutate(variant_emergence_day = as_factor(variant_emergence_day), 
           cross_protection_w = as_factor(cross_protection_w)
           ) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity, cross_protection_m, cross_protection_w) %>% 
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
              size = 0.5, 
              show.legend = TRUE
              ) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.10, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(title = 'Sensitivity to cross protection assumptions',
         subtitle = paste('Strategies with cumulative cases up to', 
                          total_cases_thresholds[[cases_threhold]]*100, 
                          "% of total population"
                          ),
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day',
        linetype = 'Cross protection'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'right') +
    theme_bw(base_size = 14) 


cases_isocline_plot_list[[cases_threhold]] <- outbreak_size_cp_isocline_sensitivity

file_name <-  paste0('./figures/sensitivity_analyses/cross_protection/cp_sensitivity_isocline_outbreak_size', total_cases_thresholds[cases_threhold]*100,'.png')
}

#Save the files to a pdf
pdf(width = 12, file = "./figures/sensitivity_analyses/cross_protection/cp_isoclines_outbreak_size.pdf")

for (cases_threhold in seq_along(total_cases_thresholds)) {
    print(cases_isocline_plot_list[[cases_threhold]])
}

dev.off()

# ggsave(outbreak_size_cp_isocline_sensitivity,
#        filename = 'sensitivity_analyses/outbreak_size_isocline_cp_sensitivity_analysis.png',
#        path = git_plot_path,
#        width = 23.76,
#        height = 17.86,
#        units = 'cm')

# ggsave(outbreak_size_cp_isocline_sensitivity,
#        filename = 'outbreak_size_isocline_cp_sensitivity_analysis.eps',
#        path = thesis_plot_path,
#        width = 23.76,
#        height = 17.86,
#        units = 'cm')



# Peak prevalence ====
peak_prevalence_thresholds <- total_cases_thresholds/100

peak_prevalence_isocline_plot_list <- list()

for (peak_prevalence_threhold in seq_along(peak_prevalence_thresholds)) {
peak_prevalence_cp_isocline_df <- cp_sensitivity_analysis_df %>% 
    filter(cross_protection_w %in% c(0.5, 1), 
           variant_emergence_day %in% c(1, 61, 121, 151, max_time),
           npi_intensity %in% c(0.0, 0.1, 0.2, 0.3)
           ) %>% 
    filter(peak_prevalence <= peak_prevalence_thresholds[[peak_prevalence_threhold]]) %>%  #keep scenarios with this level of total cases
    mutate(variant_emergence_day = as_factor(variant_emergence_day), 
           cross_protection_w = as_factor(cross_protection_w)
    ) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity, cross_protection_m, cross_protection_w) %>% 
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

peak_prevalence_cp_isocline <- ggplot(peak_prevalence_cp_isocline_df, 
                                   aes(x = vax_coverage, 
                                       y = min_speed, 
                                       color = variant_emergence_day
                                   )) + 
    geom_line(aes(linetype = cross_protection_w), size = 0.5, show.legend = TRUE) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.10, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(title = 'Sensitivity to cross protection assumptions', 
         subtitle = paste('Strategies with peak prevalence up to', 
                          peak_prevalence_thresholds[[peak_prevalence_threhold]]*100, 
                          '% of total population'
                          ), 
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day',
        linetype = 'Cross protection'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'right') 


peak_prevalence_isocline_plot_list[[peak_prevalence_threhold]] <- peak_prevalence_cp_isocline

file_name <-  paste0('./figures/sensitivity_analyses/cross_protection/cp_sensitivity_peak_prevalence_isocline_', 
                     peak_prevalence_thresholds[[peak_prevalence_threhold]]*100,
                     '.png'
                     )
}

#Save the files to a pdf
pdf(width = 12, file = "./figures/sensitivity_analyses/cross_protection/cp_isoclines_peak_prevalence.pdf")

for (peak_prevalence_threhold in seq_along(peak_prevalence_thresholds)) {
    print(peak_prevalence_isocline_plot_list[[peak_prevalence_threhold]])
}

dev.off()

#Save the files 
# ggsave(peak_prevalence_cp_isocline,
#        filename = '/sensitivity_analyses/peak_prevalence_cp_isocline.png',
#        path = git_plot_path,
#        width = 23.76,
#        height = 17.86,
#        units = 'cm')

# ggsave(peak_prevalence_cp_isocline,
#        filename = 'peak_prevalence_cp_isocline.eps',
#        path = thesis_plot_path,
#        width = 23.76,
#        height = 17.86,
#        units = 'cm')





