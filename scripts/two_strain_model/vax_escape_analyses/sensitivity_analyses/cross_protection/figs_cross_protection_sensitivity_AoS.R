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



cp_sensitivity_analysis_aos <- cp_sensitivity_analysis_df %>% 
    mutate(outbreak_size_sufficient = as_factor(ifelse(total_cases <= 0.1, 'yes', 'no'))) %>% 
    group_by(variant_emergence_day, npi_intensity, cross_protection_m, outbreak_size_sufficient) 
