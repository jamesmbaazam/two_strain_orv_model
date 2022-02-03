# Packages ----
library(deSolve)
library(dplyr)
library(scales)
library(patchwork)
library(beepr)

# Helper scripts ----
source("./scripts/two_strain_model/simulation_functions.R")
source("./scripts/two_strain_model/two_strain_model_functions.R")

# Parameter tables ----

# Global parameters for population and simulation initialization
source("./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/global_params.R")

# Cross protection sensitivity params (specific to this simulation)
source("./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/cross_protection/dynamics_params_cross_protection_sensitivity.R")

# Global intervention params
source("./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/intervention_params_table_setup_global.R")

#' Towards the final simulation table
intervention_params_table <- simulation_table

dynamics_params_table_cp <- dynamics_params_cp_sensitivity

#' Bind the original intervention params df to itself as many times as the number
#' of rows in the original dynamics params df. We want to bind the two dfs together to
#' represent the various scenarios in the original dynamics df
intervention_params_expanded <- intervention_params_table %>%
  slice(rep(row_number(), times = nrow(dynamics_params_table_cp)))


#' Expand the original dynamics df so that each row is as many as the number
#' of rows in the intervention controls df. We want to bind the two resulting
#' data frames
cp_dynamics_parameters_expanded <- dynamics_params_table_cp %>%
  slice(rep(row_number(), each = nrow(intervention_params_table)))

#' Attach the new dynamics params df to the new interventions df to form the
#' simulation table specific for this simulation. Represents various scenarios
#' of vaccine efficacy loss against the variant
#' 
#' Full simulation table for all cross protection scenarios
cp_simulation_table <- bind_cols(intervention_params_expanded, cp_dynamics_parameters_expanded) %>%
  relocate(variant_emergence_day, .before = vax_rate)


# 
# cp_simulation_table_subset <- cp_simulation_table 
# # %>%
# #   filter(variant_emergence_day %in% c(1, max_time), sigma_w %in% c(0.5, 1))

# Controlled epidemic ----

#############################

# Simulate the cross protection sensitivity analysis
#############################
cp_sensitivity_analysis_summaries <- run_batch_sims_vax_escape_model(sim_table = cp_simulation_table, get_summaries = TRUE)


# save the results
saveRDS(object = cp_sensitivity_analysis_summaries, file = "./model_output/sensitivity_analyses/cross_protection/cp_sensitivity_analysis_d1_dmaxtime_summaries")

beep(sound = 3)
