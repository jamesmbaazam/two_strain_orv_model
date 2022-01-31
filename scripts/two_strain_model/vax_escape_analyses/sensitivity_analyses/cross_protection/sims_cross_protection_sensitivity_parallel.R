
# Packages ----
library(doParallel)
library(doSNOW)
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)
library(beepr)

# Helper scripts ----
source("./scripts/two_strain_model/simulation_functions.R")
source("./scripts/two_strain_model/two_strain_model_functions.R")

#Parameter tables ----

# Global parameters for population and simulation initialization
source("./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/global_params.R")

#Cross protection sensitivity params (specific to this simulation)
source("./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/cross_protection/dynamics_params_cross_protection_sensitivity.R")

#Global intervention params
source("./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/intervention_params_table_setup_global.R")

#' Towards the final simulation table
intervention_params_table <- simulation_table

dynamics_params_table_cp <- dynamics_params_cp_sensitivity

#'Bind the original intervention params df to itself as many times as the number
#'of rows in the original dynamics params df. We want to bind the two dfs together to
#'represent the various scenarios in the original dynamics df
intervention_params_expanded <- intervention_params_table %>%
    slice(rep(row_number(), times = nrow(dynamics_params_table_cp)))


#'Expand the original dynamics df so that each row is as many as the number 
#'of rows in the intervention controls df. We want to bind the two resulting 
#'data frames
cp_dynamics_parameters_expanded <- dynamics_params_table_cp %>% 
    slice(rep(row_number(), each = nrow(intervention_params_table)))

#' Attach the new dynamics params df to the new interventions df to form the 
#' simulation table specific for this simulation. Represents various scenarios
#' of vaccine efficacy loss against the variant
cp_simulation_table <- bind_cols(intervention_params_expanded, cp_dynamics_parameters_expanded) %>% 
    relocate(variant_emergence_day, .before = vax_rate)


# Set up parallelization ----
# how many cores to use in the cluster? #
num_cores <- parallel::detectCores() - 1

# set up a cluster called 'cl'
cl <- makeSOCKcluster(num_cores)

# register the cluster
registerDoSNOW(cl)

## do some parallel computations with foreach
n_sims <- nrow(cp_simulation_table)

sims_per_job <- ceiling(n_sims / num_cores)

num_of_jobs <- ceiling(n_sims / sims_per_job)

# pb <- txtProgressBar(min = 1, max = 1533000, style = 3)

# progress <- function(n){setTxtProgressBar(pb, n)}

# opts <- list(progress = progress)


start_time <- Sys.time()

cp_sensitivity_summaries <- foreach(
  i = 1:num_of_jobs,
  .combine = rbind,
  .packages = c("tidyverse", "foreach", "deSolve"),
  .errorhandling = "remove"
) %dopar% {
  start_index <- 1 + (i - 1) * sims_per_job

  end_index <- min(start_index + sims_per_job - 1, n_sims)

  sim_subset <- cp_simulation_table %>%
    slice(start_index:end_index)

  # run the model over the subset of rows
  subset_run <- run_batch_sims_vax_escape_model(sim_table = sim_subset, get_summaries = TRUE)

  return(subset_run)
}

## Shut down the timer and cluster
# close(pb)
stopCluster(cl)


end_time <- Sys.time()

run_time <- end_time - start_time
print(run_time)


# save the simulation
saveRDS(object = cp_sensitivity_summaries, file = "./model_output/sensitivity_analyses/cross_protection/cp_05_to_085_sensitivity_analysis_summaries.rds")


beepr::beep(3)
