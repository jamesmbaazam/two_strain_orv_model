
#Packages ----
library(doParallel)
library(doSNOW)
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')
source('./scripts/two_strain_model/sim_config_uncontrolled_epidemic.R')
source('./scripts/two_strain_model/simulation_functions.R')


#The simulation table
#Select the control case studies
control_sim_params <- orv_npi_control_config_table %>% 
    filter(npi_intensity %in% c(0, 0.1), 
           vax_coverage == 0.80, 
           vax_speed %in% c(1, 3.75), 
           variant_emergence_day %in% c(31, max_time)
           ) %>% 
    mutate(is_control_scenario = TRUE)


#Select the case studies for no control
no_control_sim_parms <- no_control_parms_df %>% 
    rbind(no_control_parms_df) %>% 
    mutate(variant_emergence_day = c(31, 365),
           is_control_scenario = FALSE
           )

#The final simulation tables for the controlled versus uncontrolled dynamics
all_case_studies_sim_table <- rbind(control_sim_params, no_control_sim_parms)


# Set up parallelization ----
# how many cores to use in the cluster? #
num_cores <- parallel::detectCores() - 1 

# set up a cluster called 'cl'
cl = makeSOCKcluster(num_cores)

# register the cluster
registerDoSNOW(cl)

## do some parallel computations with foreach
n_sims <- nrow(all_case_studies_sim_table)

sims_per_job <- ceiling(n_sims/num_cores)

num_of_jobs <- ceiling(n_sims/sims_per_job)

#pb <- txtProgressBar(min = 1, max = 1533000, style = 3)

#progress <- function(n){setTxtProgressBar(pb, n)}

#opts <- list(progress = progress)


start_time <- Sys.time()

dynamics_case_studies <- foreach(i = 1:num_of_jobs, 
                                          .combine = rbind,
                                          .packages = c('tidyverse', 'foreach', 'deSolve'), 
                                          .errorhandling = 'remove') %dopar% 
    {
        
        start_index <- 1 + (i - 1)*sims_per_job
        
        end_index <- min(start_index + sims_per_job - 1, n_sims)
        
        sim_subset <- all_case_studies_sim_table %>% 
            slice(start_index:end_index)
        
        
        subset_run <- run_sim_all(sim_subset, get_summaries = FALSE) 
        
        return(subset_run)
    }

## Shut down the timer and cluster 
#close(pb)
stopCluster(cl)


end_time <- Sys.time()

run_time <- end_time - start_time
print(run_time)

#save the simulation
saveRDS(object = dynamics_case_studies, file = './model_output/dynamics_case_studies_controlled_vs_uncontrolled.rds')


beepr::beep(3)
