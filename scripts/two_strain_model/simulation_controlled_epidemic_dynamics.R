
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
source('./scripts/two_strain_model/simulation_functions.R')


#The simulation table
#Select scenarios
dynamics_simulation_table <- orv_npi_control_config_table %>% 
    filter(npi_intensity %in% c(0, 0.1), 
           vax_coverage == 0.80, 
           vax_speed %in% c(1, 3.75), 
           variant_emergence_day %in% c(seq(1, max_time, by = 30), max_time)
           )

# orv_npi0_20_simulation_table <- orv_npi_control_config_table %>% 
#     filter(npi_intensity %in% c(0, 0.2)) #npi of 50%

# Set up parallelization ----
# how many cores to use in the cluster? #
num_cores <- parallel::detectCores() - 1 

# set up a cluster called 'cl'
cl = makeSOCKcluster(num_cores)

# register the cluster
registerDoSNOW(cl)

## do some parallel computations with foreach
n_sims <- nrow(dynamics_simulation_table)

sims_per_job <- ceiling(n_sims/num_cores)

num_of_jobs <- ceiling(n_sims/sims_per_job)

#pb <- txtProgressBar(min = 1, max = 1533000, style = 3)

#progress <- function(n){setTxtProgressBar(pb, n)}

#opts <- list(progress = progress)


start_time <- Sys.time()

control_dynamics_case_studies <- foreach(i = 1:num_of_jobs, 
                                          .combine = rbind,
                                          .packages = c('tidyverse', 'foreach', 'deSolve'), 
                                          .errorhandling = 'remove') %dopar% 
    {
        
        start_index <- 1 + (i - 1)*sims_per_job
        
        end_index <- min(start_index + sims_per_job - 1, n_sims)
        
        sim_subset <- dynamics_simulation_table %>% 
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
saveRDS(object = control_dynamics_case_studies, file = './model_output/control_dynamics_case_studies.rds')


beepr::beep(3)
