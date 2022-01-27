
#Packages ----
library(doParallel)
library(doSNOW)
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model_functions.R')
source('./scripts/two_strain_model/simulation_functions.R')
source('./scripts/two_strain_model/vax_escape_analyses/vax_escape_sim_global_params.R')
source('./scripts/two_strain_model/vax_escape_analyses/vax_escape_sim_table_setup.R')



#The simulation table
#All scenarios
simulation_table <- vax_escape_mod_sim_table 

# Set up parallelization ----
# how many cores to use in the cluster? #
num_cores <- parallel::detectCores() - 1 

# set up a cluster called 'cl'
cl = makeSOCKcluster(num_cores)

# register the cluster
registerDoSNOW(cl)

## do some parallel computations with foreach
n_sims <- nrow(simulation_table)

sims_per_job <- ceiling(n_sims/num_cores)

num_of_jobs <- ceiling(n_sims/sims_per_job)

start_time <- Sys.time()

vax_escape_sim_output <- foreach(i = 1:num_of_jobs, 
                                          .combine = rbind,
                                          .packages = c('tidyverse', 'foreach', 'deSolve'), 
                                          .errorhandling = 'remove') %dopar% 
    {
        
        start_index <- 1 + (i - 1)*sims_per_job
        
        end_index <- min(start_index + sims_per_job - 1, n_sims)
        
        sim_subset <- simulation_table %>% 
            slice(start_index:end_index)
        
        
        subset_run <- run_batch_sims_vax_escape_model(sim_subset, get_summaries = TRUE) 
        
        return(subset_run)
    }

stopCluster(cl)


end_time <- Sys.time()

run_time <- end_time - start_time
print(run_time)

#save the simulation
saveRDS(object = vax_escape_sim_output, file = './model_output/vax_escape_analyses/vax_escape_perfect_efficacy_summaries.rds')


beepr::beep(3)