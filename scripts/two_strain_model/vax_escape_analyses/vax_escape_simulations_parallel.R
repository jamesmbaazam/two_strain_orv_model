
#Packages ----
library(doParallel)
library(doSNOW)
library(deSolve)
library(dplyr)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model_functions.R')
source('./scripts/two_strain_model/simulation_functions.R')
source('./scripts/two_strain_model/vax_escape_analyses/vax_escape_sim_global_params.R')
source('./scripts/two_strain_model/vax_escape_analyses/vax_escape_sim_table_setup.R')



#The simulation parameters
intervention_params_table <- vax_escape_mod_sim_table 

baseline_dynamics_params <- dynamics_params


#' Attach the one row dynamics params df to each row of the MxN intervention 
#' params df to form the simulation table specific for this simulation
baseline_analysis_simulation_table <- intervention_params_table %>%
    mutate(baseline_dynamics_params) 


# Set up parallelization ----
# how many cores to use in the cluster? #
num_cores <- parallel::detectCores() - 1 

# set up a cluster called 'cl'
cl = makeSOCKcluster(num_cores)

# register the cluster
registerDoSNOW(cl)

## do some parallel computations with foreach
n_sims <- nrow(baseline_analysis_simulation_table)

sims_per_job <- ceiling(n_sims/num_cores)

num_of_jobs <- ceiling(n_sims/sims_per_job)

start_time <- Sys.time()

baseline_analysis_summaries <- foreach(i = 1:num_of_jobs, 
                                          .combine = rbind,
                                          .packages = c('dplyr', 'foreach', 'deSolve'), 
                                          .errorhandling = 'remove') %dopar% 
    {
        
        start_index <- 1 + (i - 1)*sims_per_job
        
        end_index <- min(start_index + sims_per_job - 1, n_sims)
        
        sim_subset <- baseline_analysis_simulation_table %>% 
            slice(start_index:end_index)
        
        
        subset_run <- run_batch_sims_vax_escape_model(sim_subset, get_summaries = TRUE) 
        
        return(subset_run)
    }

stopCluster(cl)


end_time <- Sys.time()

run_time <- end_time - start_time
print(run_time)

#save the simulation
saveRDS(object = baseline_analysis_summaries, file = './model_output/sensitivity_analyses/baseline_analysis_vR0_30_percent/vR0_30_percent_baseline_analysis_summaries.rds')


beepr::beep(3)