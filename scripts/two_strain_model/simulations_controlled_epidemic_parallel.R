
#Packages ----
library(doParallel)
#library(microbenchmark)
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')

#function to run simulations
run_sim_all <- function(sim_table){
    res <- sim_table %>% 
        rowwise() %>% 
        do({with(., 
                 simulate_model(pop_inits = pop_inits, 
                                dynamics_parms = dynamics_params,
                                control_parms = .,
                                max_time = max_time, 
                                dt = eval_times,
                                events_table = event_df,
                                browse = F
                 )
        )
        }) %>% 
        ungroup() %>% 
        as_tibble()
    return(res)
}


# Set up parallelization ----
# how many cores to use in the cluster? #
num_cores <- detectCores() -1 

# set up a cluster called 'cl'
cl = makeCluster(num_cores)

# register the cluster
registerDoParallel(cl)

## do some parallel computations with foreach
n_sims <- nrow(simulation_table)

sims_per_job <- ceiling(n_sims/num_cores)

num_of_jobs <- ceiling(n_sims/sims_per_job)


orv_par_sim_output <- foreach(i = 1:num_of_jobs, 
                  .combine = rbind,
                  .packages = c('tidyverse', 'foreach', 'deSolve')) %dopar% 
    {
        
    start_index <- 1 + (i - 1)*sims_per_job
    
    end_index <- min(start_index + sims_per_job - 1, n_sims)
    
    sim_subset <- simulation_table %>% 
        slice(start_index:end_index)
        
        
    subset_run <- run_sim_all(sim_subset) 
        
    return(subset_run)
}

## Shut down the cluster 
stopCluster(cl)


#save the simulation
saveRDS(object = orv_par_sim_output, file = './model_output/simulation_controlled_epidemic_parallel_run.rds')