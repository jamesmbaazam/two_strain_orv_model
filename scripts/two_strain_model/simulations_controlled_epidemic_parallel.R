
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


#function to run simulations
run_sim_all <- function(sim_table){
    res <- sim_table %>% 
        rowwise() %>% 
        do({with(., 
                 simulate_raw_dynamics(pop_inits = pop_inits, 
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


#The simulation table
simulation_table <- orv_npi_control_config_table %>% 
    filter(npi_intensity %in% c(0, 0.5)) #npi levels of 0%, 10%, and 50%


# Set up parallelization ----
# how many cores to use in the cluster? #
num_cores <- parallel::detectCores() -1 

# set up a cluster called 'cl'
cl = makeSOCKcluster(num_cores)

# register the cluster
registerDoSNOW(cl)

## do some parallel computations with foreach
n_sims <- nrow(orv_npi_control_config_table)

sims_per_job <- ceiling(n_sims/num_cores)

num_of_jobs <- ceiling(n_sims/sims_per_job)

pb <- txtProgressBar(min = 1, max = 1533000, style = 3)

progress <- function(n){setTxtProgressBar(pb, n)}

opts <- list(progress = progress)


start_time <- Sys.time()

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

## Shut down the timer and cluster 
close(pb)
stopCluster(cl)


end_time <- Sys.time()

run_time <- end_time - start_time
print(run_time)

#save the simulation
saveRDS(object = orv_par_sim_output, file = './model_output/controlled_epidemic_dynamics_parallel_run.rds')


beepr::beep(3)
