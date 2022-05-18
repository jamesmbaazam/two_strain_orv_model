.args <- if(interactive()){
    c('./data/inputs/config_global_params.RData', 
      './data/inputs/two_strain_model_functions.RData', 
      './data/inputs/simulation_functions.RData', 
      "./data/inputs/config_intervention_params.RData",
    './data/outputs/all_scenarios_summaries.rds'
    )
}else{
    commandArgs(trailingOnly = TRUE)
}


#Error handling
options(error = function(){
    traceback(2, max.lines = 3);
    if(!interactive()) quit("no", status = 1, runLast = FALSE)
    })

#Packages ----
library(doParallel, quietly = TRUE)
library(doSNOW, quietly = TRUE)
library(deSolve, quietly = TRUE)
library(tidyverse, quietly = TRUE)

# Load helper functions and param definitions ----
load(.args[[1]])
load(.args[[2]])
load(.args[[3]])
load(.args[[4]])

#function to run simulations
run_sim_all <- function(sim_table, get_summaries){
    res <- sim_table %>% 
        rowwise() %>% 
        do({with(.,
                 simulate_raw_dynamics(pop_inits = pop_inits, 
                                dynamics_parms = dynamics_params,
                                control_parms = .,
                                max_time = max_time, 
                                dt = eval_times,
                                events_table = event_df,
                                get_summaries = get_summaries,
                                browse = FALSE
                                )
           )
        }) %>% 
        ungroup() %>% 
        as_tibble()
    return(res)
}


#The simulation table
#All scenarios
orv_npi_all_scenarios_simulation_table <- orv_npi_control_config_table 

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
n_sims <- nrow(orv_npi_all_scenarios_simulation_table)

sims_per_job <- ceiling(n_sims/num_cores)

num_of_jobs <- ceiling(n_sims/sims_per_job)

#pb <- txtProgressBar(min = 1, max = 1533000, style = 3)

#progress <- function(n){setTxtProgressBar(pb, n)}

#opts <- list(progress = progress)


start_time <- Sys.time()

sim_results <- foreach(i = 1:num_of_jobs, 
                  .combine = rbind,
                  .packages = c('dplyr', 'foreach', 'deSolve'), 
                  .errorhandling = 'remove') %dopar% 
    {
        
    start_index <- 1 + (i - 1)*sims_per_job
    
    end_index <- min(start_index + sims_per_job - 1, n_sims)
    
    sim_subset <- orv_npi_all_scenarios_simulation_table %>% 
        slice(start_index:end_index)
        
        
    subset_run <- run_sim_all(sim_subset, get_summaries = TRUE) 
        
    return(subset_run)
}

## Shut down the timer and cluster 
#close(pb)
stopCluster(cl)


end_time <- Sys.time()

run_time <- end_time - start_time
print(run_time)

#save the simulation
saveRDS(object = sim_results, file = tail(.args, 1))

