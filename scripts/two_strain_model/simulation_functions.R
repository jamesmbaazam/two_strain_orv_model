#Packages ----
library(doParallel)
library(deSolve)
library(tidyverse)
library(beepr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')


#' Function to simulate the parameter combination by row 
#'
#' @param sim_table The data.frame of simulation parameter combinations
#' @param output_type The type of output to return (TRUE = the raw dynamics; FALSE = the summaries)
#'
#' @return a data.frame of simulation results, one per row 
#' @export
#'
#' @examples
run_sim_rowwise <- function(sim_table, return_dynamics){
    output <- sim_table %>%
        rowwise() %>%
        do({with(.,
                 simulate_model(pop_inits = pop_inits,
                                dynamics_parms = dynamics_params,
                                control_parms = .,
                                max_time = max_time,
                                dt = eval_times,
                                events_table = event_df,
                                return_dynamics = return_dynamics,
                                browse = F
                 )
        )
        }) %>%
        ungroup() %>%
        as_tibble()
    return(output)
}


#' Function to run the complete simulation either in parallel or sequentially
#' 
#' @param return_dynamics #If true, the function returns the raw dynamics; 
#' if false, it returns the summaries
#' @param sim_table #The simulation parameter combinations df
#' @param run_type #One of c('parallel', 'sequential'). The parallel option runs the model
#' on one less core on your machine
#'
#' @return
#' @export
#'
#' @examples
run_sim <- function(sim_table, run_type, return_raw_dynamics) {
  if (run_type == "parallel") {
    
      # Set up parallelization 
      
    # how many cores to use in the cluster? #
    num_cores <- detectCores() - 1

    # set up a cluster called 'cl'
    cl <- makeCluster(num_cores)

    # register the cluster
    registerDoParallel(cl)

    ## do some parallel computations with foreach
    n_sims <- nrow(sim_table)

    sims_per_job <- ceiling(n_sims / num_cores)

    num_of_jobs <- ceiling(n_sims / sims_per_job)

    start_time <- Sys.time()

    sim_output <- foreach(
      i = 1:num_of_jobs,
      .combine = rbind,
      .packages = c("tidyverse", "foreach", "deSolve")
      ) %dopar% {
          
          #The lower index for subsetting the simulation table for the parallel run
          start_index <- 1 + (i - 1) * sims_per_job
          
          #The upper index
          end_index <- min(start_index + sims_per_job - 1, n_sims)
          
          # Split up the simulation table for a parallel run
          
         sim_subset <- sim_table %>%
              slice(start_index:end_index) 
         
         subset_run <- run_sim_rowwise(sim_subset, return_dynamics = return_raw_dynamics) 
         return(subset_run)
          } 
    
    ## Shut down the cluster
    stopCluster(cl)
    
    end_time <- Sys.time()
    
    parallel_run_time <- start_time - end_time
    
    print(paste0("Parallel run time = ", parallel_run_time))
    } else if (run_type == "sequential") {
    sim_output <- run_sim_rowwise(sim_table, return_dynamics = return_raw_dynamics)
  } else {
    stop("Run type must be either parallel or sequential")
  }
  return(sim_output)
}