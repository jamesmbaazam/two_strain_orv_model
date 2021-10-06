# Run the model for this number of days
max_time <- 365

# Evaluate the model at these time points
eval_times <- seq(0, max_time, 1) # Simulate a 2-year epidemic

# When the variant will emerge
#variant_emergence_day <- 60
#variant_emergence_day_vec <- seq(60, max_time, 1) 

#Population parameters
target_pop <- 1E6 #target population size
Iw_index <- 10 #Initial/Index wild type cases


# Population initial conditions
pop_inits <- c(S = 1 - Iw_index/target_pop, 
           Iw = Iw_index/target_pop, 
           Im = 0, 
           Iwm = 0,
           Imw = 0,
           RwSm = 0, 
           RmSw = 0,
           R = 0,
           V = 0,
           K = 0
           )


# ===============================
# Parameters for dynamics
# ===============================
dynamics_params <- data.frame(beta_w = 1.5/(7*target_pop), 
                        beta_m = 2/(7*target_pop),
                        gamma_w = 1/7,
                        gamma_m = 1/7,
                        sigma_w = 1,
                        sigma_m = 1
                        ) 


# ===============================
# Parameters for control dynamics
# ===============================
coverage_correction <- 0.999099












