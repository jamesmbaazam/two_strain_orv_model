# Run the model for this number of days
max_time <- 365*2

# Evaluate the model at these time points
eval_times <- seq(0, max_time, 1) # Simulate a 2-year epidemic

# When the variant will emerge
#variant_emergence_day <- 60
#variant_emergence_day_vec <- seq(60, max_time, 1) 

#Population parameters
target_pop <- 50E6 #target population size
Iw_index <- 50 #Initial/Index wild type cases


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
R0_w <- 2
R0_m <- 2.5
IP_w <- 7 #Infectious period (wild type)
IP_m <- 7 #Infectious period (variant)
RP_w <- 14 #recovery period (wild type)
RP_m <- 14 #recovery period (variant)

dynamics_params <- data.frame(beta_w = 0.143, 
                        beta_m = 0.193,
                        gamma_w = 1/RP_w, 
                        gamma_m = 1/RP_m, 
                        sigma_w = 1, #cross-protection provided by wildtype
                        sigma_m = 1 #cross-protection provided by variant
                        ) 


# ===============================
# Parameters for control dynamics
# ===============================
coverage_correction <- 0.999099












