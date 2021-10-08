# Run the model for this number of days
max_time <- 365

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
R0_w <- 2.0 #R0 of the wild type
R0_m <- 3.0 #R0 of the variant
IP_w <- 7 #Infectious period (wild type)
IP_m <- 7 #Infectious period (variant)
RP_w <- 14 #recovery period (wild type)
RP_m <- 14 #recovery period (variant)

dynamics_params <- data.frame(beta_w = R0_w/RP_w, 
                        beta_m = R0_m/RP_m,
                        gamma_w = 1/RP_w, 
                        gamma_m = 1/RP_m, 
                        sigma_w = 1, #cross-protection provided by wildtype
                        sigma_m = 1 #cross-protection provided by variant
                        ) 



# ===============================
# Parameters for control dynamics
# ===============================
coverage_correction <- 0.999099



#' Event data frame for introducing the variant into the model dynamics 
#' (check ?deSolve::event for more on the structure of the event_df below)
#' 

# Variant emergence times ===
variant_emergence_times <- seq(1, max_time, 1)

event_df <- data.frame(var = c('S', 'Im'), #Compartments to change at a set time
                       value = c(-50/target_pop, 50/target_pop), #introduce 10 variant cases
                       method = c('add', 'replace')
                       ) #operation on state variables









