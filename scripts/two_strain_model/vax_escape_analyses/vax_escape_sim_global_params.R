# Run the model for this number of days
max_time <- 365

# Show the model results at these time points
eval_times <- seq(0, max_time, 1)

# When the variant will emerge
# variant_emergence_day <- 60
# variant_emergence_day_vec <- seq(60, max_time, 1)

# Population parameters
target_pop <- 50E6 # target population size
Iw_index_cases <- 50 # Initial/Index wild type cases
Im_index_cases <- 50 # Initial/Index variant cases


# Population initial conditions
pop_inits <- c(
  S = 1 - Iw_index_cases / target_pop,
  Iw = Iw_index_cases / target_pop,
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
R0_w <- 2.0 # R0 (wild type)
R0_m <- R0_w * 1.3 # R0 (variant) 30% more infectious
IP_w <- 14 # infectious period (wild type)
IP_m <- 14 # infectious period (variant)

dynamics_params <- data.frame(
  beta_w = R0_w / IP_w,
  beta_m = R0_m / IP_m,
  gamma_w = 1 / IP_w,
  gamma_m = 1 / IP_m,
  sigma_w = 1, # cross-protection provided by wild-type
  sigma_m = 1, # cross-protection provided by variant
  vax_efficacy_w = 1, #perfect vaccine efficacy against wild-type
  vax_efficacy_m = 1 #perfect vaccine efficacy against variant
)



# ===============================
# Parameters for control dynamics
# ===============================
coverage_correction <- 0.9991



#' Event data frame for introducing the variant into the model dynamics
#' (check ?deSolve::event for more on the structure of the event_df below)
#'

# Variant emergence times ===
variant_emergence_times <- c(seq(1, 61, 121, 151), max_time)

event_df <- data.frame(
  var = c("S", "Im"), # Compartments to change at a set time
  value = c(-Im_index_cases / target_pop, Im_index_cases / target_pop), # introduce 10 variant cases
  method = c("add", "replace")
) # operation on state variables
