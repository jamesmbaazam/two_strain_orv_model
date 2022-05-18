#Model evaluation params ----
# Run the model for this number of days
max_time <- 365*5

# Show the model results at these time points
eval_times <- seq(0, max_time, 1)

# When the variant will emerge
# variant_emergence_day <- 60
# variant_emergence_day_vec <- seq(60, max_time, 1)

# Population parameters
target_pop <- 50E6 # target population size
Iw_index_cases <- 50 # Initial/Index wild type cases
Im_index_cases <- 50 # Initial/Index variant cases


# Population initial conditions ----
pop_inits_vax_escape <- c(
    S = 1 - Iw_index_cases / target_pop,
    Iw = Iw_index_cases / target_pop,
    Im = 0,
    Iwm = 0,
    Imw = 0,
    RwSm = 0,
    RmSw = 0,
    RwSmV = 0,
    RmSwV = 0,
    VIw = 0,
    VIm = 0,
    R = 0,
    V = 0,
    K = 0
)

# Control parameters ----

coverage_correction <- 0.9991



#' Event data frame for introducing the variant into the model dynamics
#' (check ?deSolve::event for more on the structure of the event_df below)
#'

# Variant emergence times ===
variant_emergence_times <- c(c(1, 61, 121, 151), max_time)

event_df <- data.frame(
    var = c("S", "Im"), # Compartments to change at a set time
    value = c(-Im_index_cases / target_pop, Im_index_cases / target_pop), # introduce 10 variant cases
    method = c("add", "replace")
) # operation on state variables






