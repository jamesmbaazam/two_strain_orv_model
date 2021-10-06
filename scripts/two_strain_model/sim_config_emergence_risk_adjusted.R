#Packages ----
library(deSolve)
library(tidyverse)
#library(patchwork)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')


# Simulation set up ---- 

# Event times ====
# NPI ####
# Best case scenario ####
npi_start <- 1:max_time
npi_duration <- max_time - npi_start
npi_intensity <- seq(0, 1, 0.2) #Five levels of npi intensity (could correspond to the five stages in South Africa, for e.g)

# Vaccination ====
# Best case scenario ####
vax_start <- 1
campaign_duration <- max_time - vax_start
vax_cov <- seq(0, 1, 0.25)

# Variant emergence ====
variant_emergence_times <- seq(1, max_time, 1)



# Complete set of scenarios ====
# All combinations of vaccination coverage and NPI intensity  
control_scenarios <- expand.grid(npi_intensity = npi_intensity, vax_coverage = vax_cov) 

control_scenarios_rep <- control_scenarios %>% 
    slice(rep(1:n(), times = length(variant_emergence_times))) # repeat the control_scenarios df many times


# Repeat the variant emergence times for binding
variant_emergence_times_rep <- data.frame(variant_emergence_day = rep(variant_emergence_times, each = nrow(control_scenarios))) 

# Full simulation table ====
simulation_table <- cbind(control_scenarios_rep, variant_emergence_times_rep) %>% 
    mutate(vax_start = vax_start, 
           campaign_duration = campaign_duration,
           npi_start = npi_start,
           npi_duration = npi_duration
           )



# Simulations ----
# Events ====
#' Event data frame for introducing mutant strain into model dynamics 
#' (check ?deSolve::event for more on the structure of the event_df below)

event_df <- data.frame(var = c('S', 'Im'), #Compartments to change at a set time
                       value = c(-1/target_pop, 1/target_pop), #index number of variant cases to introduce
                       method = c('add', 'replace')
                       ) #operation on state variables


# Uncontrolled epidemic ====
no_control_parms_df <- data.frame(npi_intensity = 0, vax_coverage = 0, variant_emergence_day = 1, 
                            vax_start = 1, campaign_duration = 0, 
                            npi_start = 1, npi_duration = 0
                            )

no_control_epidemic <- simulate_model(pop_inits = pop_inits, 
                                      dynamics_parms = dynamics_params,
                                      control_parms = no_control_parms_df,
                                      max_time = max_time, 
                                      dt = eval_times,
                                      events_table = event_df,
                                      return_dynamics = TRUE,
                                      browse = FALSE
                                      )


no_control_epidemic_long <- no_control_epidemic %>%
    select(time, Iw, Im, Iwm, Imw) %>% 
    pivot_longer(cols = -time, names_to = 'state', values_to = 'value')


ggplot(data = no_control_epidemic_long, aes(x = time, y = value, color = state)) +
    geom_line(size = 1)




# Controlled epidemic ====
orv_full_simulation <- simulation_table %>% 
    filter(variant_emergence_day %in% seq(1, 150, 7)) %>% 
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


# Select the columns for the plot
orv_full_simulation_plot <- ggplot(data = orv_full_simulation) + 
    geom_tile(aes(x = vax_coverage, 
                  y = npi_intensity, 
                  fill = total_cases
    ), 
    stat = 'identity'
    ) +
    scale_x_continuous(labels = scales::percent_format(), 
                       expand = c(0,0)
    ) +
    scale_y_continuous(labels = scales::percent_format(), 
                       expand = c(0,0)
    ) +
    scale_fill_viridis_b(direction = -1) +
    labs(x = 'Vaccination coverage', 
         y = 'NPI intensity', 
         fill = 'Cumulative incidence'
    ) +
    facet_wrap('variant_emergence_day') + 
    theme_minimal()

plot(orv_full_simulation_plot)

