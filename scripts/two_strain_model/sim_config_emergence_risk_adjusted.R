#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')


# Uncontrolled epidemic simulation ----

# Controlled epidemic simulations set up ---- 

# Event times ====
# NPI ####
# Best case scenario ####
npi_start <- 1:max_time
npi_duration <- max_time - npi_start
npi_intensity <- seq(0, 1, 0.2) #Five levels of npi intensity (could correspond to the five stages in South Africa, for e.g)

# Vaccination ====
# Best case scenario ####
vax_start <- 1:max_time
campaign_duration <- max_time - vax_start
vax_cov <- seq(0, 1, 0.25)

# Variant emergence ====
variant_emergence_times <- seq(1, max_time, 5)



# Complete set of scenarios ====
# All combinations of vaccination coverage and NPI intensity  
control_scenarios <- expand.grid(npi_intensity = npi_intensity, vax_coverage = vax_cov) 

control_scenarios_rep <- control_scenarios %>% 
    slice(rep(1:n(), times = length(variant_emergence_times))) # repeat the control_scenarios df many times


# Repeat the variant emergence times for binding
variant_emergence_times_df <- data.frame(variant_emergence_day = rep(variant_emergence_times, each = nrow(control_scenarios))) 

control_and_emergence_scenarios <- cbind(control_scenarios_rep, variant_emergence_times_df)


# Full simulation table ====
simulation_table <- control_and_emergence_scenarios %>% 
    mutate(campaign_duration = campaign_duration,
           npi_duration = npi_duration
           )


# Run simulations ====
# No control ----
no_control_parms_df <- data.frame(npi_intensity = 0, vax_coverage = 0, 
                                  vax_start = 1, campaign_duration = 0, 
                                  npi_start = 1, npi_duration = 0
                                  )

no_control_epidemic <- variant_emergence_times %>% 
    purrr::map_df(function(x){simulate_model(pop_inits = pop_inits, 
                                      dynamics_parms = dynamics_params,
                                      control_parms = cbind(no_control_parms_df, data.frame(variant_emergence_day = x)),
                                      max_time = max_time, 
                                      dt = eval_times,
                                      events_table = event_df,
                                      return_dynamics = FALSE,
                                      browse = FALSE)}
                  )

#Rescale the population proportions to total sizes
no_control_epidemic_mod <- no_control_epidemic %>%
    mutate(total_cases = total_cases*target_pop,
           peak_cases = peak_cases*target_pop
           ) 


#Plot the final size
no_control_epidemic_final_size_plot <- ggplot(data = no_control_epidemic_mod) + 
    geom_line(aes(x = variant_emergence_day, 
                  y = total_cases
                  ),
              size = 1.2
              ) +
    scale_y_continuous(labels = comma) +
    labs(x = 'Variant emergence day', 
         y = 'Final size'
         ) +
    theme_minimal(base_size = 12)

plot(no_control_epidemic_final_size_plot)


#Plot the peak cases
no_control_epidemic_peak_cases_plot <- ggplot(data = no_control_epidemic_mod) + 
    geom_line(aes(x = variant_emergence_day, 
                  y = peak_cases
    ),
    size = 1.2
    ) +
    scale_y_continuous(labels = comma) +
    labs(x = 'Variant emergence day', 
         y = 'Peak cases'
    ) +
    theme_minimal(base_size = 12)

plot(no_control_epidemic_peak_cases_plot)

#Plot the two size by size
plot(no_control_epidemic_final_size_plot | no_control_epidemic_peak_cases_plot)

# Controlled epidemic
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

