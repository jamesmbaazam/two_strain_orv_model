#Running the model
library(deSolve)
library(tidyverse)
#library(patchwork)

# Source the helper scripts
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')

#NPI and vaccination start times
npi_implementation_day <- 30
campaign_start <- 180



#NPI and vaccination levels
npi_intensity <- seq(0, 1, 0.1)
vax_cov <- seq(0, 1, 0.1)


# Combine the params (Fixed campaign )
control_threshold_scenarios <- expand.grid(phi = npi_intensity, vax_coverage = vax_cov) 


# Complete set of scenarios ====
sim_table <- control_threshold_scenarios %>% 
    mutate(variant_emergence_day = rep(variant_emergence_day_vec, 
                                       each = nrow(control_threshold_scenarios)
                                       )
           )



# simulation_table <- scenario_params %>% 
#     mutate(variant_emergence_day = variant_emergence_day, 
#            npi_implementation_day = npi_implementation_day,
#            npi_duration = npi_duration,
#            campaign_start = campaign_start,
#            campaign_duration = campaign_duration
#            ) %>% 
#     relocate(variant_emergence_day, .before = 'phi')



#' Simulations ====

#' Toy simulation #### 
toy_simulation_table <- sim_table %>% 
    filter(npi_implementation_day == 46, campaign_start == 74) 

#' Event data frame for introducing mutant strain into model dynamics 
#' (check ?deSolve::event for more on the structure of the event_df below)

event_df <- data.frame(var = c('S', 'Im'), #Compartments to change at a set time
                       value = c(-0.01, 0.01), #index number of variant cases to introduce
                       method = c('add', 'replace') #operation on state variables
)


# Run the model with the toy simulation
toy_simulation_results <- toy_simulation_table %>% 
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


#Select the columns for the plot
toy_sim_hm_df <- toy_simulation_results %>% select(phi, vax_coverage, total_cases)

toy_sim_hm_plot <- ggplot(data = toy_sim_hm_df) + 
    geom_tile(aes(x = vax_coverage, 
                  y = phi, 
                  fill = total_cases
    ), 
    #  color = 'black',
    #  size = 0.01,
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
    theme_minimal()


print(toy_sim_hm_plot)



#' Full simulation #### 
full_simulation_results <- sim_table %>% 
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
full_sim_hm_plot <- ggplot(data = full_simulation_results) + 
    geom_tile(aes(x = vax_coverage, 
                  y = phi, 
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
    facet_wrap(variant_emergence_day ~ npi_implementation_day + campaign_start) + 
    theme_minimal()

plot(full_sim_hm_plot)

