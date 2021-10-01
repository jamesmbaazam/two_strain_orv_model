
#Running the model
library(deSolve)
library(tidyverse)
#library(patchwork)

# Source the helper scripts
source('./scripts/two_strain_model/two_strain_model.R')

# Run the model for this number of days
max_time <- 365*2

# Evaluate the model at these time points
eval_times <- seq(0, max_time, 1) # Simulate a 2-year epidemic

# When the variant will emerge
#variant_emergence_day <- 60
variant_emergence_day <- seq(60, max_time/2, 30) 

# Population initial conditions
pop_inits <- c(S = 0.99, 
           Iw = 0.01, 
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
dynamics_params <- data.frame(beta_w = 1.5/7, 
                        beta_m = 2/7,
                        gamma_w = 1/14,
                        gamma_m = 1/14,
                        sigma_w = 1,
                        sigma_m = 1
                        ) 


# ===============================
# Parameters for control dynamics
# ===============================
npi_implementation_lag14 <- variant_emergence_day - 14 #Implement NPIs 14 days before a new variant emerges
npi_implementation_lead14 <- variant_emergence_day + 14 #Implement NPIs 14 days before a new variant emerges
npi_duration <- 60
campaign_start_lag14 <- variant_emergence_day - 14 #vaccination campaign starts 30 days before more infectious variant emerges
campaign_start_lead14 <- variant_emergence_day + 14
campaign_duration <- 365
coverage_correction <- 0.999099



#NPI and vaccination
npi_intensity <- seq(0, 1, 0.1)
vax_cov <- seq(0, 1, 0.1)


# Combine the params (Fixed campaign )
control_threshold_scenarios <- expand.grid(phi = npi_intensity, 
                                 vax_coverage = vax_cov
                                 ) 

control_times_scenarios <- data.frame(npi_implementation_day = rep(c(npi_implementation_lag14, 
                                                            npi_implementation_lead14
                                                            ), 
                                                            each = nrow(control_threshold_scenarios)
                                                          ),
                                     campaign_start = rep(c(campaign_start_lead14, 
                                                            campaign_start_lag14
                                                            ),
                                                          each = nrow(control_threshold_scenarios)
                                                          )
                             ) 

control_threshold_scenarios_rep <- do.call('rbind', 
                                           replicate(length(variant_emergence_day) * 2, #unique control start scenarios (lag and lead) = 2 
                                                     control_threshold_scenarios, 
                                                     simplify = F
                                                     )
                                           )

# Complete set of scenarios ====
sim_table <- bind_cols(control_threshold_scenarios_rep, control_times_scenarios) %>% 
    mutate(campaign_duration = campaign_duration, 
           npi_duration = npi_duration,
           variant_emergence_day = rep(variant_emergence_day, 
                                       each = length(variant_emergence_day) * length(variant_emergence_day) * 2)
           ) %>% 
    relocate(variant_emergence_day, .before = 'phi')



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


