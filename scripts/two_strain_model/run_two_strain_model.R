
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


# Run the model
# no_vax_dynamics <- as.data.frame(lsoda(inits, eval_times, two_strain_model, parms = parms_no_vax_model,
#                                        events = list(data = event_df)
#                                        )
#                                  ) 

# no_vax_dynamics_with_event_func <- as.data.frame(lsoda(inits, eval_times, two_strain_model, parms = parms_no_vax_model,
#                                        events = list(func = event_function, time = 200)
#                                        )
#                                        ) 



# ========================
# Controlled epidemic (Vaccination and NPIs)
# ========================
# parms_vax_model <- c(beta_w = 1.5/7, 
#                      beta_m = 2/7,
#                      phi = 0.6,
#                      gamma_w = 1/14,
#                      gamma_m = 1/36,
#                      sigma_w = 0,
#                      sigma_m = 0,
#                      #variant_emergence_day = 150,
#                      npi_implementation_day = 21,
#                      npi_duration = 180,
#                      vax_day = 180,
#                      campaign_duration = 180,
#                      vax_coverage = 0.5, #Specify as a number between 0 and 1
#                      coverage_correction = 0.99999
#                      ) 

# Run model
# vax_dynamics <- as.data.frame(lsoda(inits, eval_times, two_strain_model, parms = parms_vax_model,
#                                     events = list(data = event_df)))
# 
# # vax_dynamics %>% View()
# 
# #Plot the uncontrolled and controlled epidemic on the same frames
# epidemic_complete <- left_join(no_vax_dynamics %>% select(time, K),
#                                vax_dynamics %>% select(time, K), 
#                                by = 'time', 
#                                suffix = c('_uncontrolled','_controlled'))
# 
# 
# epidemic_complete_long <- epidemic_complete %>% 
#     pivot_longer(cols = starts_with('K'), 
#                  names_to = 'epi_type', 
#                  values_to = 'Icum'
#     )
# 
# controlled_vrs_uncontrolled_epi_plot <- ggplot(data = epidemic_complete_long, 
#                                                aes(x = time, y = Icum)
# ) + geom_line(aes(color = epi_type), 
#               size = 1.2) +
#     labs(title = paste('NPIs implemented on day', 
#                        as.numeric(parms_vax_model["npi_implementation_day"]), 
#                        'for', as.numeric(parms_vax_model['npi_duration']), 
#                        'days AND\n vaccination implemented on day', 
#                        as.numeric(parms_vax_model['vax_day']),
#                        'for', 
#                        as.numeric(parms_vax_model['campaign_duration']),
#                        'days'
#     ))
# 
# print(controlled_vrs_uncontrolled_epi_plot)


# # ===========================================================
# # Reshape model results
# # ===========================================================
# 
# vax_dynamics_long <- vax_dynamics %>% 
#     pivot_longer(cols = S:V, 
#                  names_to = 'state', 
#                  values_to = 'number'
#     ) 

# ===========================================================
# Plot model dynamics
# ===========================================================

# Full dynamics
# two_strain_dynamics_controlled_plot <- ggplot(vax_dynamics_long, 
#                                    aes(x = time, y = number)
# ) + 
#     geom_line(aes(color = state), 
#               size = 1.2
#     ) +
#     theme_minimal()
# 
# print(two_strain_dynamics_controlled_plot)
# 
# # Cumulative incidence (controlled epidemic)
# 
# Icum_controlled_epidemic_plot <- ggplot(data = vax_dynamics, 
#        aes(x = time, 
#            y = K
#        )
#        ) + 
#     geom_line(color = 'red', size = 1.2) +
#     scale_y_continuous(breaks = seq(0, ceiling(max(vax_dynamics$K)), 0.2),
#                        labels = seq(0, ceiling(max(vax_dynamics$K)), 0.2)
#     ) +
#     labs(y = 'Cumulative incidence (with intervention)') +
#     theme_minimal()
# 
# print(Icum_controlled_epidemic_plot)


# no_vax_dynamics %>% View()

# ===========================================================
# Reshape model results
# ===========================================================

# no_vax_dynamics_long <- no_vax_dynamics %>% 
#     pivot_longer(cols = S:V, 
#                  names_to = 'state', 
#                  values_to = 'number'
#                  ) 

# ===========================================================
# Plot model dynamics
# ===========================================================

# Full dynamics
# two_strain_dynamics_uncontrolled_plot <- ggplot(no_vax_dynamics_long, 
#        aes(x = time, y = number)
#        ) + 
#     geom_line(aes(color = state), 
#               size = 1.2
#               ) +
#     theme_minimal()
# 
# print(two_strain_dynamics_uncontrolled_plot)
# 
# 
# 
# # Cumulative incidence (no vaccination, no NPI)
# Icum_uncontrolled_epidemic_plot <- ggplot(data = no_vax_dynamics, 
#        aes(x = time, 
#            y = K
#        )
# ) + 
#     geom_line(color = 'red', size = 1.2) +
#     scale_y_continuous(breaks = seq(0, ceiling(max(no_vax_dynamics$K)), 0.2),
#                        labels = seq(0, ceiling(max(no_vax_dynamics$K)), 0.2)
#     ) +
#     labs(y = 'Cumulative incidence') +
#     theme_minimal()
# 
# print(Icum_uncontrolled_epidemic_plot)


#function to calculate vaccinate rates (epsilon) from vaccination coverage and campaign duration
# epsilon <- function(vax_coverage,
#                     campaign_duration,
#                     coverage_correction = 0.99999
#                     ){
#     sim_epsilon <- -log(1 - vax_coverage*coverage_correction) / campaign_duration
# 
#     results <- data.frame(vax_cov = vax_coverage,
#                       campaign_duration = campaign_duration,
#                       vax_rate = sim_epsilon
#                       )
#     return(results)
#     }
# 
# epsilon(vax_coverage = 0.5, campaign_duration = 180)

# campaign_objectives_df <- tibble(campaign_duration = (170 + seq(10, 180, 30)), 
#                                      vax_cov = rep(0.75, times = length(campaign_duration))
#                                      )
# 
# 
# vax_rate_df <- campaign_objectives_df %>% 
#     rowwise() %>% 
#     with(., {
#     epsilon(vax_coverage = vax_cov, campaign_duration = campaign_duration)
# })

#NPI stringency
# phi <- seq(0, 0.75, length.out = 4)


# intervention_params_df <- expand.grid(phi, vax_rate_df)



#Check the infection dynamics
# infection_dynamics <- no_vax_dynamics %>% filter(time < 101)
# plot(infection_dynamics$time, infection_dynamics$Iw,type = 'b', col = 'red')
# lines(infection_dynamics$time, infection_dynamics$Im,type = 'b', col = 'purple')




#' npi_annotation_df <- data.frame(npi_start = rep(parms_vax_model['npi_implementation_day'], 1000),
#'                                 npi_end = rep(parms_vax_model['npi_implementation_day'], 1000) + parms_vax_model['npi_duration'],
#'                                 y_coord = 0:inits['S']
#' )
#' 
#' #complete dynamics
#' sirv_plot <- ggplot(data = no_vax_dynamics) + 
#'     geom_line(aes(x = time,
#'                   y = S,
#'                   color = 'S'
#'     ),
#'     size = 1.2
#'     ) +
#'     geom_line(data = no_vax_dynamics, 
#'               aes(x = time, 
#'                   y = I_w,
#'                   color = 'I_w'
#'               ), 
#'               size = 1.2
#'     ) +
#'     geom_line(data = no_vax_dynamics, 
#'               aes(x = time, 
#'                   y = I_m,
#'                   color = 'I_m'
#'               ), 
#'               size = 1.2
#'     ) +
#'     geom_line(data = no_vax_dynamics, 
#'               aes(x = time, 
#'                   y = V,
#'                   color = 'V'
#'               ), 
#'               size = 1.2
#'     ) +
#'     geom_rect(data = npi_annotation_df, 
#'               aes(xmin = npi_start,
#'                   xmax = npi_end,
#'                   ymin = 0,
#'                   ymax = y_coord
#'               ),
#'               fill = 'turquoise2'
#'     ) +
#'     scale_color_manual(values = c('S' = 'blue', 
#'                                   'I_w' = 'red',
#'                                   'I_m' = 'yellow',
#'                                   'V' = 'purple'
#'     ))
#' 
#' plot(sirv_plot)
#' 
#' 
#' 
#' #' cumulative incidence
#' 
#' cum_inc_plot <- ggplot(data = no_vax_dynamics) +
#'     geom_line(aes(x = time,
#'                   y = Icum
#'     ),
#'     color = 'red',
#'     linetype = 'dashed',
#'     size = 1.2
#'     ) +
#'     geom_line(data = vax_dynamics,
#'               aes(x = time,
#'                   y = Icum
#'               )
#'     )
#' 
#' plot(cum_inc_plot)


