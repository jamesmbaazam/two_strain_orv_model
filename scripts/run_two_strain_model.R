
#Running the model
library(deSolve)
library(tidyverse)


#source the scripts
source('./scripts/two_strain_model.R')


### INITIALIZE PARAMETER SETTINGS

parms_no_vax_model <- c(beta_w = 1.5/7, 
                        beta_m = 2/7,
                        phi = 0,
                        gamma_w = 1/14,
                        gamma_m = 1/36,
                        epsilon = 0, 
                        sigma_w = 0,
                        sigma_m = 0,
                        variant_emergence_day = 30,
                        npi_implementation_day = 0,
                        npi_duration = 0,
                        vax_day = 0,
                        campaign_duration = 0,
                        vax_coverage = 0, 
                        coverage_correction = 0.99999
                        ) 

parms_vax_model <- c(beta_w = 1.5/7, 
                        beta_m = 2/7,
                        phi = 0,
                        gamma_w = 1/14,
                        gamma_m = 1/36,
                        epsilon = 0, 
                        sigma_w = 0,
                        sigma_m = 0,
                        variant_emergence_day = 30,
                        npi_implementation_day = 0,
                        npi_duration = 0,
                        vax_day = 35,
                        campaign_duration = 30,
                        vax_coverage = 80, 
                        coverage_correction = 0.99999
                        ) 

dt <- seq(0, 365, 1)      # set the time points for evaluation


# Initial conditions

inits <- c(S = 0.99, 
           Iw = 0.01, 
           Im = 0.01, 
           Iwm = 0,
           Imw = 0,
           RwSm = 0, 
           RmSw = 0,
           R = 0,
           V = 0
           )

### Simulation
no_vax_dynamics <- as.data.frame(lsoda(inits, dt, two_strain_model, parms = parms_no_vax_model)) 

no_vax_dynamics %>% View()

vax_dynamics <- as.data.frame(lsoda(inits, dt, two_strain_model, parms = parms_vax_model))

plot(no_vax_dynamics$time, no_vax_dynamics$S,type = 'b', col = 'green4')
lines(no_vax_dynamics$time, no_vax_dynamics$Iw,type = 'b', col = 'red')
lines(no_vax_dynamics$time, no_vax_dynamics$Im,type = 'b', col = 'purple')
lines(no_vax_dynamics$time, no_vax_dynamics$R,type = 'b', col = 'green4')
lines(no_vax_dynamics$time, no_vax_dynamics$V,type = 'b', col = 'green')
# 
# #Cumulative incidence
# no_vax_dynamics <- no_vax_dynamics %>% 
#     mutate(I = Iw + Im + Iwm + Imw)
# 
# plot(no_vax_dynamics$time, cumsum(no_vax_dynamics$I), type = 'b', col = 'tomato3')


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


