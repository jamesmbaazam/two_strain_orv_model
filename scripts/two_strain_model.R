#' two strain model A model of two co-existing strains of the same virus
#'
#' @param t 
#' @param y 
#' @param parms 
#'
#' @return
#' @export
#'
#' @examples
two_strain_model <- function(t, y, parms){
    
    with(c(as.list(y),parms),{
        
        #Set the beta for the variant at the time point it emerges
        beta_m <- ifelse(t < variant_emergence_day, 0, beta_m)
        
        #implement an NPI between a time period
        npi_intensity <- ifelse(t < npi_implementation_day | t > npi_implementation_day + npi_duration, 
                                0, 
                                npi_intensity
                                )
        
        epsilon <- ifelse(t < vax_day | t > vax_day + campaign_duration, 0, epsilon)
        
        dSdt <- - (1-npi_intensity)*beta_w*S*I_w - (1-npi_intensity)*beta_m*S*I_m - epsilon*S 
        dIwdt <- (1-npi_intensity)*beta_w*S*I_w - gamma_w * I_w
        dImdt <- (1-npi_intensity)*beta_m*S*I_m - gamma_m * I_m
        dRwdt <- gamma_w * I_w - epsilon*R_w
        dRmdt <- gamma_m * I_m - epsilon*R_m
        dVdt <- epsilon*S + epsilon*R_w + epsilon*R_m 
        dIcum <- (1-npi_intensity)*beta_w*S*I_w + (1-npi_intensity)*beta_m*S*I_m
        
        mod_result <- list(c(dSdt, dIwdt, dImdt, dRwdt, dRmdt, dVdt, dIcum))
        return(mod_result)
    })
}



#Running the model
library(deSolve)
library(tidyverse)


### INITIALIZE PARAMETER SETTINGS

parms_no_vax_model <- c(beta_w = 0.00002, 
           beta_m = 0.00004,
           npi_intensity = 0.00005,
           gamma_w = 0.00125,
           gamma_m = 0.00125,
           epsilon = 0, 
           variant_emergence_day = 200,
           npi_implementation_day = 90,
           npi_duration = 90,
           vax_day = 0,
           campaign_duration = 0
           ) 

parms_vax_model <- c(beta_w = 0.00002, 
                     beta_m = 0.00004,
                     npi_intensity = 0.00005,
                  gamma_w = 0.00125,
                  gamma_m = 0.00125,
                  epsilon = 0.001, 
                  variant_emergence_day = 50, 
                  npi_implementation_day = 90,
                  npi_duration = 90,
                  vax_day = 365,
                  campaign_duration = 365
                  )    

dt <- seq(0, 1825, 1)      # set the time points for evaluation


# Initial conditions

inits <- c(S = 999, 
           I_w = 1, 
           I_m = 1, 
           R_w = 0, 
           R_m = 0, 
           V = 0,
           Icum = 0
           )     


### Simulation
no_vax_dynamics <- as.data.frame(lsoda(inits, dt, two_strain_model, parms = parms_no_vax_model)) 
vax_dynamics <- no_vax_dynamics <- as.data.frame(lsoda(inits, dt, two_strain_model, parms = parms_vax_model)) 


npi_annotation_df <- data.frame(npi_start = rep(parms_vax_model['npi_implementation_day'], 1000),
                                npi_end = rep(parms_vax_model['npi_implementation_day'], 1000) + parms_vax_model['npi_duration'],
                                y_coord = 0:inits['S']
                                )

#complete dynamics
sirv_plot <- ggplot(data = no_vax_dynamics) + 
    geom_line(aes(x = time,
                  y = S,
                  color = 'S'
                  ),
              size = 1.2
              ) +
    geom_line(data = no_vax_dynamics, 
              aes(x = time, 
                  y = I_w,
                  color = 'I_w'
              ), 
              size = 1.2
    ) +
    geom_line(data = no_vax_dynamics, 
              aes(x = time, 
                  y = I_m,
                  color = 'I_m'
              ), 
              size = 1.2
    ) +
    geom_line(data = no_vax_dynamics, 
              aes(x = time, 
                  y = V,
                  color = 'V'
                  ), 
    size = 1.2
    ) +
    geom_rect(data = npi_annotation_df, 
               aes(xmin = npi_start,
                   xmax = npi_end,
                   ymin = 0,
                   ymax = y_coord
                   ),
             fill = 'turquoise2'
               ) +
    scale_color_manual(values = c('S' = 'blue', 
                                  'I_w' = 'red',
                                  'I_m' = 'yellow',
                                  'V' = 'purple'
                                  ))
    
plot(sirv_plot)



#' cumulative incidence

cum_inc_plot <- ggplot(data = no_vax_dynamics) +
    geom_line(aes(x = time,
                  y = Icum
                  ),
              color = 'red',
              linetype = 'dashed',
              size = 1.2
              ) +
    geom_line(data = vax_dynamics,
              aes(x = time,
                  y = Icum
                  )
              )

plot(cum_inc_plot)
