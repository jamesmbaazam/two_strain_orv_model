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
        beta_m <- ifelse(t < variant_emergence_time, 0, beta_m)
        
        #implement an NPI between a time period
        npi_intensity <- ifelse(t < npi_implementation_time | t > npi_implementation_time + npi_duration, 
                                0, 
                                npi_intensity
                                )
        
        dSdt <- - (1-npi_intensity)*beta_w*S*I_w - (1-npi_intensity)*beta_m*S*I_w - epsilon*S 
        dIwdt <- (1-npi_intensity)*beta_w*S*I_w - gamma_w * I_w
        dImdt <- (1-npi_intensity)*beta_m*S*I_w - gamma_m * I_m
        dRwdt <- gamma_w * I_w - epsilon*R_w
        dRmdt <- gamma_m * I_m - epsilon*R_m
        dVdt <- epsilon*S + epsilon*R_w + epsilon*R_m 
        
        mod_result <- list(c(dSdt, dIwdt, dImdt, dRwdt, dRmdt, dVdt))
        return(mod_result)
    })
}



#Running the model
library(deSolve)
library(tidyverse)


### INITIALIZE PARAMETER SETTINGS

parms_no_vax_model <- c(R0_w = 1, 
           R0_m = 2,
           gamma_w = 1/8,
           gamma_m = 1/8,
           epsilon = 0, 
           alpha_w = 0,
           alpha_m = 0,
           variant_emergence_time = 50
           ) 

parms_vax_model <- c(R0_w = 1, 
                  R0_m = 2,
                  gamma_w = 1/8,
                  gamma_m = 1/8,
                  epsilon = 1/14, 
                  alpha_w = 0,
                  alpha_m = 0,
                  variant_emergence_time = 50
)    

dt <- seq(0, 365, 1)      # set the time points for evaluation


# Initial conditions

inits <- c(S = 999, 
           I_w = 1, 
           I_m = 0, 
           R_w = 0, 
           R_m = 0, 
           V = 0
           )     

# Calculation of the total number of individuals in each sub-population

N <- sum(inits)
N
### SIMULATION OF THE MODEL

## Use lsoda to solve the differential equations numerically. 

no_vax_dynamics <- as.data.frame(lsoda(inits, dt, two_strain_model, parms = parms_no_vax_model)) 
vax_dynamics <- no_vax_dynamics <- as.data.frame(lsoda(inits, dt, two_strain_model, parms = parms_vax_model)) 

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
    scale_color_manual(values = c('S' = 'blue', 
                                  'I_w' = 'red',
                                  'I_m' = 'yellow',
                                  'V' = 'purple'
                                  ))
    
print(sirv_plot)



#' cumulative incidence
no_vax_dynamics_with_cum_incidence <- no_vax_dynamics %>% mutate(I = I_w + I_m, 
                                   I_cum_no_vax = cumsum(I)
                                   ) 

vax_dynamics_with_cum_incidence <- vax_dynamics %>% mutate(I = I_w + I_m, 
                                                           I_cum_vax = cumsum(I)
                                                           ) 
    

cum_inc_plot <- ggplot(data = no_vax_dynamics_with_cum_incidence) + 
    geom_line(aes(x = time,
                  y = I_cum_no_vax
                  ),
              color = 'red', 
              size = 1.2
              ) +
    geom_line(data = vax_dynamics_with_cum_incidence, 
              aes(x = time, 
                  y = I_cum_vax
                  )
              )
    
print(cum_inc_plot)
