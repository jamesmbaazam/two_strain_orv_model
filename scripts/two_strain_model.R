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
        
        beta_w <- R0_w*gamma_w
        beta_m <- ifelse(t < variant_emergence_time, 0, R0_m*gamma_m)
        
        
        dSdt <- - beta_w*S - beta_m*S - epsilon*S + alpha_w*R_w + alpha_m*R_m
        dIwdt <- beta_w*S - gamma_w * I_w
        dImdt <- beta_m*S - gamma_m * I_m
        dRwdt <- gamma_w * I_w - alpha_w*R_w
        dRmdt <- gamma_m * I_m - alpha_m*R_m
        dVdt <- epsilon*S
        
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
