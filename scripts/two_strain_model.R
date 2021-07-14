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
        beta_m <- ifelse(t < 25, beta_w, R0_m*gamma_m)
        
        
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
library(ggplot2)
library(tibble)

### INITIALIZE PARAMETER SETTINGS

parms <- c(R0_w = 1, 
           R0_m = 2,
           gamma_w = 1/8,
           gamma_m = 1/8,
           epsilon = 1/45, 
           alpha_w = 1/360,
           alpha_m = 1/360
           )    

dt    <- seq(0, 100, 1)      # set the time points for evaluation


# Initial conditions

inits <- c(S = 1000, 
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

out <- as_tibble(lsoda(inits, dt, two_strain_model, parms = parms)) 


two_strain_model_plot <- ggplot(data = out) + 
    geom_line(aes(x = time,
                  y = S,
                  color = 'S'
                  ),
              size = 1.2
              ) +
    geom_line(data = out, 
              aes(x = time, 
                  y = I_w,
                  color = 'I_w'
              ), 
              size = 1.2
    ) +
    geom_line(data = out, 
              aes(x = time, 
                  y = I_m,
                  color = 'I_m'
              ), 
              size = 1.2
    ) +
    geom_line(data = out, 
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
    
print(two_strain_model_plot)
