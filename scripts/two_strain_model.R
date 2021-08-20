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
        
        #Introduce the index variant case after t=variant_emergence_day
        beta_m <- ifelse(t < variant_emergence_day, 0, beta_m)
        
        #implement NPIs with intensity phi between a time period
        phi <- ifelse(t < npi_implementation_day | t > npi_implementation_day + npi_duration, 
                                0, 
                                phi
                                )
        
        epsilon <- ifelse(t < vax_day | t > vax_day + campaign_duration, 0, epsilon)
        
        dSdt = - (1 -  phi)*beta_w*S*Iw - (1 - phi)*beta_m*S*Im - epsilon*S  
        dIwdt = (1 - phi)*beta_w*S*Iw - gamma_w*Iw  
        dImdt = (1 -  phi)*beta_m*S*Im - gamma_m*Im  
        dIwmdt = (1 - phi)*(1 - sigma_w)*beta_m*RwSm*Im - gamma_m*Iwm  
        dImwdt = (1 - phi)*(1 - sigma_m)*beta_w*RmSw*Im - gamma_w*Imw  
        dRwSmdt = gamma_w*Iw - epsilon*RwSm - (1 - phi)*(1 - sigma_w)*beta_m*RwSm*Im 
        dRmSwdt = gamma_m*Im - epsilon*RmSw - (1 - phi)*(1 - sigma_w)*beta_w*RmSw*Iw 
        dRdt = gamma_m*Iwm + gamma_w*Imw  
        dVdt =   epsilon*(S +  RwSm + RmSw)
        
        mod_result <- list(c(dSdt, dIwdt, dImdt, dIwmdt, dImwdt, dRwSmdt, dRmSwdt, dRdt, dVdt))
        return(mod_result)
    })
}

