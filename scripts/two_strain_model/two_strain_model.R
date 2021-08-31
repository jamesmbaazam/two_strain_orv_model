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
#' 
#' 


# ==============================================================================
#' The two strain model
# ==============================================================================

two_strain_model <- function(t, y, parms, browse = F) {
  
    if(browse) browser()
  
    with(c(as.list(y), parms), {
      
    # ========================================  
    # Force of infection of the wild-type virus
    # ========================================
    FOI_Iw <- beta_w * (Iw + Imw)
    
    # ======================================== 
    # Introduce transmission variant transmission after t >= variant_emergence_day
    # ======================================== 
    FOI_Im <- beta_m * (Im + Iwm)
    
    
    # ======================================== 
    # Implement NPIs with intensity phi between a time period
    # ======================================== 
    
    phi <- ifelse(t < npi_implementation_day | t > npi_implementation_day + npi_duration, 0, phi)
    
    # ======================================== 
    # Convert the vaccination coverage to a rate 
    # ======================================== 
    
    epsilon <- ifelse(t < campaign_start | t > campaign_start + campaign_duration | vax_coverage == 0, 0, (-log(1 - vax_coverage*coverage_correction) / campaign_duration))
    
    
    
    # ======================================== 
    # The two strain model
    # ======================================== 
    
    dSdt <- -(1 - phi) * FOI_Iw * S - (1 - phi) * FOI_Im * S - epsilon * S
    dIwdt <- (1 - phi) * FOI_Iw * S - gamma_w * Iw
    dImdt <- (1 - phi) * FOI_Im * S - gamma_m * Im
    dIwmdt <- (1 - phi) * (1 - sigma_w) * FOI_Im * RwSm - gamma_m * Iwm
    dImwdt <- (1 - phi) * (1 - sigma_m) * FOI_Iw * RmSw - gamma_w * Imw
    dRwSmdt <- gamma_w * Iw - epsilon * RwSm - (1 - phi) * (1 - sigma_w) * FOI_Im * RwSm
    dRmSwdt <- gamma_m * Im - epsilon * RmSw - (1 - phi) * (1 - sigma_m) * FOI_Iw * RmSw
    dRdt <- gamma_m * Iwm + gamma_w * Imw
    dVdt <- epsilon * (S + RwSm + RmSw)
    dKdt <- (1 - phi) * FOI_Iw * S + (1 - phi) * FOI_Im * S + (1 - phi) * (1 - sigma_w) * FOI_Im * RwSm + (1 - phi) * (1 - sigma_m) * FOI_Iw * RmSw #Cumulative incidence

    mod_result <- list(c(dSdt, dIwdt, dImdt, dIwmdt, dImwdt, dRwSmdt, dRmSwdt, dRdt, dVdt, dKdt))
    return(mod_result)
  })
}

# ==============================================================================
#' The simulation function
# ==============================================================================

simulate_model <- function(pop_inits, parms, max_time, dt, browse = F){
  
  if(browse)browser()
  # =============================================
  # Initial conditions
  # =============================================
  inits <- pop_inits

  # =============================================
  # Simulation time
  # =============================================
  model_time <- 1:max_time
  
  
  # =============================================
  # Model run
  # =============================================  
  
  sim_results <- as.data.frame(lsoda(inits, dt, two_strain_model, 
                                    parms = parms,
                                    events = list(data = event_df)
                                    )
                              )
  
  total_cases <- max(sim_results$K)
  
  # =============================================
  # Final results time
  # =============================================   
  results_df <- parms %>% mutate(total_cases = total_cases)
  
  
  return(results_df)
}