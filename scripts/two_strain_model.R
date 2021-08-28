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
#' Event data frame for introducing mutant strain into model dynamics 
#' (check ?deSolve::event for the explanation fo the df below)
# ==============================================================================
event_df <- data.frame(var = c('S', 'Im'), 
                       time = c(200, 200), #time of mutant introduction
                       value = c(-0.01, 0.01), #index number of mutant cases
                       method = c('add', 'replace') #operation on state variables
                       )

# ==============================================================================
#' Event function for introducing mutant strain into model dynamics 
#' (check ?deSolve::event for the explanation fo the df below)
# ==============================================================================

event_function <- function(t, y, parms){
  with(c(as.list(y)), {
    S <- S - 0.01
    Im <- Im + 0.01
    
  return(list(c(S, Iw, Im, Iwm, Imw, RwSm, RmSw, R, V, K)))
    })
}



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
    
    epsilon <- ifelse(t < vax_day | t > vax_day + campaign_duration | vax_coverage == 0, 0, (-log(1 - vax_coverage*coverage_correction) / campaign_duration))
    
    
    
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

simulate_ts_model <- function(pop_inits, parms, max_time, dt){
  
  # =============================================
  # Initial conditions
  # =============================================
  inits <- pop_inits
  
  # =============================================
  # Simulation parameters
  # =============================================
  parms <- c(beta_w = 1.5/7, 
                          beta_m = 2/7,
                          phi = 0,
                          gamma_w = 1/14,
                          gamma_m = 1/36,
                          epsilon = 0, 
                          sigma_w = 0,
                          sigma_m = 0,
                          variant_emergence_day = 5,
                          npi_implementation_day = 0,
                          npi_duration = 0,
                          vax_day = 0,
                          campaign_duration = 0,
                          vax_coverage = 0, 
                          coverage_correction = 0.99999
                          )
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
  
  # =============================================
  # Final results time
  # =============================================   
  results_df <- cbind(results_df, sim_results)
  
  
  return(sim_results)
}