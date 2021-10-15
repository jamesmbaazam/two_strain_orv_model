#' two strain model A model of two co-existing strains of the same virus
#' @param t 
#' @param y 
#' @param parms 
#' @param browse 


two_strain_model <- function(t, y, parms, browse = FALSE) {
  
    if(browse) browser()
  
    with(c(as.list(y), parms), {
      
    
    # Force of infection of the wild-type strain ####
    
    FOI_Iw <- beta_w * (Iw + Imw)
    
    # Force of infection of the variant strain ####
    FOI_Im <- beta_m * (Im + Iwm)
    
    
    # NPIs control #### 
    
    npi_intensity <- ifelse(t < npi_start | t > npi_start + npi_duration | npi_intensity == 0, 0, npi_intensity)
    

    # Vaccination coverage to a rate ####
    #coverage_correction <- 0.9909 #Use this to prevent the epsilon conversion from going to infinity when coverage = 1
    
    epsilon <- ifelse(t < vax_start | t > vax_start + campaign_duration | vax_coverage == 0, 0, vax_rate)
    #vax_eligible_pop <- S + RwSm + RmSw
    # 
    # vax_rate <- vax_rate*vax_eligible_pop
    
    # epsilon <- ifelse(t >= vax_start, vax_rate, 0)
    
    
   # vax_campaign_objective <- vax_coverage*vax_eligible_pop
      
   # epsilon <- ifelse((t >= vax_start) & (vax_campaign_objective > 0), min(vax_rate, vax_campaign_objective), 0)
    
    # Model equations ####
    
    dSdt <- -(1 - npi_intensity) * FOI_Iw * S - (1 - npi_intensity) * FOI_Im * S - epsilon * S
    dIwdt <- (1 - npi_intensity) * FOI_Iw * S - gamma_w * Iw
    dImdt <- (1 - npi_intensity) * FOI_Im * S - gamma_m * Im
    dIwmdt <- (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm - gamma_m * Iwm
    dImwdt <- (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw - gamma_w * Imw
    dRwSmdt <- gamma_w * Iw - epsilon * RwSm - (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm
    dRmSwdt <- gamma_m * Im - epsilon * RmSw - (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw
    dRdt <- gamma_m * Iwm + gamma_w * Imw - epsilon * R
    dVdt <- epsilon * (S + RwSm + RmSw + R)
    dKdt <- (1 - npi_intensity) * FOI_Iw * S + (1 - npi_intensity) * FOI_Im * S + (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm + (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw #Cumulative incidence

    mod_result <- list(c(dSdt, dIwdt, dImdt, dIwmdt, dImwdt, dRwSmdt, dRmSwdt, dRdt, dVdt, dKdt))
    return(mod_result)
  })
}


#' The simulation function ====
#'
#' @param pop_inits 
#' @param dynamics_parms 
#' @param control_parms 
#' @param max_time 
#' @param dt 
#' @param events_table 
#' @param return_dynamics 
#' @param browse 

simulate_model <- function(pop_inits, dynamics_parms, 
                           control_parms, max_time, 
                           dt, events_table, 
                           return_dynamics = FALSE, browse = FALSE
                           ){
  
  if(browse)browser()

  
  # Simulation time ####
 # model_time <- 1:max_time
  
  events_table <- events_table %>% 
    mutate(time = rep(control_parms$variant_emergence_day, 
                      nrow(events_table)
                      )
           )
  
  sim_parms <- bind_cols(dynamics_parms, control_parms)
  # Model run ####
  sim_results <- as.data.frame(lsoda(pop_inits, dt, 
                                     two_strain_model, 
                                    parms = sim_parms,
                                    events = list(data = events_table)
                                    )
                              )
  
  incidence <- sim_results$Iw + sim_results$Im + sim_results$Iwm + sim_results$Imw
  peak_magnitude <- max(incidence)
  total_cases <- max(sim_results$K)
  total_vaccinated <- max(sim_results$V)
  
  
  if(return_dynamics){
    return(sim_results)
  }else{
    control_parms_df <- as_tibble(control_parms)
    
    results_df <- control_parms_df %>% mutate(total_cases = total_cases, 
                                              peak_cases = peak_magnitude, 
                                              total_vaccinated = total_vaccinated
                                              )
    return(results_df)
    }
}

# calc_campaign_duration <- function(vax_coverage, vax_rate, coverage_correction){
#   duration <- round(-log(1 - vax_coverage*coverage_correction) / vax_rate)
#   res <- data.frame(vax_rate = vax_rate, vax_coverage = vax_coverage, campaign_duration = duration)
#   return(res)
# }
