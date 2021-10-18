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

simulate_raw_dynamics <- function(pop_inits, dynamics_parms, 
                           control_parms, max_time, 
                           dt, events_table, 
                           browse = FALSE
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
  
  #Add the controls as id cols to the dynamics 
  sim_results_with_controls <- sim_results %>% 
    mutate(variant_emergence_day = control_parms$variant_emergence_day, 
           vax_coverage = control_parms$vax_cov, 
           vax_rate = control_parms$vax_rate, 
           vax_speed = control_parms$vax_speed
           )
  
    return(sim_results_with_controls)
}


#' Extract certain summaries from the two strain SIR model output
#'
#' @param dynamics_df 
#'
#' @return
#' @export
#'
#' @examples
extract_model_summaries <- function(dynamics_df, control_parms_df) {
  incidence <- dynamics_df$Iw + dynamics_df$Im + dynamics_df$Iwm + dynamics_df$Imw
  peak_magnitude <- max(incidence)
  total_cases <- max(dynamics_df$K)
  total_vaccinated <- max(dynamics_df$V)
  
  
#The final summaries
  results_df <- control_parms_df %>% mutate(
    total_cases = total_cases,
    peak_cases = peak_magnitude,
    total_vaccinated = total_vaccinated
  )
  return(results_df)
}


#' Calculate vaccination hazards from coverage and campaign duration 
#'
#' @param vax_coverage Vaccination coverage
#' @param coverage_correction A value for preventing the result from going to zero (see details) 
#' @param max_time #The total time of the simulation
#' @param campaign_start #The commencement time of the campaign
#'
#' @return A data.frame of the vaccination coverage, vax rate, and campaign duration trio
#' @export
#'
#' @examples calc_vax_rate(1, 0.9999, 365, 1)
calc_vax_rate <- function(vax_coverage, coverage_correction = 0.9999, max_time, campaign_start){
  
  vax_rate <- -log(1 - vax_coverage*coverage_correction) / (max_time - campaign_start)
  vax_controls_df <- data.frame(vax_rate = vax_rate, vax_coverage = vax_coverage, 
                                campaign_duration = max_time - campaign_start
                                )
  return(vax_controls_df)
}
