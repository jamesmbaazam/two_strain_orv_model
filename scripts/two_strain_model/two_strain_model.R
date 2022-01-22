#' Two strain model A model of two co-existing strains of the same virus
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



#' Two strain model with vaccine escape
#' @param t 
#' @param y 
#' @param parms 
#' @param browse 


ts_model_vax_escape <- function(t, y, parms, browse = FALSE) {
    
    if(browse) browser()
    
    with(c(as.list(y), parms), {
        
        
        # Force of infection of the wild-type strain ####
        
        FOI_Iw <- beta_w * (Iw + Imw + VIw)
        
        # Force of infection of the variant strain ####
        FOI_Im <- beta_m * (Im + Iwm + VIm)
        
        
        # NPIs control #### 
        
        npi_intensity <- ifelse(t < npi_start | t > npi_start + npi_duration | npi_intensity == 0, 0, npi_intensity)
        
        
        # Vaccination coverage to a rate ####
        #coverage_correction <- 0.9909 #Use this to prevent the epsilon conversion from going to infinity when coverage = 1
        
        epsilon <- ifelse(t < vax_start | t > vax_start + campaign_duration | vax_coverage == 0, 0, vax_rate)
        
        # Model equations ####
        
        dSdt <- -(1 - npi_intensity) * FOI_Iw * S - (1 - npi_intensity) * FOI_Im * S - epsilon * S
        dIwdt <- (1 - npi_intensity) * FOI_Iw * S - gamma_w * Iw
        dImdt <- (1 - npi_intensity) * FOI_Im * S - gamma_m * Im
        dIwmdt <- (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm - gamma_m * Iwm
        dImwdt <- (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw - gamma_w * Imw
        dRwSmdt <- gamma_w * Iw - epsilon * RwSm - (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm
        dRmSwdt <- gamma_m * Im - epsilon * RmSw - (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw
        dRdt <- gamma_m * Iwm + gamma_w * Imw #those who have recovered from both infections are fully protected
        dVdt <- epsilon * S - (1 - vax_efficacy_w)*(1 - npi_intensity) * FOI_Iw * V - (1 - vax_efficacy_m)*(1 - npi_intensity) * FOI_Im * V #Vaccinated individuals can still catch either infections but at a reduced rate
        dVIwdt <- (1 - vax_efficacy_w)*(1 - npi_intensity) * FOI_Iw * V - gamma_m * VIw
        dVImdt <- (1 - vax_efficacy_m)*(1 - npi_intensity) * FOI_Im * V - gamma_m * VIm
        dRwSmVdt <- epsilon * RwSm + gamma_m * VIw #those vax'd and also recovered from one infection are fully protected
        dRmSwVdt <- epsilon * RmSw + gamma_m * VIm #those vax'd and also recovered from one infection are fully protected
        
        
        #Track infection outflows to calculate cumulative incidence
        Iw_incidence <- (1 - npi_intensity) * FOI_Iw * S
        Im_incidence <- (1 - npi_intensity) * FOI_Im * S
        RwSm_incidence <- (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm
        RmSw_incidence <- (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw
        VIw_incidence <- (1 - vax_efficacy_w)*(1 - npi_intensity) * FOI_Iw * V
        VIm_incidence <- (1 - vax_efficacy_m)*(1 - npi_intensity) * FOI_Im * V
        
        
        
        #Cumulative incidence
        dKdt <- Iw_incidence + Im_incidence + RwSm_incidence + RmSw_incidence + VIw_incidence + VIm_incidence 
        
        mod_result <- list(c(dSdt, dIwdt, dImdt, dIwmdt, dImwdt, dRwSmdt, dRmSwdt, dRdt, dVdt, dVIwdt, dVImdt, dRwSmVdt, dRmSwVdt, dKdt))
        return(mod_result)
    })
}



#' Extract certain summaries from the two strain no vax escape 
#'
#' @param dynamics_df 
#'
#' @return
#' @export
#'
#' @examples
extract_model_summaries <- function(dynamics_df) {
  
  prevalence <- dynamics_df$Iw + dynamics_df$Im + dynamics_df$Iwm + dynamics_df$Imw
  
  #The final summaries
  results_df <- data.frame(variant_emergence_day = unique(dynamics_df$variant_emergence_day),
                           vax_coverage = unique(dynamics_df$vax_coverage), 
                           vax_rate = unique(dynamics_df$vax_rate), 
                           vax_speed = unique(dynamics_df$vax_speed),
                           npi_intensity = unique(dynamics_df$npi_intensity),
                           npi_duration = unique(dynamics_df$npi_duration),
                           total_cases = max(dynamics_df$K), 
                           peak_cases = max(prevalence), 
                           total_vaccinated = max(dynamics_df$V)
  )
  return(results_df)
}



#' Extract certain summaries from the two strain with vax escape
#'
#' @param dynamics_df 
#'
#' @return
#' @export
#'
#' @examples
extract_vax_escape_mod_summaries <- function(dynamics_df) {
    
    prevalence <- dynamics_df$Iw + dynamics_df$Im + dynamics_df$Iwm + dynamics_df$Imw + dynamics_df$VIw + dynamics_df$VIm
    
    # prevalence_wildtype <- dynamics_df$Iw + dynamics_df$Imw + dynamics_df$VIw
    # 
    # prevalence_variant <- dynamics_df$Im + dynamics_df$Iwm + dynamics_df$VIm
    
    
    dynamics_df$all_vaxed <- dynamics_df$V + dynamics_df$RwSmV + dynamics_df$RmSwV 
    
    #The final summaries
    results_df <- data.frame(variant_emergence_day = unique(dynamics_df$variant_emergence_day),
                             vax_coverage = unique(dynamics_df$vax_coverage), 
                             vax_rate = unique(dynamics_df$vax_rate), 
                             vax_speed = unique(dynamics_df$vax_speed),
                             npi_intensity = unique(dynamics_df$npi_intensity),
                             npi_duration = unique(dynamics_df$npi_duration),
                             total_cases = max(dynamics_df$K), 
                             peak_cases = max(prevalence), 
                             total_vaccinated = max(dynamics_df$all_vaxed)
    )
    return(results_df)
}




#' Function to run one instance of the two strain model ====
#'
#' @param pop_inits 
#' @param dynamics_parms 
#' @param control_parms 
#' @param max_time 
#' @param dt 
#' @param events_table 
#' @param return_dynamics 
#' @param browse 

simulate_raw_dynamics <- function(model_func = two_strain_model, 
                                  pop_inits, dynamics_parms, 
                                  control_parms, max_time, 
                           dt, events_table, 
                           get_summaries,
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
                                     model_func, 
                                    parms = sim_parms,
                                    events = list(data = events_table)
                                    )
                              )
  
  #Add the controls as id cols to the dynamics 
  sim_results_with_controls <- sim_results %>% 
    mutate(variant_emergence_day = control_parms$variant_emergence_day, 
           vax_coverage = control_parms$vax_cov, 
           vax_rate = control_parms$vax_rate, 
           vax_speed = control_parms$vax_speed,
           npi_intensity = control_parms$npi_intensity,
           npi_duration = control_parms$npi_duration
           )
  
  if(get_summaries){
    model_summaries <- extract_model_summaries(sim_results_with_controls)
    return(model_summaries)
  }else{
    return(sim_results_with_controls)
  }
  }


#' Function to run all the scenarios provided in a df ====
#'
#' @param sim_table #The simulation table
#'
#' @return
#' @export
#'
#' @examples
run_sim_all <- function(sim_table, model_func, get_summaries = TRUE){
    res <- sim_table %>% 
        rowwise() %>% 
        do({with(.,
                 simulate_raw_dynamics(model_func = model_func, 
                                       pop_inits = pop_inits, 
                                       dynamics_parms = dynamics_params,
                                       control_parms = .,
                                       max_time = max_time, 
                                       dt = eval_times,
                                       events_table = event_df,
                                       get_summaries = get_summaries,
                                       browse = FALSE
                 )
        )
        }) %>% 
        ungroup() %>% 
        as_tibble()
    return(res)
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
