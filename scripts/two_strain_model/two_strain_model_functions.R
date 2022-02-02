#' Two strain model A model of two co-existing strains of the same virus
#' @param t
#' @param y
#' @param parms
#' @param browse


ts_model <- function(t, y, parms, browse = FALSE) {
  if (browse) browser()

  with(c(as.list(y), parms), {


    # Force of infection of the wild-type strain ####

    FOI_Iw <- beta_w * (Iw + Imw)

    # Force of infection of the variant strain ####
    FOI_Im <- beta_m * (Im + Iwm)


    # NPIs control ####

    npi_intensity <- ifelse(t < npi_start | t > npi_start + npi_duration | npi_intensity == 0, 0, npi_intensity)


    # Vaccination coverage to a rate ####
    # coverage_correction <- 0.9909 #Use this to prevent the epsilon conversion from going to infinity when coverage = 1

    epsilon <- ifelse(t < vax_start | t > vax_start + campaign_duration | vax_coverage == 0, 0, vax_rate)
    # vax_eligible_pop <- S + RwSm + RmSw
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
    dKdt <- (1 - npi_intensity) * FOI_Iw * S + (1 - npi_intensity) * FOI_Im * S + (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm + (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw # Cumulative incidence

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
  if (browse) browser()

  with(c(as.list(y), parms), {


    # Force of infection of the wild-type strain ####

    FOI_Iw <- beta_w * (Iw + Imw + VIw)

    # Force of infection of the variant strain ####
    FOI_Im <- beta_m * (Im + Iwm + VIm)


    # NPIs control ####

    npi_intensity <- ifelse(t < npi_start | t > npi_start + npi_duration | npi_intensity == 0, 0, npi_intensity)


    # Vaccination coverage to a rate ####
    # coverage_correction <- 0.9909 #Use this to prevent the epsilon conversion from going to infinity when coverage = 1

    epsilon <- ifelse(t < vax_start | t > vax_start + campaign_duration | vax_coverage == 0, 0, vax_rate)

    # Model equations ####

    dSdt <- -(1 - npi_intensity) * FOI_Iw * S - (1 - npi_intensity) * FOI_Im * S - epsilon * S
    dIwdt <- (1 - npi_intensity) * FOI_Iw * S - gamma_w * Iw
    dImdt <- (1 - npi_intensity) * FOI_Im * S - gamma_m * Im
    dIwmdt <- (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm - gamma_m * Iwm
    dImwdt <- (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw - gamma_w * Imw
    dRwSmdt <- gamma_w * Iw - epsilon * RwSm - (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm
    dRmSwdt <- gamma_m * Im - epsilon * RmSw - (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw
    dRdt <- gamma_m * Iwm + gamma_w * Imw # those who have recovered from both infections are fully protected
    dVdt <- epsilon * S - (1 - vax_efficacy_w) * (1 - npi_intensity) * FOI_Iw * V - (1 - vax_efficacy_m) * (1 - npi_intensity) * FOI_Im * V # Vaccinated individuals can still catch either infections but at a reduced rate
    dVIwdt <- (1 - vax_efficacy_w) * (1 - npi_intensity) * FOI_Iw * V - gamma_w * VIw
    dVImdt <- (1 - vax_efficacy_m) * (1 - npi_intensity) * FOI_Im * V - gamma_m * VIm
    dRwSmVdt <- epsilon * RwSm + gamma_w * VIw # those vax'd and also recovered from one infection are fully protected
    dRmSwVdt <- epsilon * RmSw + gamma_m * VIm # those vax'd and also recovered from one infection are fully protected


    # Track infection outflows to calculate cumulative incidence
    Iw_incidence <- (1 - npi_intensity) * FOI_Iw * S
    Im_incidence <- (1 - npi_intensity) * FOI_Im * S
    RwSm_incidence <- (1 - npi_intensity) * (1 - sigma_w) * FOI_Im * RwSm
    RmSw_incidence <- (1 - npi_intensity) * (1 - sigma_m) * FOI_Iw * RmSw
    VIw_incidence <- (1 - npi_intensity) * (1 - vax_efficacy_w) * FOI_Iw * V
    VIm_incidence <- (1 - npi_intensity) * (1 - vax_efficacy_m) * FOI_Im * V



    # Cumulative incidence
    dKdt <- Iw_incidence + Im_incidence + RwSm_incidence + RmSw_incidence + VIw_incidence + VIm_incidence

    mod_result <- list(c(dSdt, dIwdt, dImdt, dIwmdt, dImwdt, 
                         dRwSmdt, dRmSwdt, dRdt, dVdt, dVIwdt, 
                         dVImdt, dRwSmVdt, dRmSwVdt, dKdt
                         )
                       )
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
extract_summaries_ts_model <- function(dynamics_df) {
  prevalence <- dynamics_df$Iw + dynamics_df$Im + dynamics_df$Iwm + dynamics_df$Imw

  # The final summaries
  results_df <- data.frame(
    variant_emergence_day = unique(dynamics_df$variant_emergence_day),
    vax_coverage = unique(dynamics_df$vax_coverage),
    vax_rate = unique(dynamics_df$vax_rate),
    vax_speed = unique(dynamics_df$vax_speed),
    npi_intensity = unique(dynamics_df$npi_intensity),
    npi_duration = unique(dynamics_df$npi_duration),
    R0m = unique(dynamics_df$R0_m),
    vax_efficacy_w = unique(dynamics_params$vax_efficacy_w),
    vax_efficacy_m = unique(dynamics_params$vax_efficacy_m),
    cross_protection_w = unique(dynamics_params$sigma_w),
    cross_protection_m = unique(dynamics_params$sigma_m),
    total_cases = max(dynamics_df$K),
    peak_prevalence = max(prevalence),
    total_vaccinated = max(dynamics_df$V)
  )
  return(results_df)
}



#' Extract certain summaries from the full time series of dynamics from the two strain with vax escape
#'
#' @param dynamics_df
#'
#' @return
#' @export
#'
#' @examples
extract_summaries_vax_escape_ts_model <- function(dynamics_df, browse = FALSE) {
    if(browse) browser()
    
  prevalence <- dynamics_df$Iw + dynamics_df$Im + dynamics_df$Iwm + dynamics_df$Imw + dynamics_df$VIw + dynamics_df$VIm

  # prevalence_wildtype <- dynamics_df$Iw + dynamics_df$Imw + dynamics_df$VIw
  #
  # prevalence_variant <- dynamics_df$Im + dynamics_df$Iwm + dynamics_df$VIm

  #Collapse all the vaccinated classes
  dynamics_df$all_vaxed <- dynamics_df$V +  dynamics_df$VIw + dynamics_df$VIm + dynamics_df$RwSmV + dynamics_df$RmSwV

  # The final summaries
  results_df <- data.frame(
    variant_emergence_day = unique(dynamics_df$variant_emergence_day),
    vax_coverage = unique(dynamics_df$vax_coverage),
    vax_rate = unique(dynamics_df$vax_rate),
    vax_speed = unique(dynamics_df$vax_speed),
    npi_intensity = unique(dynamics_df$npi_intensity),
    npi_duration = unique(dynamics_df$npi_duration),
    R0m = unique(dynamics_df$R0m),
    vax_efficacy_w = unique(dynamics_df$vax_efficacy_w),
    vax_efficacy_m = unique(dynamics_df$vax_efficacy_m),
    cross_protection_w = unique(dynamics_df$cross_protection_w),
    cross_protection_m = unique(dynamics_df$cross_protection_m),
    total_cases = max(dynamics_df$K),
    peak_prevalence = max(prevalence),
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

simulate_dynamics_ts_model <- function(pop_inits, dynamics_parms,
                                       control_parms, max_time,
                                       dt, events_table,
                                       get_summaries,
                                       browse = FALSE) {
  if (browse) browser()


  # Simulation time ####
  # model_time <- 1:max_time

  events_table <- events_table %>%
    mutate(time = rep(
      control_parms$variant_emergence_day,
      nrow(events_table)
    ))

  sim_parms <- bind_cols(dynamics_parms, control_parms)
  # Model run ####
  sim_results <- as.data.frame(lsoda(pop_inits, dt, ts_model,
    parms = sim_parms,
    events = list(data = events_table)
  ))

  # Add the controls as id cols to the dynamics
  sim_results_with_controls <- sim_results %>%
    mutate(
      variant_emergence_day = control_parms$variant_emergence_day,
      vax_coverage = control_parms$vax_cov,
      vax_rate = control_parms$vax_rate,
      vax_speed = control_parms$vax_speed,
      npi_intensity = control_parms$npi_intensity,
      npi_duration = control_parms$npi_duration,
      R0m = dynamics_parms$R0m,
      vax_efficacy_w = dynamics_parms$vax_efficacy_w,
      vax_efficacy_m = dynamics_parms$vax_efficacy_m,
      cross_protection_w = dynamics_parms$sigma_w,
      cross_protection_m = dynamics_parms$sigma_m
    )

  if (get_summaries) {
    model_summaries <- extract_summaries_ts_model(sim_results_with_controls, dynamics_parms = dynamics_params)
    return(model_summaries)
  } else {
    return(sim_results_with_controls)
  }
}


#' Function to run one instance of the two strain model with vax escape ====
#'
#' @param pop_inits
#' @param dynamics_parms could be NULL or have elements which are binded to control_parms_df along the way
#' @param control_parms
#' @param max_time
#' @param dt
#' @param events_table #introduces index variant cases into the model dynamics at set time
#' @param get_summaries # True or False. If False, raw time series is returned
#' @param browse #If True, opens the browser for debugging

simulate_dynamics_vax_escape_ts_model <- function(pop_inits, all_parms, 
                                                  max_time, dt, 
                                                  events_table,
                                                  get_summaries, 
                                                  browse = FALSE
                                                  ) {
  if (browse) browser()


  # Simulation time ####
  # model_time <- 1:max_time

  events_table <- events_table %>%
    mutate(time = rep(all_parms$variant_emergence_day, nrow(events_table)))
  
  # Model run ####
  sim_results <- as.data.frame(lsoda(pop_inits, dt, ts_model_vax_escape, parms = all_parms, events = list(data = events_table)))

  # Add the controls as id cols to the dynamics
  sim_results_with_controls <- sim_results %>%
    mutate(
      variant_emergence_day = all_parms$variant_emergence_day,
      vax_coverage = all_parms$vax_coverage,
      vax_rate = all_parms$vax_rate,
      vax_speed = all_parms$vax_speed,
      npi_intensity = all_parms$npi_intensity,
      npi_duration = all_parms$npi_duration,
      R0m = all_parms$R0m,
      vax_efficacy_w = all_parms$vax_efficacy_w,
      vax_efficacy_m = all_parms$vax_efficacy_m,
      cross_protection_w = all_parms$sigma_w,
      cross_protection_m = all_parms$sigma_m
    )

  if (get_summaries) {
    model_summaries <- extract_summaries_vax_escape_ts_model(dynamics_df = sim_results_with_controls)
    return(model_summaries)
  } else {
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
# function to run simulations
run_batch_sims_ts_model <- function(sim_table, get_summaries) {
  res <- sim_table %>%
    rowwise() %>%
    do({
      with(
        .,
        simulate_dynamics_ts_model(
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



#' Run the model with vax escape over all rows of an NxM simulation df
#'
#' @param sim_table #A data.frame of inputs with each row representing a scenario
#' @param get_summaries # If True, returns the model summaries. If False, returns the time series per row of inputs
#'
#' @return
#' @export
#'
#' @examples
run_batch_sims_vax_escape_model <- function(sim_table, get_summaries) {
  res <- sim_table %>%
    rowwise() %>%
    do({
      with(
        .,
        simulate_dynamics_vax_escape_ts_model(
          pop_inits = pop_inits_vax_escape,
          all_parms = .,
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
