# setup ########

library(deSolve) # ODE numerical integration

# Default model (fixed NPIs, +60% transmissibility/40% immune escape, leaky immunity) ########

  ## model function ########

leaky_model <- function(pop_size, vax_daily, R0, NPI_levels, threshold, lag, NPI_off, recov_rate, transm_adv, cross_protect, vax_dist_1, vax_dist_2, max_effic, infected_elig, recovered_elig, max_coverage, tau_1, tau_2, tau_vax, duration, delta){ # NPI_levels = efficacy for three states: level 1, level 2, OFF
  
  #-----------------------------------
  
  model<- function(t ,start, pars){
    
    with(as.list( c(start, pars)),{
      
      E <- S + kI*(I1S + I2S) + kR*(R1 + R2 + R12) + kI*kR*(I1R + I2R)
      
      # E = eligible for vaccination
      
      p <- if( (E > (1 - coverage)*N) & (t < (tstart + coverage/v)) ){v*N/E} else{0}
      
      # p = per capita rate of vaccination among eligible compartments
      
      I1 <- I1S + I1S_V + I1R + I1R_V + I1V_2 + I1V_0
      
      I2 <- I2S + I2S_V + I2R + I2R_V + I2V_1 + I2V_0
      
      dSdt <- -p*S - (1-r)*beta*S*I1 - (1-r)*(1 + lambda)*beta*S*I2
      
      dI1Sdt <- (1-r)*beta*S*I1 - kI*p*I1S - gamma*I1S
      
      dI1S_Vdt <- kI*p*I1S - gamma*I1S_V
      
      dI2Sdt <- (1-r)*(1 + lambda)*beta*S*I2 - kI*p*I2S - gamma*I2S
      
      dI2S_Vdt <- kI*p*I2S - gamma*I2S_V
      
      dR1dt <- gamma*I1S - kR*p*R1 - (1-r)*(1 - phi)*(1 + lambda)*beta*R1*I2
      
      dR2dt <- gamma*I2S - kR*p*R2 - (1-r)*(1 - phi)*beta*R2*I1
      
      dR12dt <- gamma*(I1R + I2R) - kR*p*R12
      
      dV1dt <- kR*p*R1 + gamma*I1S_V + gamma*I1V_0 - (1-r)*(1 - phi)*(1 - omega*(1 - alpha2))*(1 + lambda)*beta*V1*I2
      
      dV2dt <- kR*p*R2 + gamma*I2S_V + gamma*I2V_0 - (1-r)*(1 - phi)*(1 - omega*(1 - alpha1))*beta*V2*I1
      
      dV12dt <- kR*p*R12 + gamma*(I1R_V + I2R_V + I1V_2 + I2V_1)
      
      dV0dt <- p*S - (1-r)*(1 - omega*(1 - alpha1))*beta*V0*I1 - (1-r)*(1 - omega*(1 - alpha2))*(1 + lambda)*beta*V0*I2
      
      dI1Rdt <- (1-r)*(1 - phi)*beta*R2*I1 - kI*kR*p*I1R - gamma*I1R
      
      dI1R_Vdt <- kI*kR*p*I1R - gamma*I1R_V
      
      dI2Rdt <- (1-r)*(1 - phi)*(1 + lambda)*beta*R1*I2 - kI*kR*p*I2R - gamma*I2R
      
      dI2R_Vdt <- kI*kR*p*I2R - gamma*I2R_V
      
      dI1V_2dt <- (1-r)*(1 - phi)*(1 - omega*(1 - alpha1))*beta*V2*I1 - gamma*I1V_2
      
      dI2V_1dt <- (1-r)*(1 - phi)*(1 - omega*(1 - alpha2))*(1 + lambda)*beta*V1*I2 - gamma*I2V_1
      
      dI1V_0dt <- (1-r)*(1 - omega*(1 - alpha1))*beta*V0*I1 - gamma*I1V_0
      
      dI2V_0dt <- (1-r)*(1 - omega*(1 - alpha2))*(1 + lambda)*beta*V0*I2 - gamma*I2V_0
      
      dY1dt <- (1-r)*beta*I1*(S + (1 - phi)*R2 + (1 - omega*(1 - alpha1))*(V0 + (1 - phi)*V2))
      
      # Y1 = total strain 1 infections (active + recovered)
      
      dY2dt <- (1-r)*(1 + lambda)*beta*I2*(S + (1 - phi)*R1 + (1 - omega*(1 - alpha2))*(V0 + (1 - phi)*V1))
      
      # Y2 = total strain 2 infections (active + recovered)
      
      return(list(c(dSdt, dI1Sdt, dI1S_Vdt, dI2Sdt, dI2S_Vdt, dR1dt, dR2dt, dR12dt, dV1dt, dV2dt, dV12dt, dV0dt, dI1Rdt, dI1R_Vdt, dI2Rdt, dI2R_Vdt, dI1V_2dt, dI2V_1dt, dI1V_0dt, dI2V_0dt, dY1dt, dY2dt)))
      
    })
  }
  
  #-----------------------------------
  
  transm_rate <- R0*recov_rate/pop_size
  
  start <- c(S = pop_size, I1S = 0, I1S_V = 0, I2S = 0, I2S_V = 0, R1 = 0, R2 = 0, R12 = 0, V1 = 0, V2 = 0, V12 = 0, V0 = 0, I1R = 0, I1R_V = 0, I2R = 0, I2R_V = 0, I1V_2 = 0, I2V_1 = 0, I1V_0 = 0, I2V_0 = 0, Y1 = 0, Y2 = 0)
  
  pars <- c(N = pop_size, v = 0, beta = transm_rate, gamma = recov_rate, lambda = transm_adv, phi = cross_protect, alpha1 = vax_dist_1, alpha2 = vax_dist_2, omega = max_effic, kI = infected_elig, kR = recovered_elig, coverage = max_coverage, tstart = 0, r = NPI_levels[1]) # r = reduction in transmission by control measures
  
  NPI <- rep(1, duration + 1) # record of NPI "levels" over time
  
  switch <- 0 # most recent NPI on/off switch 
  
  t_off <- if(is.na(tau_vax)|is.na(NPI_off)){duration}else{min(duration, tau_vax + ceiling(NPI_off/vax_daily) + 1)} # the +1 here is not for indexing purposes - it is so that t_off = the first day with NPIs turned off, as opposed to the last day with NPIs on
  
  out <- matrix(NA, nrow = duration + 1, ncol = length(start) + 1)
  
  colnames(out) <- c("time", "S", "I1S", "I1S_V", "I2S", "I2S_V", "R1", "R2", "R12", "V1", "V2", "V12", "V0", "I1R", "I1R_V", "I2R", "I2R_V", "I1V_2", "I2V_1", "I1V_0", "I2V_0", "Y1", "Y2")
  
  #----------------------------------- 
  
  update_conditions <- function(t){
    
    if(!is.na(tau_1) & tau_1 == t){
      
      start[c("S","I1S","Y1")] <<- start[c("S","I1S","Y1")] + c(-1, 1, 1)
      
    }
    
    if(!is.na(tau_2) & tau_2 == t){
      
      start[c("S","I2S","Y2")] <<- start[c("S","I2S","Y2")] + c(-1, 1, 1)
      
    } 
    
    if(!is.na(tau_vax) & tau_vax == t){
      
      pars["v"] <<- vax_daily
      
      pars["tstart"] <<- t
      
    }
    
    lag_state <- if(t < lag){NULL}else{
      
      if(lag == 0){start}else{out[t - lag + 1,]}
      
    }
    
    I_sum <- sum(lag_state[c("I1S", "I1S_V", "I2S", "I2S_V", "I1R", "I1R_V", "I2R", "I2R_V", "I1V_2", "I2V_1", "I1V_0", "I2V_0")])
    
    if(t >= t_off){ NPI[t + 1] <<- 3}else{
      
      NPI_last <- if(t == 0){ 1 }else{ NPI[t] }
      
      NPI[t+1] <<- NPI_last
      
      if((t >= (switch + lag)) & (NPI_last == 1) & (I_sum > threshold*pars["N"])){
        
        NPI[t+1] <<- 2; switch <<- t
        
      }
      
      if((t >= (switch + lag)) & (NPI_last == 2) & (I_sum < threshold*pars["N"])){
        
        NPI[t+1] <<- 1; switch <<- t
        
      }
      
    }
    
    pars["r"] <<- NPI_levels[NPI[t + 1]]
    
  }
  
  #---------------------------------
  
  update_conditions(t = 0); out[1, 2:ncol(out)] <- start
  
  for(i in 1:duration){
    
    t <- seq(i - 1, i, by = delta)
    
    update_conditions(t = i - 1)
    
    output <- ode(start, t, model, pars, method = "adams"); rm(t)
    
    out[i+1,] <- output[nrow(output),]
    
    start <- output[nrow(output), 2:ncol(output)]; rm(output)
    
  }
  
  
  #-----------------------------------
  
  return(list(out, NPI))
  
}

  ## iterative model function ########

iterate_leaky_model <- function(pop_size = 1e8, vax_daily_vals, R0 = 2.5, NPI_levels, threshold = 0.01, lag = 14, NPI_off = NA, recov_rate = 1/14, transm_adv, cross_protect, vax_dist_1 = 0, vax_dist_2, max_effic = 0.95, infected_elig = 1, recovered_elig = 1, max_coverage = 1, tau_1 = 0, tau_2_vals, tau_vax_vals, duration, delta = 0.1){
  
  out_array <- array(NA, dim = c(length(vax_daily_vals), length(tau_2_vals), length(tau_vax_vals), duration + 1, 23), dimnames = list(NULL,NULL,NULL,NULL,c("time","S", "I1S", "I1S_V", "I2S", "I2S_V", "R1", "R2", "R12", "V1", "V2", "V12", "V0", "I1R", "I1R_V", "I2R", "I2R_V", "I1V_2", "I2V_1", "I1V_0", "I2V_0", "Y1", "Y2")))
  
  NPI_array <- array(NA, dim = c(length(vax_daily_vals), length(tau_2_vals), length(tau_vax_vals), duration + 1))
  
  for(i in 1:length(vax_daily_vals)){
    
    for(j in 1:length(tau_2_vals)){
      
      for(k in 1:length(tau_vax_vals)){
        
        results <- leaky_model(pop_size = pop_size, vax_daily = vax_daily_vals[i], R0 = R0, NPI_levels = NPI_levels, threshold = threshold, lag = lag, NPI_off = NPI_off, recov_rate = recov_rate, transm_adv = transm_adv, cross_protect = cross_protect, vax_dist_1 = vax_dist_1, vax_dist_2 = vax_dist_2, max_effic = max_effic, infected_elig = infected_elig, recovered_elig = recovered_elig, max_coverage = max_coverage, tau_1 = tau_1, tau_2 = tau_2_vals[j], tau_vax = tau_vax_vals[k], duration = duration, delta = delta)
        
        out_array[i,j,k,,] <- results[[1]]
        
        NPI_array[i,j,k,] <- results[[2]]
        
      }
      
    }
    
  }
  
  outputs <- list(out_array, NPI_array)
  
  return(outputs)
  
}

  ## simulation parameters ########

vax_daily_vals <- 12/(365*(3:12))

tau_2_vals <- round(9*365/12, digits = 0)

tau_vax_vals <- round((9:18)*365/12, digits = 0)

  ## Simulations - full control measures ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Default model/Full control measures/")

    ### variant strains alone ########

W <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_2_vals = NA, tau_vax_vals = NA, duration = 1095) # W and V0 are equivalent

#save(W, file = "W"); rm(W)

V1 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V1, file = "V1"); rm(V1)

V2 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V2, file = "V2"); rm(V2)

V3 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V3, file = "V3"); rm(V3)

    ### WT + variant strains ########

W_V0 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V0, file = "W_V0"); rm(W_V0)

W_V1 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V1, file = "W_V1"); rm(W_V1)

W_V2 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V2, file = "W_V2"); rm(W_V2)

W_V3 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V3, file = "W_V3"); rm(W_V3)

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - NPIs lifted at 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Default model/NPIs stopped early (50% vax coverage)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Default model/Incomplete vax coverage (50%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 70% vax efficacy ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Default model/Lower vaccine efficacy (70%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)
  ## Simulations - triple threat ########
  
  setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Default model/Triple threat/")
  
    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

# Decreased immune escape ########

  ## Simulations - full control measures ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Decreased escape/Full control measures/")

    ### variant strains alone ########

W <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(W, file = "W"); rm(W)

V1 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V1, file = "V1"); rm(V1)

V2 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.8, vax_dist_2 = 0.2, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V2, file = "V2"); rm(V2)

V3 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.8, vax_dist_2 = 0.2, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V3, file = "V3"); rm(V3)

    ### WT + variant strains ########

W_V0 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V0, file = "W_V0"); rm(W_V0)

W_V1 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V1, file = "W_V1"); rm(W_V1)

W_V2 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.8, vax_dist_2 = 0.2, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V2, file = "W_V2"); rm(W_V2)

W_V3 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.8, vax_dist_2 = 0.2, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V3, file = "W_V3"); rm(W_V3)

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.8, vax_dist_2 = 0.2, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.8, vax_dist_2 = 0.2, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - NPIs lifted at 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Decreased escape/NPIs stopped early (50% vax coverage)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.8, vax_dist_2 = 0.2, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.8, vax_dist_2 = 0.2, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Decreased escape/Incomplete vax coverage (50%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.8, vax_dist_2 = 0.2, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.8, vax_dist_2 = 0.2, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 70% vax efficacy ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Decreased escape/Lower vaccine efficacy (70%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.8, max_effic = 0.7, vax_dist_2 = 0.2, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.8, max_effic = 0.7, vax_dist_2 = 0.2, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - triple threat ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Decreased escape/Triple threat/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.8, max_effic = 0.7, vax_dist_2 = 0.2, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.8, max_effic = 0.7, vax_dist_2 = 0.2, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

# Increased immune escape ########

  ## Simulations - full control measures ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Increased escape/Full control measures/")

    ### variant strains alone ########

W <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 2190)

#save(W, file = "W"); rm(W)

V1 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 2190)

#save(V1, file = "V1"); rm(V1)

V2 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.2, vax_dist_2 = 0.8, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 2190)

#save(V2, file = "V2"); rm(V2)

V3 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.2, vax_dist_2 = 0.8, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 2190)

#save(V3, file = "V3"); rm(V3)

    ### WT + variant strains ########

W_V0 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 2190)

#save(W_V0, file = "W_V0"); rm(W_V0)

W_V1 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 2190)

#save(W_V1, file = "W_V1"); rm(W_V1)

W_V2 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.2, vax_dist_2 = 0.8, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 2190)

#save(W_V2, file = "W_V2"); rm(W_V2)

W_V3 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.2, vax_dist_2 = 0.8, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 2190)

#save(W_V3, file = "W_V3"); rm(W_V3)

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.2, vax_dist_2 = 0.8, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.2, vax_dist_2 = 0.8, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - NPIs lifted at 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Increased escape/NPIs stopped early (50% vax coverage)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.2, vax_dist_2 = 0.8, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.2, vax_dist_2 = 0.8, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Increased escape/Incomplete vax coverage (50%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.2, vax_dist_2 = 0.8, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.2, vax_dist_2 = 0.8, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 70% vax efficacy ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Increased escape/Lower vaccine efficacy (70%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.2, max_effic = 0.7, vax_dist_2 = 0.8, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.2, max_effic = 0.7, vax_dist_2 = 0.8, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - triple threat ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Increased escape/Triple threat/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.2, max_effic = 0.7, vax_dist_2 = 0.8, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.2, max_effic = 0.7, vax_dist_2 = 0.8, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 2190)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

# Rolling lockdowns ########

  ## Simulations - full control measures ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Rolling lockdowns/Full control measures/")

    ### variant strains alone ########

W <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_2_vals = NA, tau_vax_vals = NA, duration = 1095) # W and V0 are equivalent

#save(W, file = "W"); rm(W)

V1 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V1, file = "V1"); rm(V1)

V2 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V2, file = "V2"); rm(V2)

V3 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V3, file = "V3"); rm(V3)

    ### WT + variant strains ########

W_V0 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V0, file = "W_V0"); rm(W_V0)

W_V1 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V1, file = "W_V1"); rm(W_V1)

W_V2 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V2, file = "W_V2"); rm(W_V2)

W_V3 <- iterate_leaky_model(vax_daily_vals = NA, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V3, file = "W_V3"); rm(W_V3)

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - NPIs lifted at 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Rolling lockdowns/NPIs stopped early (50% vax coverage)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Rolling lockdowns/Incomplete vax coverage (50%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 70% vax efficacy ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Rolling lockdowns/Lower vaccine efficacy (70%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), transm_adv = 0.6, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - triple threat ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/Rolling lockdowns/Triple threat/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.3, 0.7, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

# All-or-nothing immunity ########

  ## model function ########

alt_model <- function(pop_size, vax_daily, R0, NPI_levels, threshold, lag, NPI_off, recov_rate, transm_adv, cross_protect, vax_dist_1, vax_dist_2, max_effic, infected_elig, recovered_elig, max_coverage, tau_1, tau_2, tau_vax, duration, delta){ # NPI_levels = efficacy for three states: level 1, level 2, OFF
  
  #-----------------------------------
  
  model<- function(t ,start, pars){
    
    with(as.list( c(start, pars)),{
      
      E <- S + kI*(I1S + I2S) + kR*(R1 + R2 + R12) + kI*kR*(I1R + I2R)
      
      # E = eligible for vaccination
      
      p <- if( (E > (1 - coverage)*N) & (t < (tstart + coverage/v)) ){v*N/E} else{0}
      
      # p = per capita rate of vaccination among eligible compartments
      
      I1 <- I1S + I1S_V + I1R + I1R_V + I1V_2 + I1V_0
      
      I2 <- I2S + I2S_V + I2R + I2R_V + I2V_1 + I2V_0
      
      dSdt <- -p*S - (1-r)*beta*S*I1 - (1-r)*(1 + lambda)*beta*S*I2
      
      dI1Sdt <- (1-r)*beta*S*I1 - kI*p*I1S - gamma*I1S
      
      dI1S_Vdt <- kI*p*I1S - gamma*I1S_V
      
      dI2Sdt <- (1-r)*(1 + lambda)*beta*S*I2 - kI*p*I2S - gamma*I2S
      
      dI2S_Vdt <- kI*p*I2S - gamma*I2S_V
      
      dR1dt <- (1 - phi)*gamma*I1S - kR*p*R1 - (1-r)*(1 + lambda)*beta*R1*I2
      
      dR2dt <- (1 - phi)*gamma*I2S - kR*p*R2 - (1-r)*beta*R2*I1
      
      dR12dt <- gamma*(I1R + I2R + phi*I1S + phi*I2S) - kR*p*R12
      
      dV1dt <- omega*(1 - alpha1)*alpha2*p*S + (1 - omega*(1 - alpha2))*kR*p*R1 + (1 - phi)*(1 - omega*(1 - alpha2))*gamma*I1S_V + (1 - phi)*gamma*I1V_0 - (1-r)*(1 + lambda)*beta*V1*I2
      
      dV2dt <- omega*alpha1*(1 - alpha2)*p*S + (1 - omega*(1 - alpha1))*kR*p*R2 + (1 - phi)*(1 - omega*(1 - alpha1))*gamma*I2S_V + (1 - phi)*gamma*I2V_0 - (1-r)*beta*V2*I1
      
      dV12dt <- omega*(1 - alpha1)*(1 - alpha2)*p*S + omega*(1 - alpha2)*kR*p*R1 + omega*(1 - alpha1)*kR*p*R2 + kR*p*R12 + (phi + (1 - phi)*omega*(1 - alpha2))*gamma*I1S_V + (phi + (1 - phi)*omega*(1 - alpha1))*gamma*I2S_V + gamma*(I1R_V + I2R_V + I1V_2 + I2V_1 + phi*I1V_0 + phi*I2V_0)
      
      dV0dt <- (1 - omega*(1 - alpha1*alpha2))*p*S - (1-r)*beta*V0*I1 - (1-r)*(1 + lambda)*beta*V0*I2
      
      dI1Rdt <- (1-r)*beta*R2*I1 - kI*kR*p*I1R - gamma*I1R
      
      dI1R_Vdt <- kI*kR*p*I1R - gamma*I1R_V
      
      dI2Rdt <- (1-r)*(1 + lambda)*beta*R1*I2 - kI*kR*p*I2R - gamma*I2R
      
      dI2R_Vdt <- kI*kR*p*I2R - gamma*I2R_V
      
      dI1V_2dt <- (1-r)*beta*V2*I1 - gamma*I1V_2
      
      dI2V_1dt <- (1-r)*(1 + lambda)*beta*V1*I2 - gamma*I2V_1
      
      dI1V_0dt <- (1-r)*beta*V0*I1 - gamma*I1V_0
      
      dI2V_0dt <- (1-r)*(1 + lambda)*beta*V0*I2 - gamma*I2V_0
      
      dY1dt <- (1-r)*beta*(S + R2 + V2 + V0)*I1 
      
      # Y1 = total strain 1 infections (active + recovered)
      
      dY2dt <- (1-r)*(1 + lambda)*beta*(S + R1 + V1 + V0)*I2
      
      # Y2 = total strain 2 infections (active + recovered)
      
      return(list(c(dSdt, dI1Sdt, dI1S_Vdt, dI2Sdt, dI2S_Vdt, dR1dt, dR2dt, dR12dt, dV1dt, dV2dt, dV12dt, dV0dt, dI1Rdt, dI1R_Vdt, dI2Rdt, dI2R_Vdt, dI1V_2dt, dI2V_1dt, dI1V_0dt, dI2V_0dt, dY1dt, dY2dt)))
      
    })
  }
  
  #-----------------------------------
  
  transm_rate <- R0*recov_rate/pop_size
  
  start <- c(S = pop_size, I1S = 0, I1S_V = 0, I2S = 0, I2S_V = 0, R1 = 0, R2 = 0, R12 = 0, V1 = 0, V2 = 0, V12 = 0, V0 = 0, I1R = 0, I1R_V = 0, I2R = 0, I2R_V = 0, I1V_2 = 0, I2V_1 = 0, I1V_0 = 0, I2V_0 = 0, Y1 = 0, Y2 = 0)
  
  pars <- c(N = pop_size, v = 0, beta = transm_rate, gamma = recov_rate, lambda = transm_adv, phi = cross_protect, alpha1 = vax_dist_1, alpha2 = vax_dist_2, omega = max_effic, kI = infected_elig, kR = recovered_elig, coverage = max_coverage, tstart = 0, r = NPI_levels[1]) # r = reduction in transmission by control measures
  
  NPI <- rep(1, duration + 1) # record of NPI "levels" over time
  
  switch <- 0 # most recent NPI on/off switch 
  
  t_off <- if(is.na(tau_vax)|is.na(NPI_off)){duration}else{min(duration, tau_vax + ceiling(NPI_off/vax_daily) + 1)} # the +1 here is not for indexing purposes - it is so that t_off = the first day with NPIs turned off, as opposed to the last day with NPIs on
  
  out <- matrix(NA, nrow = duration + 1, ncol = length(start) + 1)
  
  colnames(out) <- c("time", "S", "I1S", "I1S_V", "I2S", "I2S_V", "R1", "R2", "R12", "V1", "V2", "V12", "V0", "I1R", "I1R_V", "I2R", "I2R_V", "I1V_2", "I2V_1", "I1V_0", "I2V_0", "Y1", "Y2")
  
  #----------------------------------- 
  
  update_conditions <- function(t){
    
    if(!is.na(tau_1) & tau_1 == t){
      
      start[c("S","I1S","Y1")] <<- start[c("S","I1S","Y1")] + c(-1, 1, 1)
      
    }
    
    if(!is.na(tau_2) & tau_2 == t){
      
      start[c("S","I2S","Y2")] <<- start[c("S","I2S","Y2")] + c(-1, 1, 1)
      
    } 
    
    if(!is.na(tau_vax) & tau_vax == t){
      
      pars["v"] <<- vax_daily
      
      pars["tstart"] <<- t
      
    }
    
    lag_state <- if(t < lag){NULL}else{
      
      if(lag == 0){start}else{out[t - lag + 1,]}
      
    }
    
    I_sum <- sum(lag_state[c("I1S", "I1S_V", "I2S", "I2S_V", "I1R", "I1R_V", "I2R", "I2R_V", "I1V_2", "I2V_1", "I1V_0", "I2V_0")])
    
    if(t >= t_off){ NPI[t + 1] <<- 3}else{
      
      NPI_last <- if(t == 0){ 1 }else{ NPI[t] }
      
      NPI[t+1] <<- NPI_last
      
      if((t >= (switch + lag)) & (NPI_last == 1) & (I_sum > threshold*pars["N"])){
        
        NPI[t+1] <<- 2; switch <<- t
        
      }
      
      if((t >= (switch + lag)) & (NPI_last == 2) & (I_sum < threshold*pars["N"])){
        
        NPI[t+1] <<- 1; switch <<- t
        
      }
      
    }
    
    pars["r"] <<- NPI_levels[NPI[t + 1]]
    
  }
  
  #---------------------------------
  
  update_conditions(t = 0); out[1, 2:ncol(out)] <- start
  
  for(i in 1:duration){
    
    t <- seq(i - 1, i, by = delta)
    
    update_conditions(t = i - 1)
    
    output <- ode(start, t, model, pars, method = "adams"); rm(t)
    
    out[i+1,] <- output[nrow(output),]
    
    start <- output[nrow(output), 2:ncol(output)]; rm(output)
    
  }
  
  
  #-----------------------------------
  
  return(list(out, NPI))
  
}

  ## iterative model function ########

iterate_altmodel <- function(pop_size = 1e8, vax_daily_vals, R0 = 2.5, NPI_levels, threshold = 0.01, lag = 14, NPI_off = NA, recov_rate = 1/14, transm_adv, cross_protect, vax_dist_1 = 0, vax_dist_2, max_effic = 0.95, infected_elig = 1, recovered_elig = 1, max_coverage = 1, tau_1 = 0, tau_2_vals, tau_vax_vals, duration, delta = 0.1){
  
  out_array <- array(NA, dim = c(length(vax_daily_vals), length(tau_2_vals), length(tau_vax_vals), duration + 1, 23), dimnames = list(NULL,NULL,NULL,NULL,c("time","S", "I1S", "I1S_V", "I2S", "I2S_V", "R1", "R2", "R12", "V1", "V2", "V12", "V0", "I1R", "I1R_V", "I2R", "I2R_V", "I1V_2", "I2V_1", "I1V_0", "I2V_0", "Y1", "Y2")))
  
  NPI_array <- array(NA, dim = c(length(vax_daily_vals), length(tau_2_vals), length(tau_vax_vals), duration + 1))
  
  for(i in 1:length(vax_daily_vals)){
    
    for(j in 1:length(tau_2_vals)){
      
      for(k in 1:length(tau_vax_vals)){
        
        results <- alt_model(pop_size = pop_size, vax_daily = vax_daily_vals[i], R0 = R0, NPI_levels = NPI_levels, threshold = threshold, lag = lag, NPI_off = NPI_off, recov_rate = recov_rate, transm_adv = transm_adv, cross_protect = cross_protect, vax_dist_1 = vax_dist_1, vax_dist_2 = vax_dist_2, max_effic = max_effic, infected_elig = infected_elig, recovered_elig = recovered_elig, max_coverage = max_coverage, tau_1 = tau_1, tau_2 = tau_2_vals[j], tau_vax = tau_vax_vals[k], duration = duration, delta = delta)
        
        out_array[i,j,k,,] <- results[[1]]
        
        NPI_array[i,j,k,] <- results[[2]]
        
      }
      
    }
    
  }
  
  outputs <- list(out_array, NPI_array)
  
  return(outputs)
  
}

  ## Simulations - full control measures ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/All-or-nothing immunity/Full control measures/")

    ### variant strains alone ########

W <- iterate_altmodel(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_2_vals = NA, tau_vax_vals = NA, duration = 1095)

#save(W, file = "W"); rm(W)

V1 <- iterate_altmodel(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V1, file = "V1"); rm(V1)

V2 <- iterate_altmodel(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V2, file = "V2"); rm(V2)

V3 <- iterate_altmodel(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = NA, tau_2_vals = 0, tau_vax_vals = NA, duration = 1095)

#save(V3, file = "V3"); rm(V3)

    ### WT + variant strains ########

W_V0 <- iterate_altmodel(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V0, file = "W_V0"); rm(W_V0)

W_V1 <- iterate_altmodel(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V1, file = "W_V1"); rm(W_V1)

W_V2 <- iterate_altmodel(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V2, file = "W_V2"); rm(W_V2)

W_V3 <- iterate_altmodel(vax_daily_vals = NA, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = NA, duration = 1095)

#save(W_V3, file = "W_V3"); rm(W_V3)

    ### WT + variants + vax ########

W_V0_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - NPIs lifted at 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/All-or-nothing immunity/NPIs stopped early (50% vax coverage)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 50% vax coverage ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/All-or-nothing immunity/Incomplete vax coverage (50%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - 70% vax efficacy ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/All-or-nothing immunity/Lower vaccine efficacy (70%)/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_altmodel(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), transm_adv = 0.6, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

  ## Simulations - triple threat ########

setwd("C:/Users/mab4629/OneDrive - Harvard University/COVID/Two-strain vaccine model/Results/All-or-nothing immunity/Triple threat/")

    ### WT + variants + vax ########

W_V0_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V0_vax, file = "W_V0_vax"); rm(W_V0_vax)

W_V1_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 1, max_effic = 0.7, vax_dist_2 = 0, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V1_vax, file = "W_V1_vax"); rm(W_V1_vax)

W_V2_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V2_vax, file = "W_V2_vax"); rm(W_V2_vax)

W_V3_vax <- iterate_leaky_model(vax_daily_vals = vax_daily_vals, NPI_levels = c(0.4, 0.4, 0), NPI_off = 0.5, transm_adv = 0.6, cross_protect = 0.6, max_effic = 0.7, vax_dist_2 = 0.4, max_coverage = 0.5, tau_1 = 0, tau_2_vals = tau_2_vals, tau_vax_vals = tau_vax_vals, duration = 1095)

#save(W_V3_vax, file = "W_V3_vax"); rm(W_V3_vax)

