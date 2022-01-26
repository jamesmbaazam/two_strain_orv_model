covid_model <- function(t, y, parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    
    FOI <- gamma * R_0(t, maxR_1, f_1, phi_1) * I

    dSdt <- mu * N0 - (FOI + mu) * S + sigma_1 * R
    dEdt <- FOI * S - (nu + mu) * E
    dIdt <- nu * E - (gamma + mu) * I
    dRdt <- gamma * I - (sigma_1 + mu) * R

    list(c(dSdt, dEdt, dIdt, dRdt))
  })
}

covid2_model <- function(t, y, parms){
  # First infection is different from subsequent infections
  with(c(as.list(y),parms),{
    
    FOI <- gamma * R_0(t, maxR_1, f_1, phi_1) * (I + Ix)
    
    dSdt <- mu * N0 - (FOI + mu) * S
    dEdt <- FOI * S - (nu + mu) * E
    dIdt <- nu * E - (gamma + mu) * I
    dRdt <- gamma * I - (sigma_1 + mu) * R
    dSxdt <- sigma_1 * (R + Rx) - (FOI + mu) * Sx
    dExdt <- FOI * Sx - (nu + mu) * Ex
    dIxdt <- nu * Ex - (gamma + mu) * Ix
    dRxdt <- gamma * Ix - (sigma_1 + mu) * Rx
    
    list(c(dSdt, dEdt, dIdt, dRdt, dSxdt, dExdt, dIxdt, dRxdt))
  })
}

R_0 <- function(t, maxR = 2, f = 0.21, phi = 28){
  x <- (2*pi / 52) * (t - phi)
  maxR * ((f / 2) * cos(x) + (1 - f/2))
}

corona_model <- function(t, y, parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    
    I1 <- I1S2 + I1E2 + I1I2 + I1R2
    I2 <- S1I2 + E1I2 + I1I2 + R1I2
    FOI_1 <- gamma * R_0(t, maxR_1, f_1, phi_1) * I1
    FOI_2 <- gamma * R_0(t, maxR_2, f_2, phi_2) * I2

    dS1S2dt <- mu * N0 - (FOI_1 + FOI_2 + mu) * S1S2 + sigma_1 * R1S2 + sigma_2 * S1R2
    dE1S2dt <- FOI_1 * S1S2 + sigma_2 * E1R2 - (nu + mu) * E1S2 - (1-chi12) * FOI_2 * E1S2
    dS1E2dt <- FOI_2 * S1S2 + sigma_1 * R1E2 - (nu + mu) * S1E2 - (1-chi21) * FOI_1 * S1E2
    dI1S2dt <- nu * E1S2 + sigma_2 * I1R2 - (gamma + mu) * I1S2 - (1-chi12) * FOI_2 * I1S2
    dE1E2dt <- (1-chi12) * FOI_2 * E1S2 + (1-chi21) * FOI_1 * S1E2 - (2*nu + mu) * E1E2
    dS1I2dt <- nu * S1E2 + sigma_1 * R1I2 - (gamma + mu) * S1I2 - (1-chi21) * FOI_1 * S1I2
    dR1S2dt <- gamma * I1S2 + sigma_2 * R1R2 - (sigma_1 + mu) * R1S2 - (1-chi12) * FOI_2 * R1S2
    dI1E2dt <- (1-chi12) * FOI_2 * I1S2 + nu * E1E2 - (gamma + nu + mu) * I1E2
    dE1I2dt <- (1-chi21) * FOI_1 * S1I2 + nu * E1E2 - (gamma + nu + mu) * E1I2
    dS1R2dt <- gamma * S1I2 + sigma_1 * R1R2 - (sigma_2 + mu) * S1R2 - (1-chi21) * FOI_1 * S1R2
    dR1E2dt <- (1-chi12) * FOI_2 * R1S2 + gamma * I1E2 - (sigma_1 + nu + mu) * R1E2
    dI1I2dt <- nu * I1E2 + nu * E1I2 - (2*gamma + mu) * I1I2
    dE1R2dt <- gamma * E1I2 + (1-chi21) * FOI_1 * S1R2 - (sigma_2 + nu + mu) * E1R2
    dR1I2dt <- nu * R1E2 + gamma * I1I2 - (gamma + sigma_1 + mu) * R1I2
    dI1R2dt <- nu * E1R2 + gamma * I1I2 - (gamma + sigma_2 + mu) * I1R2
    dR1R2dt <- gamma * R1I2 + gamma * I1R2 - (sigma_1 + sigma_2 + mu) * R1R2
    
    list(c(dS1S2dt, dE1S2dt, dS1E2dt, dI1S2dt, dE1E2dt, dS1I2dt, dR1S2dt, dI1E2dt
           , dE1I2dt, dS1R2dt, dR1E2dt, dI1I2dt, dE1R2dt, dR1I2dt, dI1R2dt, dR1R2dt))
  })
}
