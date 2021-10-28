#Packages ----
library(deSolve)
library(dplyr)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')


# control parameters (turned off for uncontrolled epidemic)
no_control_parms_df <- data.frame(vax_coverage = 0, vax_rate = 0, 
                                  vax_speed = 0, vax_start = max_time, 
                                  campaign_duration = 0, npi_intensity = 0, 
                                  npi_start = max_time, npi_duration = 0
                                  )
