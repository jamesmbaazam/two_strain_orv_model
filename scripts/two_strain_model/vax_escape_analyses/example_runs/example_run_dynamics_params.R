# Source the global inputs
source('./scripts/two_strain_model/sim_config_global_params.R')

# control parameters 
control_parms_df_example_run <- data.frame(vax_coverage = 0, 
                               vax_rate = 0, 
                               vax_speed = 0, 
                               vax_start = max_time, 
                               campaign_duration = 0, 
                               npi_intensity = 0.1, 
                               npi_start = 1, 
                               npi_duration = max_time
                               )