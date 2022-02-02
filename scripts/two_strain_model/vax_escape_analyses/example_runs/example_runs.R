library(deSolve)
library(dplyr)
library(ggplot2)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model_functions.R')
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/vax_escape_analyses/example_runs/example_run_dynamics_params.R')

# Simulation table set up ====

variant_emergence_df <- data.frame(variant_emergence_day = c(1, 61, 121, 151, max_time)) #variant can emerge any day of the epidemic

#'Bind the original intervention params df to itself as many times as the number
#'of rows in the original dynamics params df. We want to bind the two dfs together to
#'represent the various scenarios in the original dynamics df
intervention_params_expanded <- control_parms_df_example_run %>%
    slice(rep(row_number(), times = nrow(dynamics_params_vax_escape)))


#'Expand the original dynamics df so that each row is as many as the number 
#'of rows in the intervention controls df. We want to bind the two resulting 
#'data frames
dynamics_params_expanded <- dynamics_params_vax_escape %>% 
    slice(rep(row_number(), each = nrow(control_parms_df_example_run)))

#' Attach the new dynamics params df to the new interventions df to form the 
#' simulation table specific for this simulation. Represents various scenarios
#' of vaccine efficacy loss against the variant
sim_df_example_run <- bind_cols(intervention_params_expanded, dynamics_params_expanded) 


sim_table_example_run <- variant_emergence_df %>% 
    mutate(sim_df_example_run)


output_example_run <- run_batch_sims_vax_escape_model(sim_table = sim_table_example_run,
                                                      get_summaries = FALSE #return the dynamics
                                                      )



#Plot

output_example_run_mod <- output_example_run %>% 
    mutate(prevalence = Iw + Imw + VIw + Im + Iwm + VIw)


example_run_plot <- ggplot(output_example_run_mod) + 
    geom_line(aes(x = time, y = prevalence, linetype=as.factor(variant_emergence_day)
                  ), size = 1
              )  
print(example_run_plot)
