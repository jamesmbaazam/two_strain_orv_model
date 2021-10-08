#Packages ----
library(deSolve)
library(scales)
library(patchwork)
library(tidyverse)

# Helper scripts ----
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')






# Controlled epidemic

#Baseline; no variant emerges and only wild type prevails
baseline_no_variant <- baseline_params %>% 
    rowwise() %>% 
    do({with(., 
             simulate_model(pop_inits = pop_inits, 
                            dynamics_parms = dynamics_params,
                            control_parms = .,
                            max_time = max_time, 
                            dt = eval_times,
                            events_table = event_df,
                            browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()  

#Scale the proportions to totals
baseline_no_variant_rescaled <- baseline_no_variant %>% 
    mutate(total_cases = total_cases*target_pop, 
           peak_cases = peak_cases*target_pop
           ) %>% 
    select(-c(vax_start, vax_coverage, starts_with('npi_')))


#Simulations with variant emergence
orv_full_simulation <- simulation_table %>% 
    rowwise() %>% 
    do({with(., 
             simulate_model(pop_inits = pop_inits, 
                            dynamics_parms = dynamics_params,
                            control_parms = .,
                            max_time = max_time, 
                            dt = eval_times,
                            events_table = event_df,
                            browse = F
             )
    )
    }) %>% 
    ungroup() %>% 
    as_tibble()

#Rescale the population proportions to the actual sizes
controlled_epidemic_rescaled <- orv_full_simulation %>%
    mutate(total_cases = total_cases*target_pop, 
           peak_cases = peak_cases*target_pop,
           emergence_emerges = ifelse(variant_emergence_day < max_time, 
                             'yes',
                             'no'
                             )
           ) %>% 
    select(-c(vax_start, vax_coverage, starts_with('npi_')))


#Line plot of total cases per vaccination rate
total_cases_line_plot <- ggplot(data = controlled_epidemic_rescaled) + 
    geom_line(aes(x = vax_rate, 
                  y = total_cases,
                  group = variant_emergence_day,
                  color = emergence_emerges,
                  linetype = emergence_emerges
                  ), size = 1
              ) +
    geom_text(data = controlled_epidemic_rescaled %>% 
                  filter(variant_emergence_day %in% c(1, max_time)), 
              aes(x = vax_rate, 
                  y = total_cases, 
                  group = variant_emergence_day, 
                  label = variant_emergence_day
                  ),
              size = 3,
              color = 'red',
              check_overlap = TRUE
              ) + 
    scale_color_viridis_d() + 
    scale_x_continuous(labels = scales::percent_format()) +
    scale_y_continuous(labels = comma) +
    expand_limits(x = min(vax_rate_vec)) +
    labs(title = 'Total cases per variant emergence day', 
         subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
         x = 'Vaccination rate',
         y = 'Total cases',
         color = 'Variant emerges',
         linetype = 'Variant emerges'
         ) +
    theme_minimal(base_size = 12)

print(total_cases_line_plot)



#Line plot of peak daily cases per vaccination rate
peak_daily_cases_line_plot <- ggplot(data = controlled_epidemic_rescaled) + 
    geom_line(aes(x = vax_rate, 
                  y = peak_cases,
                  group = variant_emergence_day,
                  color = emergence_emerges,
                  linetype = emergence_emerges
                  ), size = 1
    ) +
    geom_text(data = controlled_epidemic_rescaled %>% 
                  filter(variant_emergence_day %in% c(1, max_time)), 
              aes(x = vax_rate, 
                  y = peak_cases, 
                  group = variant_emergence_day, 
                  label = variant_emergence_day
              ),
              size = 3,
              color = 'red',
              check_overlap = TRUE
    ) + 
    scale_color_viridis_d() + 
    scale_x_continuous(labels = scales::percent_format()) +
    scale_y_continuous(labels = comma) +
    expand_limits(x = min(vax_rate_vec)) +
    labs(title = 'Peak daily cases per variant emergence day', 
         subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
         x = 'Vaccination rate',
         y = 'Peak daily cases',
         color = 'Variant emerges',
         linetype = 'Variant emerges'
    ) +
#    coord_flip() +
    theme_minimal(base_size = 12)

print(peak_daily_cases_line_plot)

#Contour plot
# ggplot(data = controlled_epidemic_rescaled) + 
#     geom_contour(aes(x = vax_rate, 
#                   y = variant_emergence_day,
#                   z = total_cases
#                   ),
#                  size = 1
#               ) +
#     scale_x_continuous(labels = scales::percent_format(),
#                        expand = c(0,0)
#     ) +
#     theme_minimal(base_size = 12)


# Select the columns for the plot
# orv_full_simulation_plot <- ggplot(data = controlled_epidemic_rescaled) + 
#     geom_tile(aes(x = vax_rate, 
#                   y = campaign_duration, 
#                   fill = total_cases
#     ), 
#     stat = 'identity'
#     ) +
#     scale_x_continuous(labels = scales::percent_format(), 
#                        expand = c(0,0)
#     ) +
#     scale_y_continuous(labels = scales::percent_format(), 
#                        expand = c(0,0)
#     ) +
#     scale_fill_viridis_b(direction = -1) +
#     labs(x = 'Vaccination coverage', 
#          y = 'NPI intensity', 
#          fill = 'Cumulative incidence'
#     ) +
#     #facet_wrap(variant_emergence_day) + 
#     theme_minimal()
# 
# plot(orv_full_simulation_plot)