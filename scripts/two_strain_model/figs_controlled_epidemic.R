#Packages ----
library(scales)
library(patchwork)
library(tidyverse)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')

#Read the model output
controlled_epidemic <- readRDS('./model_output/simulation_controlled_epidemic_parallel_run.rds')


#Rescale the population proportions to the actual sizes
controlled_epidemic_rescaled <- controlled_epidemic %>%
    mutate(total_cases = total_cases*target_pop, 
           peak_cases = peak_cases*target_pop,
           total_vaccinated = total_vaccinated*target_pop,
           variant_emerges = ifelse(variant_emergence_day < max_time, 
                                      'yes',
                                      'no'
           )
    ) %>% 
    select(-c(vax_start, starts_with('npi_')))
    # select(-c(vax_start, vax_coverage, starts_with('npi_')))



#Line plot of total cases per vaccination rate
total_cases_line_plot <- ggplot(data = controlled_epidemic_rescaled %>% 
                                    filter(variant_emergence_day %in% seq(1, max_time, 40)
                                           )
                                ) + 
    geom_line(aes(x = vax_speed, 
                  y = total_cases,
                  group = variant_emergence_day,
                  color = variant_emerges,
                  linetype = variant_emerges
    ), size = 1
    ) +
    geom_text(data = controlled_epidemic_rescaled %>% 
                  filter(variant_emergence_day %in% c(1, max_time)), 
              aes(x = vax_speed, 
                  y = total_cases, 
                  group = variant_emergence_day, 
                  label = variant_emergence_day
              ),
              size = 3,
              color = 'red',
              check_overlap = TRUE
    ) + 
    scale_color_viridis_d() + 
    # scale_x_continuous(labels = scales::percent_format()) +
    scale_y_log10(labels = comma) +
    facet_wrap('vax_coverage') +
   # expand_limits(x = min(vax_rate_vec)) +
    labs(title = 'Total cases per vaccination coverage level', 
        # subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
         x = 'Vaccination speed',
         y = 'Total cases',
         color = 'Variant emerges',
         linetype = 'Variant emerges'
    ) +
    theme_minimal(base_size = 12)

print(total_cases_line_plot)

ggsave(plot = total_cases_line_plot,
       filename = './figures/total_cases_line_plot.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
)


#Line plot of peak daily cases per vaccination rate
peak_daily_cases_line_plot <- ggplot(data = controlled_epidemic_rescaled %>% 
                                         filter(variant_emergence_day %in% seq(1, max_time, 40))
                                     ) + 
    geom_line(aes(x = vax_speed, 
                  y = peak_cases,
                  group = variant_emergence_day,
                  color = variant_emerges,
                  linetype = variant_emerges
    ), size = 1
    ) +
    geom_text(data = controlled_epidemic_rescaled %>% 
                  filter(variant_emergence_day %in% c(1, max_time)), 
              aes(x = vax_speed, 
                  y = peak_cases, 
                  group = variant_emergence_day, 
                  label = variant_emergence_day
              ),
              size = 3,
              color = 'red',
              check_overlap = TRUE
    ) + 
    scale_color_viridis_d() + 
    # scale_x_continuous(labels = scales::percent_format()) +
    scale_y_log10(labels = comma) +
    facet_wrap('vax_coverage') +
   # expand_limits(x = min(vax_rate_vec)) +
    labs(title = 'Peak daily cases per vaccination coverage level', 
        # subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
         x = 'Vaccination speed',
         y = 'Peak daily cases',
         color = 'Variant emerges',
         linetype = 'Variant emerges'
    ) +
    #    coord_flip() +
    theme_minimal(base_size = 12)

print(peak_daily_cases_line_plot)

ggsave(plot = peak_daily_cases_line_plot,
       filename = './figures/peak_daily_cases_line_plot.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
)



total_vaccinated_line_plot <- ggplot(data = controlled_epidemic_rescaled) + 
    geom_line(aes(x = vax_speed, 
                  y = total_vaccinated,
                  group = variant_emergence_day,
                  color = variant_emerges,
                  linetype = variant_emerges
    ), size = 1
    ) +
    geom_text(data = controlled_epidemic_rescaled %>% 
                  filter(variant_emergence_day %in% c(1, max_time)), 
              aes(x = vax_speed, 
                  y = total_vaccinated, 
                  group = variant_emergence_day, 
                  label = variant_emergence_day
              ),
              size = 3,
              color = 'red',
              check_overlap = TRUE
    ) + 
    scale_color_viridis_d() + 
    # scale_x_continuous(labels = scales::percent_format()) +
    scale_y_continuous(labels = comma) +
    facet_wrap('vax_coverage') +
    labs(title = 'Total vaccinated individuals per vaccination coverage level', 
         # subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
         x = 'Vaccination speed',
         y = 'Total vaccinations',
         color = 'Variant emerges',
         linetype = 'Variant emerges'
    ) +
    #    coord_flip() +
    theme_minimal(base_size = 12)

print(total_vaccinated_line_plot)



ggplot(data = controlled_epidemic_rescaled) + 
    geom_contour(aes(x = vax_rate, y = variant_emergence_day, z = total_cases))
