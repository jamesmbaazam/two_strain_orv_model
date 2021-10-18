#Packages ----
library(scales)
library(patchwork)
library(mdthemes)
library(tidyverse)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')

#Read the model output
# controlled_epidemic <- readRDS('./model_output/controlled_epidemic_dynamics.rds')
controlled_epidemic_npi0_20 <- readRDS('./model_output/orv_npi0_20_simulation_dynamics_parallel.rds')
controlled_epidemic_npi50 <- readRDS('./model_output/orv_npi50_dynamics_parallel.rds')
controlled_epidemic <- rbind(controlled_epidemic_npi0_20, controlled_epidemic_npi50)

#Rescale the population proportions to the actual sizes
# controlled_epidemic_rescaled <- controlled_epidemic %>%
#     mutate(across(.cols = S:K, .fns = ~ .x*target_pop)) #rescale the population proportions to total sizes

controlled_epidemic_rescaled <- controlled_epidemic %>%
    mutate(across(.cols = c(total_cases, peak_cases, total_vaccinated), 
                  .fns = ~ .x*target_pop
                  ), 
           variant_emerges = ifelse(variant_emergence_day < max_time, 'yes', 'no')
           ) #rescale the population proportions to total sizes


#Line plot of total cases per vaccination rate
outbreak_size_line_plot <- ggplot(data = controlled_epidemic_rescaled) + 
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
         x = '**Vaccination speed**',
         y = '**Total cases (log-transformed)**',
         color = '**Variant emerges**',
         linetype = '**Variant emerges**'
    ) +
    md_theme_minimal(base_size = 12)


#Plot in the viewer
print(outbreak_size_line_plot)


#Save the plot to file
ggsave(plot = outbreak_size_line_plot,
       filename = './figures/outbreak_size_line_plot.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Line plot of peak incidence per vaccination rate
peak_incidence_line_plot <- ggplot(data = controlled_epidemic_rescaled) + 
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
    labs(title = 'Peak incidence per vaccination coverage level', 
        # subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
         x = '**Vaccination speed**',
         y = '**Peak incidence (log-transformed)**',
         color = '**Variant emerges**',
         linetype = '**Variant emerges**'
    ) +
    md_theme_minimal(base_size = 12)

#Plot in the viewer
print(peak_incidence_line_plot)

#Save the plot to file
ggsave(plot = peak_incidence_line_plot,
       filename = './figures/peak_incidence_line_plot.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


#Total vaccinations
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
         x = '**Vaccination speed**',
         y = '**Total vaccinations (log-transformed)**',
         color = '**Variant emerges**',
         linetype = '**Variant emerges**'
    ) +
    md_theme_minimal(base_size = 12)

#Plot in the viewer
print(total_vaccinated_line_plot)

#Save the plot to file
ggsave(plot = total_vaccinated_line_plot,
       filename = './figures/total_vaccinated_line_plot.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )

#Contour plots ----

#Outbreak size
outbreak_size_contour <- ggplot(data = controlled_epidemic_rescaled %>% 
                                  filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
                              ) + 
    geom_contour_filled(aes(x = vax_speed,
                  y = vax_coverage,
                  z = log10(total_cases)
                  ), binwidth  = 1
              ) +
    # scale_y_log10() +
    # scale_fill_viridis_d() +
    scale_x_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
                       labels = seq(1, max(controlled_epidemic$vax_speed), 1)
                       ) +
    scale_y_continuous(labels = scales::percent_format()) +
    facet_wrap('variant_emergence_day') +
    expand_limits(x = c(0, 0)) +
    labs(title = 'Outbreak size (log-transformed) per variant emergence day', 
         x = '**Vaccination speed**',
         y = '**Vaccination coverage**',
         fill = '**Outbreak size**'
         ) +
    md_theme_bw(base_size = 12)

#Plot in the viewer
print(outbreak_size_contour)

#Save the plot to file
ggsave(plot = outbreak_size_contour,
       filename = './figures/outbreak_size_contour.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )


# Peak incidence
peak_incidence_contour <- ggplot(data = controlled_epidemic_rescaled %>% 
                                  filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
                             ) + 
    geom_contour_filled(aes(x = vax_speed,
                            y = vax_coverage,
                            z = peak_cases
    ), bins = 10
    ) +
    # scale_fill_viridis_d() +
    scale_x_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
                       labels = seq(1, max(controlled_epidemic$vax_speed), 1)
    ) +
    scale_y_continuous(labels = scales::percent_format(), trans = 'log') +
    facet_wrap('variant_emergence_day') +
    expand_limits(x = c(0, 0)) +
    labs(title = 'Peak incidence (log-transformed) per variant emergence day', 
         x = '**Vaccination speed**',
         y = '**Vaccination coverage**',
         fill = '**Peak incidence**'
    ) +
    md_theme_bw(base_size = 12)

#Plot in the viewer
print(peak_incidence_contour)

#Save the plot to file
ggsave(plot = peak_incidence_contour,
       filename = './figures/peak_incidence_contour.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )



#heatmaps
# Outbreak size
outbreak_size_heatmap <- ggplot(data = controlled_epidemic_rescaled %>% 
                                     filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
                                ) + 
    geom_tile(aes(x = vax_speed,
                  y = vax_coverage,
                  fill = total_cases
    )) +
    # scale_fill_viridis_d() +
    scale_x_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
                       labels = seq(1, max(controlled_epidemic$vax_speed), 1)
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_continuous(trans = 'log', labels = unit_format(unit = 'M', scale = 1E-6)) +
    facet_wrap('variant_emergence_day') +
    # expand_limits(x = c(0, 0)) +
    labs(title = 'Outbreak size (log-transformed) per variant emergence day', 
         x = 'Vaccination speed',
         y = 'Vaccination coverage',
         fill = 'Outbreak size'
    ) + 
    theme_bw(base_size = 12)

#Plot in the viewer
print(outbreak_size_heatmap)

#Save the plot to file
ggsave(plot = outbreak_size_heatmap,
       filename = './figures/outbreak_size_heatmap.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
       )



# Peak incidence
peak_incidence_heatmap <- ggplot(data = controlled_epidemic_rescaled %>% 
                                     filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
                                 ) + 
    geom_tile(aes(x = vax_speed,
                            y = vax_coverage,
                            fill = peak_cases
                            )) +
    # scale_fill_viridis_d() +
    scale_x_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
                       labels = seq(1, max(controlled_epidemic$vax_speed), 1)
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_continuous(trans = 'log', labels = unit_format(unit = 'M', scale = 1E-6)) +
    facet_wrap('variant_emergence_day') +
    # expand_limits(x = c(0, 0)) +
    labs(title = 'Peak incidence (log-transformed) per variant emergence day', 
         x = 'Vaccination speed',
         y = 'Vaccination coverage',
         fill = 'Peak incidence'
    ) + 
    theme_bw(base_size = 12)

#Plot in the viewer
print(peak_incidence_heatmap)


#Save the plot to file
ggsave(plot = peak_incidence_heatmap,
       filename = './figures/peak_incidence_heatmap.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
)


#
outbreak_size_by_vax_coverage_contour <- ggplot(data = controlled_epidemic_rescaled %>% 
           filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
       ) + 
    geom_contour_filled(aes(x = variant_emergence_day,
                            y = vax_speed,
                            z = log10(total_cases)
    ), binwidth  = 1
    ) +
    scale_y_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
                       labels = seq(1, max(controlled_epidemic$vax_speed), 1)
    ) +
    facet_wrap('vax_coverage') +
    labs(title = 'Outbreak size (log-transformed) per vaccination coverage', 
         x = 'Variant emergence day',
         y = 'Vaccination coverage',
         fill = 'Outbreak size'
    ) + 
    theme_bw(base_size = 12)

#Plot in the viewer
print(outbreak_size_by_vax_coverage_contour)

#Save the plot to file
ggsave(plot = outbreak_size_by_vax_coverage_contour,
       filename = './figures/outbreak_size_by_vax_coverage_contour.png',
       width = 23.76,
       height = 17.86,
       units = 'cm'
)
