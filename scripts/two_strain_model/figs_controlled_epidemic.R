#Packages ----
library(scales)
library(patchwork)
library(mdthemes)
library(tidyverse)
library(metR)

#helper scripts
source('./scripts/two_strain_model/sim_config_global_params.R')
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')

#plot paths for thesis and git
git_plot_path <- './figures/'
thesis_plot_path <- "C:/Users/JAMESAZAM/Dropbox/My Academic Repository/_SACEMA/Academic/_PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/"

#Load the model output
controlled_epidemic <- readRDS('./model_output/orv_npi_all_scenarios_dynamics_parallel.rds')

#Rescale the population proportions to total sizes
controlled_epidemic_rescaled <- controlled_epidemic %>%
    mutate(across(.cols = c(total_cases, peak_cases, total_vaccinated), 
                  .fns = ~ .x*target_pop #rescale the population proportions to total sizes
                  ),
           variant_emerges = ifelse(variant_emergence_day < max_time, 'yes', 'no')
           ) %>% 
    select(-c(npi_duration, total_vaccinated))





#' Towards the outbreak size isoclines

outbreak_size_isocline_df <- controlled_epidemic_rescaled %>% 
    filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time), 
           npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           total_cases <= 1000
           ) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity) %>% # add npi_intensity
    mutate(min_speed = min(vax_speed)) %>% 
    ungroup()

#Outbreak size isoclines ----
outbreak_size_isocline <- ggplot(outbreak_size_isocline_df, 
                                 aes(x = vax_coverage, 
                                     y = min_speed, 
                                     color = as.factor(variant_emergence_day)
                                     )) + 
    geom_line(size = 1, show.legend = TRUE) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.30, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(#title = paste('Cumulative cases threshold <= 1000 at various NPI intensity levels'), 
         x = 'Vaccination coverage', 
         y = 'Vaccination speed', 
         color = 'Variant emergence day'
         ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(outbreak_size_isocline)

#Save the files 
ggsave(outbreak_size_isocline,
       filename = 'outbreak_size_isocline_summary.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

ggsave(outbreak_size_isocline,
       filename = 'outbreak_size_isocline_summary.eps',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

#' Loop through the npi levels and facet by total cases target
npi_levels <- npi_intensity

plot_list <- list()

for (npi_level in 1:length(npi_levels)) {
    #subset the data
    isocline_df <- controlled_epidemic_rescaled %>% 
        filter(npi_intensity == npi_levels[npi_level], total_cases <= 1000) %>% 
        group_by(variant_emergence_day, vax_coverage) %>% # add npi_intensity
        mutate(min_speed = min(vax_speed))

    isocline_plot <- ggplot(isocline_df %>% 
                                filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time)), 
                            aes(x = vax_coverage, 
                                y = min_speed, 
                                color = as.factor(variant_emergence_day)
                                )
                            ) + 
        geom_line(size = 1, show.legend = FALSE) + 
    scale_x_continuous(labels = percent_format(), limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    labs(title = paste('Cumulative cases threshold <= 1000, NPI = ', unique(isocline_df$npi_intensity)), 
         x = 'Vaccination coverage', 
         y = 'Vaccination speed', 
         color = 'Variant emergence day'
    ) + 
        theme_bw(base_size = 14) +
        theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 
    


    plot_list[[npi_level]] <- isocline_plot

    # file_name <-  paste0('./figures/isocline_plot_npi_', npi_levels[npi_level],'.png')
    # #Save the files to a pdf
    # ggsave(isocline_plot,
    #        file = file_name,
    #        width = 23.76,
    #        height = 17.86,
    #        units = 'cm')
}


pdf("./figures/isoclines_outbreak_size.pdf")
for (npi_level in 1:length(npi_levels)) {
    print(plot_list[[npi_level]])
}
dev.off()



#Peak incidence isoclines ----

peak_incidence_isocline_df <- controlled_epidemic_rescaled %>% 
    filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time), 
           npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           peak_cases <= 300
           ) %>% 
    group_by(variant_emergence_day, vax_coverage) %>% # add npi_intensity
    mutate(min_speed = min(vax_speed))

peak_incidence_isocline <- ggplot(peak_incidence_isocline_df %>% 
                                     filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time)), 
                                 aes(x = vax_coverage, 
                                     y = min_speed, 
                                     color = as.factor(variant_emergence_day)
                                 )) + 
    geom_line(size = 1, show.legend = TRUE) + 
    scale_x_continuous(labels = percent_format()) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(#title = paste('Peak incidence <= 100 at various NPI intensity levels'), 
         x = 'Vaccination coverage', 
         y = 'Vaccination speed', 
         color = 'Variant emergence day'
    ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(peak_incidence_isocline)


#Save the files 
ggsave(peak_incidence_isocline,
       filename = 'peak_incidence_isocline_summary.png',
       path = git_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')

ggsave(peak_incidence_isocline,
       filename = 'peak_incidence_isocline_summary.eps',
       path = thesis_plot_path,
       width = 23.76,
       height = 17.86,
       units = 'cm')


#' Loop through the npi levels and facet by peak incidence target
npi_levels <- npi_intensity

plot_list <- list()

for (npi_level in 1:length(npi_levels)) {
    #subset the data
    isocline_df <- controlled_epidemic_rescaled %>% 
        filter(npi_intensity == npi_levels[npi_level], peak_cases <= 100) %>% 
        group_by(variant_emergence_day, vax_coverage) %>% # add npi_intensity
        mutate(min_speed = min(vax_speed))
    
    isocline_plot <- ggplot(isocline_df %>% 
                                filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time)), 
                            aes(x = vax_coverage, 
                                y = min_speed, 
                                color = as.factor(variant_emergence_day)
                            )
    ) + 
        geom_line(size = 1, show.legend = FALSE) + 
        scale_x_continuous(labels = percent_format(), limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1)) +
        scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
        labs(#title = paste('Peak incidence <= 100, NPI = ', unique(isocline_df$npi_intensity)), 
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day'
        ) +
        theme_bw(base_size = 14) +
        theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 
    
    
    
    plot_list[[npi_level]] <- isocline_plot
    
    # file_name <-  paste0('./figures/isocline_plot_npi_', npi_levels[npi_level],'.png')
    #Save the files to a pdf
    # ggsave(isocline_plot,
    #        file = file_name,
    #        width = 23.76,
    #        height = 17.86,
    #        units = 'cm')
}


pdf("./figures/isoclines_peak_incidence.pdf")
for (npi_level in 1:length(npi_levels)) {
    print(plot_list[[npi_level]])
}
dev.off()





#Contour plots
#' Loop through the npi levels and facet by total cases target

npi_levels <- npi_intensity

plot_list <- list()

for (npi_level in 1:length(npi_levels)) {
    plot_df <- controlled_epidemic_rescaled %>%
        filter(npi_intensity == npi_levels[npi_level])

    contour_plot <- ggplot(data = plot_df) +
        stat_contour(aes(x = vax_coverage,
                         y = vax_speed,
                         z = total_cases_log10,
                         color = variant_emerges,
                         group = variant_emergence_day
        ),
        size = 0.5
        ) +
        scale_color_viridis_d(option = "cividis") +
        scale_x_continuous(labels = percent_format()) +
        scale_y_continuous(breaks = seq(0, max(controlled_epidemic$vax_speed), 2),
                           labels = seq(0, max(controlled_epidemic$vax_speed), 2)
                           ) +
        labs(title = paste('NPI = ', npi_levels[npi_level]),
             x = 'Vaccination coverage',
             y = 'Vaccination campaign speed'
             ) +
        facet_wrap('total_cases_log10_target') +
        theme_bw(base_size = 14) +
        theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 
    


    plot_list[[npi_level]] <- contour_plot

    # file_name <-  paste0('./figures/contour_plot_npi_', npi_levels[npi_level],'.png')
    # #Save the files to a pdf
    # ggsave(contour_plot,
    #        file = file_name,
    #        width = 23.76,
    #        height = 17.86,
    #        units = 'cm')
}


pdf("./figures/summary_contour_plots.pdf")
for (npi_level in 1:length(npi_levels)) {
    print(plot_list[[npi_level]])
}
dev.off()
# 
# 
# #This aspect of the data is not plotting for some reason
# controlled_epidemic_rescaled %>% 
#     filter(vax_coverage == '0.3', 
#            vax_speed == '6', 
#            npi_intensity == 0,
#            total_cases_log10_target == '(6,7]'
#            ) %>% 
#     ggplot() +
#     geom_density_2d(aes(x = vax_coverage, 
#                      y = vax_speed,
#                      # z = total_cases_log10,
#                      color = variant_emerges,
#                      group = variant_emergence_day
#     ),
#     size = 0.5
#     ) +
#     scale_color_viridis_d(option = "cividis") +
#     scale_x_continuous(labels = percent_format()) + 
#     scale_y_continuous(breaks = seq(0, max(controlled_epidemic$vax_speed), 2), 
#                        labels = seq(0, max(controlled_epidemic$vax_speed), 2)
#     ) +
#     theme_bw(base_size = 12) +
#     labs(title = paste('NPI = ', 0),
#          x = 'Vaccination coverage',
#          y = 'Vaccination campaign speed'
#     ) 


#contour plots
# outbreak_size_by_target_contour <- ggplot(data = controlled_epidemic_rescaled)+
#     geom_contour_filled(aes(x = vax_coverage, 
#                      y = vax_speed,
#                      z = total_cases_log10,
#                      group = variant_emergence_day
#     ), binwidth = 1,
#     ) +
#     scale_x_continuous(labels = percent_format()) +
#     scale_y_continuous(breaks = seq(0, max(controlled_epidemic$vax_speed), 2), 
#                        labels = seq(0, max(controlled_epidemic$vax_speed), 2)
#     ) +
#     facet_wrap('total_cases_log10_target', ncol = 4, labeller = label_both) +
#     labs(#title = 'Impact of vaccination and NPI\'s on outbreak size', 
#          x = 'Vaccination coverage',
#          y = 'Vaccination speed'
#     ) + 
#     theme_bw(base_size = 12)
# 
# #Plot in the viewer
# print(outbreak_size_by_target_contour)


#Contour plots faceted by the interventions
# outbreak_size_by_vax_coverage_contour <- ggplot(data = controlled_epidemic_rescaled %>% 
#                                                     filter(npi_intensity %in% c(0, 0.5)))+
#     geom_contour_filled(aes(x = variant_emergence_day,
#                             y = vax_speed,
#                             z = log10(total_cases)
#     ), binwidth  = 1
#     ) +
#     scale_y_continuous(breaks = seq(0, max(controlled_epidemic$vax_speed), 2), 
#                        labels = seq(0, max(controlled_epidemic$vax_speed), 2)
#     ) +
#     facet_wrap(vax_coverage ~ npi_intensity, ncol = 4, labeller = label_both) +
#     labs(title = 'Impact of vaccination and NPI\'s on outbreak size', 
#          x = 'Variant emergence day',
#          y = 'Vaccination speed',
#          fill = 'Outbreak size (log-transformed)'
#     ) + 
#     theme_bw(base_size = 12)
# 
# #Plot in the viewer
# print(outbreak_size_by_vax_coverage_contour)
# 
# #Save the plot to file
# ggsave(plot = outbreak_size_by_vax_coverage_contour,
#        filename = './figures/outbreak_size_by_vax_coverage_contour.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
# )
# 
# 
# #Peak incidence
# peak_incidence_by_vax_coverage_contour <- ggplot(data = controlled_epidemic_rescaled %>% 
#                                                     filter(npi_intensity %in% c(0, 0.5)))+
#     geom_contour_filled(aes(x = variant_emergence_day,
#                             y = vax_speed,
#                             z = log10(peak_cases)
#     ), binwidth  = 1
#     ) +
#     scale_y_continuous(breaks = seq(0, max(controlled_epidemic$vax_speed), 2), 
#                        labels = seq(0, max(controlled_epidemic$vax_speed), 2)
#     ) +
#     facet_wrap(vax_coverage ~ npi_intensity, ncol = 4, labeller = label_both) +
#     labs(title = 'Impact of vaccination and NPI\'s on peak incidence', 
#          x = 'Variant emergence day',
#          y = 'Vaccination speed',
#          fill = 'Peak incidence (log-transformed)'
#     ) + 
#     theme_bw(base_size = 12)
# 
# #Plot in the viewer
# print(peak_incidence_by_vax_coverage_contour)
# 
# #Save the plot to file
# ggsave(plot = peak_incidence_by_vax_coverage_contour,
#        filename = './figures/peak_incidence_by_vax_coverage_contour.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
# )



#Line plot of total cases per vaccination rate
# outbreak_size_line_plot <- ggplot(data = controlled_epidemic_rescaled) + 
#     geom_line(aes(x = vax_speed, 
#                   y = total_cases,
#                   group = variant_emergence_day,
#                   color = variant_emerges,
#                   linetype = variant_emerges
#     ), size = 1
#     ) +
#     geom_text(data = controlled_epidemic_rescaled %>% 
#                   filter(variant_emergence_day %in% c(1, max_time)), 
#               aes(x = vax_speed, 
#                   y = total_cases, 
#                   group = variant_emergence_day, 
#                   label = variant_emergence_day
#               ),
#               size = 3,
#               color = 'red',
#               check_overlap = TRUE
#     ) + 
#     scale_color_viridis_d() + 
#     # scale_x_continuous(labels = scales::percent_format()) +
#     scale_y_log10(labels = comma) +
#     facet_wrap('vax_coverage') +
#    # expand_limits(x = min(vax_rate_vec)) +
#     labs(title = 'Total cases per vaccination coverage level', 
#         # subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
#          x = '**Vaccination speed**',
#          y = '**Total cases (log-transformed)**',
#          color = '**Variant emerges**',
#          linetype = '**Variant emerges**'
#     ) +
#     md_theme_minimal(base_size = 12)
# 
# 
# #Plot in the viewer
# print(outbreak_size_line_plot)
# 
# 
# #Save the plot to file
# ggsave(plot = outbreak_size_line_plot,
#        filename = './figures/outbreak_size_line_plot.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
#        )
# 
# 
# #Line plot of peak incidence per vaccination rate
# peak_incidence_line_plot <- ggplot(data = controlled_epidemic_rescaled) + 
#     geom_line(aes(x = vax_speed, 
#                   y = peak_cases,
#                   group = variant_emergence_day,
#                   color = variant_emerges,
#                   linetype = variant_emerges
#     ), size = 1
#     ) +
#     geom_text(data = controlled_epidemic_rescaled %>% 
#                   filter(variant_emergence_day %in% c(1, max_time)), 
#               aes(x = vax_speed, 
#                   y = peak_cases, 
#                   group = variant_emergence_day, 
#                   label = variant_emergence_day
#               ),
#               size = 3,
#               color = 'red',
#               check_overlap = TRUE
#     ) + 
#     scale_color_viridis_d() + 
#     # scale_x_continuous(labels = scales::percent_format()) +
#     scale_y_log10(labels = comma) +
#     facet_wrap('vax_coverage') +
#    # expand_limits(x = min(vax_rate_vec)) +
#     labs(title = 'Peak incidence per vaccination coverage level', 
#         # subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
#          x = '**Vaccination speed**',
#          y = '**Peak incidence (log-transformed)**',
#          color = '**Variant emerges**',
#          linetype = '**Variant emerges**'
#     ) +
#     md_theme_minimal(base_size = 12)
# 
# #Plot in the viewer
# print(peak_incidence_line_plot)
# 
# #Save the plot to file
# ggsave(plot = peak_incidence_line_plot,
#        filename = './figures/peak_incidence_line_plot.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
#        )


#Total vaccinations
# total_vaccinated_line_plot <- ggplot(data = controlled_epidemic_rescaled) + 
#     geom_line(aes(x = vax_speed, 
#                   y = total_vaccinated,
#                   group = variant_emergence_day,
#                   color = variant_emerges,
#                   linetype = variant_emerges
#     ), size = 1
#     ) +
#     geom_text(data = controlled_epidemic_rescaled %>% 
#                   filter(variant_emergence_day %in% c(1, max_time)), 
#               aes(x = vax_speed, 
#                   y = total_vaccinated, 
#                   group = variant_emergence_day, 
#                   label = variant_emergence_day
#               ),
#               size = 3,
#               color = 'red',
#               check_overlap = TRUE
#     ) + 
#     scale_color_viridis_d() + 
#     # scale_x_continuous(labels = scales::percent_format()) +
#     scale_y_continuous(labels = comma) +
#     facet_wrap('vax_coverage') +
#     labs(title = 'Total vaccinated individuals per vaccination coverage level', 
#          # subtitle = paste0('Campaign starts on day ', vax_start, ' with ', vax_cov*100, '% coverage objective'),
#          x = '**Vaccination speed**',
#          y = '**Total vaccinations (log-transformed)**',
#          color = '**Variant emerges**',
#          linetype = '**Variant emerges**'
#     ) +
#     md_theme_minimal(base_size = 12)
# 
# #Plot in the viewer
# print(total_vaccinated_line_plot)
# 
# #Save the plot to file
# ggsave(plot = total_vaccinated_line_plot,
#        filename = './figures/total_vaccinated_line_plot.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
#        )

#Contour plots ----

#Outbreak size
# outbreak_size_contour <- ggplot(data = controlled_epidemic_rescaled %>% 
#                                   filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
#                               ) + 
#     geom_contour_filled(aes(x = vax_speed,
#                   y = vax_coverage,
#                   z = log10(total_cases)
#                   ), binwidth  = 1
#               ) +
#     # scale_y_log10() +
#     # scale_fill_viridis_d() +
#     scale_x_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
#                        labels = seq(1, max(controlled_epidemic$vax_speed), 1)
#                        ) +
#     scale_y_continuous(labels = scales::percent_format()) +
#     facet_wrap('variant_emergence_day') +
#     expand_limits(x = c(0, 0)) +
#     labs(title = 'Outbreak size (log-transformed) per variant emergence day', 
#          x = '**Vaccination speed**',
#          y = '**Vaccination coverage**',
#          fill = '**Outbreak size**'
#          ) +
#     md_theme_bw(base_size = 12)
# 
# #Plot in the viewer
# print(outbreak_size_contour)
# 
# #Save the plot to file
# ggsave(plot = outbreak_size_contour,
#        filename = './figures/outbreak_size_contour.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
#        )
# 
# 
# # Peak incidence
# peak_incidence_contour <- ggplot(data = controlled_epidemic_rescaled %>% 
#                                   filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
#                              ) + 
#     geom_contour_filled(aes(x = vax_speed,
#                             y = vax_coverage,
#                             z = peak_cases
#     ), bins = 10
#     ) +
#     # scale_fill_viridis_d() +
#     scale_x_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
#                        labels = seq(1, max(controlled_epidemic$vax_speed), 1)
#     ) +
#     scale_y_continuous(labels = scales::percent_format(), trans = 'log') +
#     facet_wrap('variant_emergence_day') +
#     expand_limits(x = c(0, 0)) +
#     labs(title = 'Peak incidence (log-transformed) per variant emergence day', 
#          x = '**Vaccination speed**',
#          y = '**Vaccination coverage**',
#          fill = '**Peak incidence**'
#     ) +
#     md_theme_bw(base_size = 12)
# 
# #Plot in the viewer
# print(peak_incidence_contour)
# 
# #Save the plot to file
# ggsave(plot = peak_incidence_contour,
#        filename = './figures/peak_incidence_contour.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
#        )



#heatmaps
# Outbreak size
# outbreak_size_heatmap <- ggplot(data = controlled_epidemic_rescaled %>% 
#                                      filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
#                                 ) + 
#     geom_tile(aes(x = vax_speed,
#                   y = vax_coverage,
#                   fill = total_cases
#     )) +
#     # scale_fill_viridis_d() +
#     scale_x_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
#                        labels = seq(1, max(controlled_epidemic$vax_speed), 1)
#     ) +
#     scale_y_continuous(labels = scales::percent_format()) +
#     scale_fill_continuous(trans = 'log', labels = unit_format(unit = 'M', scale = 1E-6)) +
#     facet_wrap('variant_emergence_day') +
#     # expand_limits(x = c(0, 0)) +
#     labs(title = 'Outbreak size (log-transformed) per variant emergence day', 
#          x = 'Vaccination speed',
#          y = 'Vaccination coverage',
#          fill = 'Outbreak size'
#     ) + 
#     theme_bw(base_size = 12)
# 
# #Plot in the viewer
# print(outbreak_size_heatmap)
# 
# #Save the plot to file
# ggsave(plot = outbreak_size_heatmap,
#        filename = './figures/outbreak_size_heatmap.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
#        )
# 
# 
# 
# # Peak incidence
# peak_incidence_heatmap <- ggplot(data = controlled_epidemic_rescaled %>% 
#                                      filter(variant_emergence_day %in% seq(1, max_time - 1, 20))
#                                  ) + 
#     geom_tile(aes(x = vax_speed,
#                             y = vax_coverage,
#                             fill = peak_cases
#                             )) +
#     # scale_fill_viridis_d() +
#     scale_x_continuous(breaks = seq(1, max(controlled_epidemic$vax_speed), 1), 
#                        labels = seq(1, max(controlled_epidemic$vax_speed), 1)
#     ) +
#     scale_y_continuous(labels = scales::percent_format()) +
#     scale_fill_continuous(trans = 'log', labels = unit_format(unit = 'M', scale = 1E-6)) +
#     facet_wrap('variant_emergence_day') +
#     # expand_limits(x = c(0, 0)) +
#     labs(title = 'Peak incidence (log-transformed) per variant emergence day', 
#          x = 'Vaccination speed',
#          y = 'Vaccination coverage',
#          fill = 'Peak incidence'
#     ) + 
#     theme_bw(base_size = 12)
# 
# #Plot in the viewer
# print(peak_incidence_heatmap)
# 
# 
# #Save the plot to file
# ggsave(plot = peak_incidence_heatmap,
#        filename = './figures/peak_incidence_heatmap.png',
#        width = 23.76,
#        height = 17.86,
#        units = 'cm'
# )



