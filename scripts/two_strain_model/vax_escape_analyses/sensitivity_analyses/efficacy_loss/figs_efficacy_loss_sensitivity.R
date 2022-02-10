#Packages ----
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(mdthemes)

# helper scripts
source("./scripts/two_strain_model/two_strain_model_functions.R")
source("./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/global_params.R")
source("./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/intervention_params_table_setup_global.R")

# plot paths for thesis and git
git_plot_path <- "./figures/"
thesis_plot_path <- "C:/Users/JAMESAZAM/Dropbox/My Academic Repository/_SACEMA/Academic/_PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/"


# Change how large numbers are printed
options(scipen = 10000)

# Load the model output

# Original analysis
#' The baseline analysis assumes perfect vaccine efficacy against the wild-type
#' and variant and perfect cross protection between the strains. The variant is
#' assumed to be 30% more transmissible than the wild-type (R0w = 1.3R0m)
# baseline_analysis_results <- readRDS('./model_output/sensitivity_analyses/baseline_analysis_vR0_30_percent/vR0_30_percent_baseline_analysis_summaries.rds') %>%
#     rename(peak_prevalence = peak_cases)

# Cross protection sensitivity analysis results
efficacy_loss_sensitivity_analysis_summaries <- readRDS('./model_output/sensitivity_analyses/efficacy_loss/efficacy_loss_sensitivity_analysis_summaries.rds')


# Remove redundant columns
efficacy_loss_sensitivity_analysis_df <- efficacy_loss_sensitivity_analysis_summaries %>%
  select(-c(
    vax_rate,
    npi_duration,
    starts_with("cross_protection"),
    R0m,
    total_vaccinated
  ))
#' Determine the subset of scenarios that meet the threshold outbreak size
#' and the minimum speed required

outbreak_size_thresholds <- c(1E3, 5E5, 5E6, 10E6)

#' Determine the subset of scenarios that meet the threshold outbreak size
#' and the minimum speed required
#'

for (i in seq_along(outbreak_size_thresholds)) {
  efficacy_loss_cases_threshold <- outbreak_size_thresholds[i] / target_pop

  outbreak_size_efficacy_loss_isocline_df <- efficacy_loss_sensitivity_analysis_df %>%
    filter(
      variant_emergence_day %in% c(1, 61, 121, 151, max_time),
      npi_intensity %in% c(0.0, 0.1, 0.2, 0.3)
    ) %>%
    filter(total_cases <= efficacy_loss_cases_threshold) %>%
    mutate(
      variant_emergence_day = as_factor(variant_emergence_day),
      vax_efficacy_m = as_factor(vax_efficacy_m)
    ) %>%
    group_by(variant_emergence_day, vax_coverage, npi_intensity, vax_efficacy_m) %>%
    mutate(min_speed = min(vax_speed)) %>%
    ungroup() 
      # mutate(min_speed = min(vax_speed),  
      #        max_speed = 10
      #        ) %>% 
      # ungroup() %>%
    # group_by(variant_emergence_day, vax_efficacy_m) %>% 
    # mutate(min_coverage = min(vax_coverage), 
    #        max_min_speed = max(min_speed)
    #        ) %>% 
    # ungroup() %>% 
    # relocate(c(min_speed, max_min_speed, max_speed), .after = vax_speed)

  # Isoclines ----
  # Outbreak size ====
  outbreak_size_efficacy_loss_isocline_sensitivity <- ggplot(
    outbreak_size_efficacy_loss_isocline_df,
    aes(
      x = vax_coverage,
      y = min_speed,
      color = variant_emergence_day
    )
  ) +
    geom_line(aes(linetype = vax_efficacy_m),
      size = 1.3,
      show.legend = TRUE
    ) +
      # geom_segment(data = outbreak_size_efficacy_loss_isocline_df, 
      #              aes(x = min_coverage, 
      #                  y = max_min_speed,
      #                  xend = min_coverage,
      #                  yend = max_speed,
      #                  linetype = vax_efficacy_m
      #              ),
      #              size = 1.3,
      #              show.legend = TRUE
      # ) +
    scale_x_continuous(labels = percent_format(
      accuracy = 1,
      suffix = ""
    ), breaks = seq(0.10, 1, 0.1)) +
    scale_y_continuous(
      breaks = seq(1, 10, 1),
      labels = seq(1, 10, 1),
      limits = c(1, 10)
    ) +
    scale_color_viridis_d(option = "viridis") +
    labs(
      title = "Sensitivity to vaccine efficacy against variant assumptions",
      subtitle = paste0(
        "Strategies with cumulative cases up to ",
        efficacy_loss_cases_threshold * 100,
        "% of total population"
      ),
      x = "Vaccination coverage (%)",
      y = "Vaccination speed",
      color = "Variant emergence day",
      linetype = "Efficacy against variant"
    ) +
    facet_wrap("npi_intensity", labeller = "label_both") +
    theme_bw(base_size = 14) +
    theme(
      strip.text.x = element_text(size = 12, face = "bold"),
      legend.position = "right"
    ) +
    guides(color = guide_legend(ncol = 1, byrow = TRUE))

  # Automate filenames
  filename_png <- paste0(
    "/sensitivity_analyses/efficacy_loss/efficacy_loss_outbreak_size_",
    efficacy_loss_cases_threshold * 100,
    "_perc_pop_isocline_sensitivity_analysis.png"
  )

  filename_eps <- paste0(
    "/sensitivity_analyses/efficacy_loss/efficacy_loss_outbreak_size_",
    efficacy_loss_cases_threshold * 100,
    "_perc_pop_isocline_sensitivity_analysis.eps"
  )


  print(outbreak_size_efficacy_loss_isocline_sensitivity)
  
  # Save the files
  ggsave(outbreak_size_efficacy_loss_isocline_sensitivity,
         filename = filename_png,
         path = git_plot_path,
         width = 9.5,
         units = "in"
  )
  
  ggsave(outbreak_size_efficacy_loss_isocline_sensitivity,
         filename = filename_eps,
         path = git_plot_path,
         width = 9.5,
         units = "in"
  )
}


# Peak prevalence ====
peak_prevalence_thresholds <- c(300, 500, 1000, 10000)

for (i in seq_along(peak_prevalence_thresholds)) {
    
    efficacy_loss_peak_prevalence_threshold <- peak_prevalence_thresholds[i] / target_pop
    
peak_prevalence_efficacy_loss_isocline_df <- efficacy_loss_sensitivity_analysis_df %>% 
    filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time),
           npi_intensity %in% c(0.0, 0.1, 0.2, 0.3)
           )%>%
    filter(peak_prevalence <= efficacy_loss_peak_prevalence_threshold) %>% 
    mutate(
        variant_emergence_day = as_factor(variant_emergence_day),
        vax_efficacy_m = as_factor(vax_efficacy_m)
        ) %>%
    group_by(variant_emergence_day, vax_coverage, npi_intensity, vax_efficacy_m) %>% 
    mutate(min_speed = min(vax_speed)) %>%
    ungroup()
    # ungroup(vax_coverage, npi_intensity, vax_efficacy_m) %>% 
    # mutate(min_coverage = min(vax_coverage), 
    #        max_min_speed = max(min_speed), 
    #        max_speed = 10
    # ) %>% 
    # ungroup() %>% 
    # relocate(c(min_speed, max_min_speed, max_speed), .after = vax_speed)

peak_prevalence_efficacy_loss_isocline <- ggplot(peak_prevalence_efficacy_loss_isocline_df, 
                                   aes(x = vax_coverage, 
                                       y = min_speed, 
                                       color = variant_emergence_day
                                   )) + 
    geom_line(aes(linetype = vax_efficacy_m), 
              size = 1.3, 
              show.legend = TRUE
              ) +
    # geom_segment(data = peak_prevalence_efficacy_loss_isocline_df, 
    #              aes(x = min_coverage, 
    #                  y = max_min_speed,
    #                  xend = min_coverage,
    #                  yend = max_speed,
    #                  linetype = vax_efficacy_m
    #              ),
    #              size = 1.3,
    #              show.legend = TRUE
    # ) +
    scale_x_continuous(labels = percent_format(
        accuracy = 1,
        suffix = ""
    ), breaks = seq(0.10, 1, 0.1)) +
    scale_y_continuous(
        breaks = seq(1, 10, 1),
        labels = seq(1, 10, 1),
        limits = c(1, 10)
    ) +
    scale_color_viridis_d(option = "viridis") +
    labs(
        title = "Sensitivity to vaccine efficacy against variant assumptions",
        subtitle = paste0(
            "Strategies with peak prevalence up to ",
            efficacy_loss_peak_prevalence_threshold * 100,
            "% of total population"
        ),
        x = "Vaccination coverage (%)",
        y = "Vaccination speed",
        color = "Variant emergence day",
        linetype = "Efficacy against variant"
    ) +
    facet_wrap("npi_intensity", labeller = "label_both") +
    theme_bw(base_size = 14) +
    theme(
        strip.text.x = element_text(size = 12, face = "bold"),
        legend.position = "right"
    ) +
    guides(color = guide_legend(ncol = 1, byrow = TRUE)) 

# Automate filenames
filename_png <- paste0(
    "/sensitivity_analyses/efficacy_loss/efficacy_loss_peak_prevalence_",
    efficacy_loss_peak_prevalence_threshold * 100,
    "_perc_pop_isocline_sensitivity_analysis.png"
)

filename_eps <- paste0(
    "/sensitivity_analyses/efficacy_loss/efficacy_loss_peak_prevalence_",
    efficacy_loss_peak_prevalence_threshold * 100,
    "_perc_pop_isocline_sensitivity_analysis.eps"
)


print(peak_prevalence_efficacy_loss_isocline)


# Save the files
ggsave(peak_prevalence_efficacy_loss_isocline,
       filename = filename_png,
       path = git_plot_path,
       width = 9.5,
       units = "in"
)

ggsave(peak_prevalence_efficacy_loss_isocline,
       filename = filename_eps,
       path = git_plot_path,
       width = 9.5,
       units = "in"
)

}

#trial

# efficacy_loss_cases_threshold <- 1E3/target_pop
# 
# outbreak_size_efficacy_loss_isocline_df_trial <- efficacy_loss_sensitivity_analysis_df %>%
#     filter(
#         variant_emergence_day %in% c(1, 61, 121, 151, max_time),
#         npi_intensity %in% c(0)
#     ) %>%
#     filter(total_cases <= efficacy_loss_cases_threshold) %>%
#     mutate(
#         variant_emergence_day = as_factor(variant_emergence_day),
#         vax_efficacy_m = as_factor(vax_efficacy_m)
#     ) %>%
#     group_by(variant_emergence_day, vax_coverage, npi_intensity, vax_efficacy_m) %>%
#     mutate(min_speed = min(vax_speed)) %>%
#     ungroup(vax_coverage, npi_intensity, vax_efficacy_m) %>% 
#     mutate(min_coverage = min(vax_coverage),
#            max_min_speed = max(min_speed),
#            max_speed = 10
#            ) %>% 
#     ungroup() %>% 
#     relocate(c(min_speed, max_min_speed, max_speed), .after = vax_speed)
# 
# 
# 
# trial_plot <- ggplot(
#     outbreak_size_efficacy_loss_isocline_df_trial,
#     aes(
#         x = vax_coverage,
#         y = min_speed,
#         color = variant_emergence_day
#     )
# ) +
#     geom_line(aes(linetype = vax_efficacy_m),
#               size = 1.3,
#               show.legend = TRUE
#     ) +
#     geom_segment(data = outbreak_size_efficacy_loss_isocline_df_trial, 
#               aes(x = min_coverage, 
#                   y = max_min_speed,
#                   xend = min_coverage,
#                   yend = max_speed,
#                   linetype = vax_efficacy_m
#                   ),
#               size = 1.3,
#               show.legend = TRUE
#     ) +
#     scale_x_continuous(labels = percent_format(
#         accuracy = 1,
#         suffix = ""
#     ), breaks = seq(0.10, 1, 0.1)) +
#     scale_y_continuous(
#         breaks = seq(1, 10, 1),
#         labels = seq(1, 10, 1),
#         limits = c(1, 10)
#     ) +
#     scale_color_viridis_d(option = "viridis") +
#     labs(
#         title = "Sensitivity to vaccine efficacy against variant assumptions",
#         subtitle = paste0(
#             "Strategies with cumulative cases up to ",
#             efficacy_loss_cases_threshold * 100,
#             "% of total population"
#         ),
#         x = "Vaccination coverage (%)",
#         y = "Vaccination speed",
#         color = "Variant emergence day",
#         linetype = "Efficacy against variant"
#     ) +
#     facet_wrap("npi_intensity", labeller = "label_both") +
#     theme_bw(base_size = 14) +
#     theme(
#         strip.text.x = element_text(size = 12, face = "bold"),
#         legend.position = "right"
#     ) +
#     guides(color = guide_legend(ncol = 1, byrow = TRUE))
# 
# print(trial_plot)
