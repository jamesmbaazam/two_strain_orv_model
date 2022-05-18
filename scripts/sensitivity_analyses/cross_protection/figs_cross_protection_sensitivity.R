# Packages ----
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
thesis_plot_path <- "C:/Users/James Azam/Dropbox/My Academic Repository/_Degrees/PhD/__PhD_Thesis/_Thesis_LaTeX/figs/chapter_5/"

# Change how large numbers are printed
options(scipen = 10000)


# Load the model output

# Original analysis
#' The baseline analysis assumes perfect vaccine efficacy against the wild-type
#' and variant and perfect cross protection between the strains. The variant is
#' assumed to be 30% more transmissible than the wild-type (R0w = 1.3R0m)
#'
#'
# Cross protection sensitivity analysis results
cp_sensitivity_analysis_summaries <- readRDS("./model_output/sensitivity_analyses/cross_protection/cp_sensitivity_analysis_summaries.rds")


# Remove redundant columns
cp_sensitivity_analysis_df <- cp_sensitivity_analysis_summaries %>%
    select(-c(
        vax_rate,
        npi_duration,
        starts_with("vax_efficacy"),
        #  cross_protection_m,
        R0m,
        total_vaccinated
    ))


#' Determine the subset of scenarios that meet the threshold outbreak size
#' and the minimum speed required
total_cases_thresholds <- c(1E3, 5E5, 5E6, 10E6)


for (cases_threhold in seq_along(total_cases_thresholds)) {
    
    cp_total_cases_threshold <- total_cases_thresholds[cases_threhold] / target_pop
    
    # Determine the isocline data
    outbreak_size_cp_isocline_df <- cp_sensitivity_analysis_df %>%
        filter(
            cross_protection_w %in% c(0.5, 1),
            variant_emergence_day %in% c(1, 61, 121, 151, max_time),
            npi_intensity %in% c(0.00, 0.10, 0.20, 0.30)
        ) %>%
        filter(total_cases <= cp_total_cases_threshold) %>% # keep scenarios with this level of total cases
        mutate(
            variant_emergence_day = as_factor(variant_emergence_day),
            cross_protection_w = as_factor(cross_protection_w)
        ) %>%
        group_by(variant_emergence_day, vax_coverage, npi_intensity, cross_protection_w) %>%
        mutate(min_speed = min(vax_speed)) %>%
        ungroup()
    
    # Isoclines ----
    # Outbreak size ====
    outbreak_size_cp_isocline_sensitivity <- ggplot(
        outbreak_size_cp_isocline_df,
        aes(
            x = vax_coverage,
            y = min_speed,
            color = variant_emergence_day
        )
    ) +
        geom_line(aes(linetype = cross_protection_w),
                  size = 1.3,
                  show.legend = TRUE
        ) +
        scale_x_continuous(
            labels = percent_format(
                accuracy = 1,
                suffix = ""
            ),
            breaks = seq(0.10, 1, 0.1)
        ) +
        scale_y_continuous(
            breaks = seq(1, 10, 1),
            labels = seq(1, 10, 1),
            limits = c(1, 10)
        ) +
        scale_color_viridis_d(option = "viridis") +
        labs(
            title = "Sensitivity to cross protection assumptions",
            subtitle = paste0(
                "Strategies with cumulative cases up to ",
                cp_total_cases_threshold * 100,
                "% of total population"
            ),
            x = "Vaccination coverage (%)",
            y = "Vaccination speed",
            color = "Variant emergence day",
            linetype = "Cross protection"
        ) +
        facet_wrap("npi_intensity", labeller = "label_both") +
        theme(
            strip.text.x = element_text(size = 12, face = "bold"),
            legend.position = "right"
        ) +
        theme_bw(base_size = 14) +
        guides(color = guide_legend(ncol = 1, byrow = TRUE))
    
    
    # Automate filenames
    filename_png <- paste0(
        "/sensitivity_analyses/cross_protection/cp_outbreak_size_",
        cp_total_cases_threshold * 100,
        "_perc_pop_isocline_sensitivity_analysis.png"
    )
    
    filename_eps <- paste0(
        "/sensitivity_analyses/cross_protection/cp_outbreak_size_",
        cp_total_cases_threshold * 100,
        "_perc_pop_isocline_sensitivity_analysis.eps"
    )
    
    print(outbreak_size_cp_isocline_sensitivity)
    
    # Save the files
    ggsave(outbreak_size_cp_isocline_sensitivity,
           filename = filename_png,
           path = git_plot_path,
           width = 9.5,
           height = 6.5,
           units = "in"
    )
    
    ggsave(outbreak_size_cp_isocline_sensitivity,
           filename = filename_eps,
           path = git_plot_path,
           width = 9.5,
           height = 6.5,
           units = "in"
    )
}



# Peak prevalence ====
peak_prevalence_thresholds <- c(300, 500, 1000, 10000)

for (i in seq_along(peak_prevalence_thresholds)) {
    
    cp_peak_prevalence_threshold <- peak_prevalence_thresholds[i] / target_pop
    
    #Determine the isocline data
    peak_prevalence_cp_isocline_df <- cp_sensitivity_analysis_df %>%
        filter(
            cross_protection_w %in% c(0.5, 1),
            variant_emergence_day %in% c(1, 61, 121, 151, max_time),
            npi_intensity %in% c(0.0, 0.1, 0.2, 0.3)
        ) %>%
        filter(peak_prevalence <= cp_peak_prevalence_threshold) %>% # keep scenarios with this level of total cases
        mutate(
            variant_emergence_day = as_factor(variant_emergence_day),
            cross_protection_w = as_factor(cross_protection_w)
        ) %>%
        group_by(variant_emergence_day, vax_coverage, npi_intensity, cross_protection_w) %>%
        mutate(min_speed = min(vax_speed)) %>%
        ungroup()
    
    # Plot the isoclines
    peak_prevalence_cp_isocline <- ggplot(
        peak_prevalence_cp_isocline_df,
        aes(
            x = vax_coverage,
            y = min_speed,
            color = variant_emergence_day
        )
    ) +
        geom_line(aes(linetype = cross_protection_w),
                  size = 1.3,
                  show.legend = TRUE
        ) +
        scale_x_continuous(
            labels = percent_format(
                accuracy = 1,
                suffix = ""
            ),
            breaks = seq(0.10, 1, 0.1)
        ) +
        scale_y_continuous(
            breaks = seq(1, 10, 1),
            labels = seq(1, 10, 1),
            limits = c(1, 10)
        ) +
        scale_color_viridis_d(option = "viridis") +
        labs(
            title = "Sensitivity to cross protection assumptions",
            subtitle = paste0(
                "Strategies with peak prevalence up to ",
                cp_peak_prevalence_threshold * 100,
                "% of total population"
            ),
            x = "Vaccination coverage (%)",
            y = "Vaccination speed",
            color = "Variant emergence day",
            linetype = "Cross protection"
        ) +
        facet_wrap("npi_intensity", labeller = "label_both") +
        theme_bw(base_size = 14) +
        theme(
            strip.text.x = element_text(size = 12, face = "bold"),
            legend.position = "right"
        ) +
        guides(color = guide_legend(ncol = 1, byrow = TRUE))
    
    # Automate filenames
    peak_prev_plot_filename_png <- paste0(
        "/sensitivity_analyses/cross_protection/cp_peak_prev_",
        cp_peak_prevalence_threshold * 100,
        "_perc_pop_isocline_sensitivity_analysis.png"
    )
    
    peak_prev_plot_filename_eps <- paste0(
        "/sensitivity_analyses/cross_protection/cp_peak_prev_",
        cp_peak_prevalence_threshold * 100,
        "_perc_pop_isocline_sensitivity_analysis.eps"
    )
    
    
    print(peak_prevalence_cp_isocline)
    
    
    # Save the files
    ggsave(peak_prevalence_cp_isocline,
           filename = peak_prev_plot_filename_png,
           path = git_plot_path,
           width = 9.5,
           units = "in"
    )
    
    ggsave(peak_prevalence_cp_isocline,
           filename = peak_prev_plot_filename_eps,
           path = git_plot_path,
           width = 9.5,
           units = "in"
    )
}
