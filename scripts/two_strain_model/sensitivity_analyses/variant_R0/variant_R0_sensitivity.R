#Packages ----
library(tidyverse)
library(patchwork)
library(scales)
library(mdthemes)
library(tidyverse)


#helper scripts
source('./scripts/two_strain_model/sensitivity_analyses/sim_config_global_params_sensitivity.R')
source('./scripts/two_strain_model/two_strain_model.R')
source('./scripts/two_strain_model/sim_config_emergence_risk_adjusted.R')

#Load the model output

#variant R0 30% more than wildtype
variant_R0_30_percent <- readRDS('./model_output/sensitivity_analyses/variant_R0/model_dynamics_R0m_30_percent.rds')

#variant R0 60% more than wildtype
variant_R0_60_percent <- readRDS('./model_output/sensitivity_analyses/variant_R0/model_dynamics_R0m_60_percent.rds')


#Combine the files
variant_R0_sensitivity <- bind_rows(variant_R0_30_percent, variant_R0_60_percent)



#Rescale the population proportions to total sizes
variant_R0_sensitivity_isocline_df <- variant_R0_sensitivity %>%
    mutate(across(.cols = c(total_cases, peak_cases, total_vaccinated), 
                  .fns = ~ .x*target_pop #rescale the population proportions to total sizes
    )
    ) %>% 
    select(-c(npi_duration, total_vaccinated)) %>% 
    filter(variant_emergence_day %in% c(1, 61, 121, 151, max_time), 
           npi_intensity %in% c(0.0, 0.1, 0.2, 0.3), 
           total_cases <= 1000
           ) %>% 
    mutate(variant_emergence_day = as_factor(variant_emergence_day)) %>% 
    group_by(variant_emergence_day, vax_coverage, npi_intensity, R0m) %>% 
    mutate(min_speed = min(vax_speed), R0m = as.factor(R0m)) %>% 
    ungroup()





# Isoclines ----
# Outbreak size ====
outbreak_size_isocline <- ggplot(variant_R0_sensitivity_isocline_df, 
       aes(x = vax_coverage, 
           y = min_speed, 
           color = variant_emergence_day,
           linetype = R0m
           )
       ) + 
    geom_line(size = 1, show.legend = TRUE) + 
    scale_x_continuous(labels = percent_format(), breaks = seq(0.30, 1, 0.1)) +
    scale_y_continuous(breaks = seq(1, 10, 1), labels = seq(1, 10, 1)) +
    scale_color_viridis_d(option = 'viridis') +
    labs(#title = paste('Cumulative cases threshold <= 1000 at various NPI intensity levels'), 
        x = 'Vaccination coverage', 
        y = 'Vaccination speed', 
        color = 'Variant emergence day',
        linetype = 'Variant R0'
        ) +
    facet_wrap('npi_intensity', labeller = 'label_both') +
    theme_bw(base_size = 14) +
    theme(strip.text.x = element_text(size = 12, face = 'bold'), legend.position = 'bottom') 


print(outbreak_size_isocline)
