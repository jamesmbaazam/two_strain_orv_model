#Packages ----
library(tidyverse)

#Read the model output
variant_R0_30_percent <- readRDS('./model_output/sensitivity_analyses/variant_R0/model_dynamics_R0m_30_percent.rds')


variant_R0_60_percent <- readRDS('./model_output/sensitivity_analyses/variant_R0/model_dynamics_R0m_60_percent.rds')



variant_R0_sensitivity <- bind_rows(variant_R0_30_percent, variant_R0_60_percent)