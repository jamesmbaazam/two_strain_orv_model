# Sensitivity analyses ----
# Source the global inputs
source('./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/global_inputs/global_params.R')

# Dynamics parameters ----
R0_w <- 2.0 # R0 (wild type)
R0_m_original <- R0_w * 1.3 # R0 (variant) 30% more infectious
IP_w <- 14 # infectious period (wild type)
IP_m <- 14 # infectious period (variant)


# Cross protection, cp (everything else remains as original) ====
#' The assumption is that cross protection, cp, is symmetric so that 
#' sigma_w = sigma_m = 0.85

dynamics_params_cp_sensitivity <- data.frame(
    R0w = R0_w,
    R0m = R0_m_original,
    beta_w = R0_w / IP_w,
    beta_m = R0_m_original / IP_m,
    gamma_w = 1 / IP_w,
    gamma_m = 1 / IP_m,
    sigma_w = 0.85, # cross-protection provided by wild-type
    sigma_m = 0.85, # cross-protection provided by variant
    vax_efficacy_w = 1, # perfect vaccine efficacy against wild-type
    vax_efficacy_m = 1 # perfect vaccine efficacy against variant
)

