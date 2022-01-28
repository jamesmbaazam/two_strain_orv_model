rm(list = ls()) #Clear the environment

# Cross protection ----
#1. Run the cross protection sensitivity analyses and save output
source('./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/cross_protection/sims_cross_protection_sensitivity_parallel.R')

# Variant transmissibility ----
#2. Run the variant R0 sensitivity analyses and save output
source('./scripts/two_strain_model/vax_escape_analyses/sensitivity_analyses/variant_R0/sims_vR0_sensitivity_parallel.R')
