#Shortcuts to directories
MODEL_SIM_FUNCTIONS := ./scripts/model_sim_functions
INPUTS_DIR := ./data/inputs
SENSITIVITY_ANALYSIS_DIR := ./scripts/sensitivity_analysis
MAIN_ANALYSIS_DIR := ./scripts/main_analysis


#Convert the scripts into Rdata objects for use in analyses
$(INPUTS_DIR)/two_strain_model_functions.RData: $(MODEL_SIM_FUNCTIONS)/two_strain_model_functions.R
	Rscript $^ $@

$(INPUTS_DIR)/simulation_functions.RData: $(MODEL_SIM_FUNCTIONS)/simulation_functions.R
	Rscript $^ $@

$(INPUTS_DIR)/config_global_params.RData: $(MAIN_ANALYSIS_DIR)/sim_config/config_global_params.R
	Rscript $^ $@

$(INPUTS_DIR)/config_intervention_params.RData: $(MAIN_ANALYSIS_DIR)/sim_config/config_intervention_params.R $(INPUTS_DIR)/config_global_params.RData $(INPUTS_DIR)/two_strain_model_functions.RData $(INPUTS_DIR)/simulation_functions.RData
	Rscript $^ $@

