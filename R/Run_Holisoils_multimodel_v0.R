
#Libraries that need to be installed and uploaded
library(SoilR)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(readxl)
library(matrixStats)

PATH_functions <- "/PATH/TO/FOLDERS" #Choose the path where you stored the documents
#---------------------------
#Call required functions
source(paste0(PATH_functions,"Holisoils_multimodel_v0.R"))
source(paste0(PATH_functions,"Holisoils_RothC_v0.R"))
source(paste0(PATH_functions,"Holisoils_ICBM_v0.R"))
source(paste0(PATH_functions,"Holisoils_Century_v0.R"))
source(paste0(PATH_functions,"Holisoils_AMG_v0.R"))
source(paste0(PATH_functions,"AMG_environmental_functions.R"))
source(paste0(PATH_functions,"Holisoils_SG_v0.R"))
source(paste0(PATH_functions,"Forcing_data_test_run_v0.R"))
#Yasso SoilR version
#source(paste0(PATH_functions,"Holisoils_Yasso07_v0_nonFIX.R"))
#source(paste0(PATH_functions,"Yasso07_environmental_functions.R"))
#Yasso Boris version
source(paste0(PATH_functions,"Holisoils_Yasso07_v0.R"))
source(paste0(PATH_functions,"yasso07Model.soilr.fi_function.r"))

#---------------------------
#USER DEFINED PARAMETERS
#Decide wether to plot default figures [TRUE] or not [FALSE]
plot_figures=TRUE

#Define time step forward run [in years]
computation_time_step_fwd = 1/12

#---------------------
#Other parameters that need to be defined
#Define spinup duration [years]
spinup_length=1000
#Define forward simulation length [years]
#as number of data/time_step
simulation_length=length(Temp_month$Date)*computation_time_step_fwd

par(mfrow=c(2,1))
test_mm <-Call_MULTIMODEL(plot_figures=plot_figures,simulation_length=simulation_length, spinup_length=spinup_length,
                          #computation_time_step_fwd=1, computation_time_step_spinup=1, tolti, perche ogni modello ha il suo gia predefinito
                          start_date_simulations=start_date_simulations_site,
                          temperature=Temp_month, precipitation=Precip_month, potential_evapotranspiration=Potevap_month,
                          SOC_0=SOC_initial,C_input_spinup=Cinput_spinup,C_input_fwd=Cinput_fwd,clay_p=clay_site,silt_p=silt_site,carbonat_p=CaCO3_site, soil_thickness=soil_OM_thick,
                          lignin_to_nitrogen=lignin_to_nitrogen_ratio,structural_in_lignin=structural_in_lignin_ratio,woodylittersize=woodylittersize_scalar,
                          CN_Ratio=CN_site, Bulk_Density=BD_site, WFPS=water_filled_pore_space, vswc=volumetric_soil_water_cont,CH4_Conc=CH4_conc_site,
                          historical_land_use=historical_LU,
                          decomposition_param_RothC=ksRothC,
                          decomposition_param_ICBM=param_ICBM,
                          decomposition_param_Century=ksCent,
                          decomposition_param_Yasso07=paramYasso07,
                          decomposition_param_AMG = ksAMG)


