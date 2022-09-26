#######################################################################
#This is the script that needs to be launched to run the simulations
#######################################################################

#---------------------------------
#To install the required packages (i.e., libraries) in R, the following command can be used: install.packages(“PACKAGE_REQUIRED”)
#Example to install SoilR:
#Run the following command
#install.packages("SoilR")
#---------------------------------
#Once the libraries installed, the user can upload them
#Libraries that need to be installed and uploaded
library(SoilR)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(readxl)
library(matrixStats)

#Define the path to the folders where the scripts are stored
PATH_functions <- "/PATH/TO/FOLDERS/" #Choose the path where you stored the documents
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
#Yasso07 SoilR version
#source(paste0(PATH_functions,"Holisoils_Yasso07_v0_nonFIX.R"))
#source(paste0(PATH_functions,"Yasso07_environmental_functions.R"))
#Yasso07 Boris version
source(paste0(PATH_functions,"Holisoils_Yasso07_v0.R"))
source(paste0(PATH_functions,"Yasso07Model_fixed.R"))

#---------------------
#Other parameters that need to be defined (although it doesn't have to be the end-user)
#Define spinup duration [years]
spinup_length=1000
#Define forward simulation length [years]
#as number of data/time_step
simulation_length=length(Temp_month$Date)/12
#Choose how many plots to draw in the same figure (num of rows, num of columns)
par(mfrow=c(2,1))
#---------------------

##############################
#Call the multi-model ensemble
###############################
#This will launch the function defined in "Holisoils_multimodel_v0.R"
test_mm <-Call_MULTIMODEL(plot_figures=plot_figures,
                          simulation_length=simulation_length, spinup_length=spinup_length,
                          start_date_simulations=start_date_simulations_site,
                          temperature=Temp_month, precipitation=Precip_month, potential_evapotranspiration=Potevap_month,
                          SOC_0=SOC_initial,
                          C_input_spinup=Cinput_spinup,C_input_fwd=Cinput_fwd,
                          clay_p=clay_site,silt_p=silt_site,carbonat_p=CaCO3_site, soil_thickness=soil_OM_thick,
                          lignin_to_nitrogen=lignin_to_nitrogen_ratio,structural_in_lignin=structural_in_lignin_ratio,
                          woodylittersize=woodylittersize_scalar,
                          CN_Ratio=CN_site, Bulk_Density=BD_site, WFPS=water_filled_pore_space, vswc=volumetric_soil_water_cont,CH4_Conc=CH4_conc_site,
                          historical_land_use=historical_LU,
                          decomposition_param_RothC=ksRothC,
                          decomposition_param_ICBM=param_ICBM,
                          decomposition_param_Century=ksCent,
                          decomposition_param_Yasso07=paramYasso07,
                          decomposition_param_AMG = ksAMG)


