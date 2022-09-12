###################################################
#Holisoils Multimodel ensemble version 0
###################################################
#This function launches an ensemble of models 
#that are either coded in SoilR, or sourced from .R files that must come with this script
#The models currently included in the ensemble are:

# Century
#' Parton, W.J, D.S. Schimel, C.V. Cole, and D.S. Ojima. 1987.
#' Analysis of factors controlling soil organic matter levels in Great Plain
#' grasslands. Soil Science Society of America Journal 51: 1173--1179.

#RothC
#' Jenkinson, D. S., S. P. S. Andrew, J. M. Lynch, M. J. Goss, and
#' P. B. Tinker. 1990. The Turnover of Organic Carbon and Nitrogen in Soil.
#' Philosophical Transactions: Biological Sciences 329:361-368.

#ICBM 
#' Andren, O. and T. Katterer. 1997. ICBM: The Introductory Carbon
#' Balance Model for Exploration of Soil Carbon Balances. Ecological
#' Applications 7:1226-1236.

#Yasso07
#' Tuomi, M., Thum, T., Jarvinen, H., Fronzek, S., Berg, B.,
#' Harmon, M., Trofymow, J., Sevanto, S., and Liski, J. (2009). Leaf litter
#' decomposition-estimates of global variability based on Yasso07 model.
#' Ecological Modelling, 220:3362 - 3371.
#' 

#AMG
#'Clivot, Hugues, Jean-Christophe Mouny, Annie Duparque,
#' Jean-Louis Dinh, Pascal Denoroy, Sabine Houot, Françoise Vertès,
#'  et al. “Modeling Soil Organic Carbon Evolution in Long-Term 
#'  Arable Experiments with AMG Model.” 
#'  Environmental Modelling & Software 118 (August 2019): 99–113.
#'   https://doi.org/10.1016/j.envsoft.2019.04.004.


#SG
#'Hashimoto, S. et al. (2011) 
#'Simple models for soil CO2, CH4, and N2O fluxes calibrated 
#'using a Bayesian approach and multi-site data. Ecological Modelling, 
#'222: 1283-1292. http://dx.doi.org/10.1016/j.ecolmodel.2011.01.013


#Other references
#' Sierra, C.A., M. Mueller, S.E. Trumbore. 2012. Models of soil organic matter
#' decomposition: the SoilR package version 1.0. Geoscientific Model
#' Development 5, 1045-1060.


########################################################################
#PREDIFINED INPUT DATA & UNITS
########################################################################

#If not specified otherwise, all soil variables should refer to the the top 0-30 cm
#         --Plot--
#'@param plot_figures a boolean (TRUE/FALSE) deciding wether to plot the default figures or not
#         --Time--
#'@param simulation_length               A scalar indicating the forward simulation lenght: [years]
#'@param spinup_length                   A scalar indicating the lenght of the spinup run:  [years]
# computation_time_step_fwd       The time step of the forward simulation:           [years]
#computation_time_step_spinup:   The time step of the spinup run:                   [years]

#'@param start_date_simulations          A `Date` object in the format "YYYY-MM-DD"" [e.g. as.Date("2021-01-01")]
#         --Meteo
#'@param temperature                     A `data.frame` object that has as first column the dates of measurments in the format "YYYY-MM" and as second column the monthly average temperatures [˚C]
#'@param precipitation                   A `data.frame` object that has as first column the dates of measurments in the format "YYYY-MM" and as second column the monthly cumulative precipitations [mm/month]
#'@param potential_evapotranspiration    A `data.frame` object that has as first column the dates of measurments in the format "YYYY-MM" and as second column the monthly potentia evapotranspiration [mm/month]
#? soil_moisture

#         --Soil--
#'@param SOC_0                           A scalar indicating the level of SOC stock at the beginning of the simulation                                             [Mg/ha]
#'@param C_input_spinup                  A scalar indicating the level of average annual C input during the spinup run                                             [Mg/ha/yr]
#'@param C_input_fwd                     A vector (or scalar if information is not available) indicating the level of annual C input during the forward simulation [Mg/ha/yr]
#'@param clay_p                          A scalar indicating the percentage concentration of clay       [%]
#'@param silt_p                          A scalar indicating the percentage concentration of silt       [%]
#'@param soil_thickness                  A scalar indicating the thikness of the organic layer topsoil  [cm] (default is 25)

#         --Litter--
#'@param lignin_to_nitrogen              A scalar indicating Lignin:Nitrogen ratio of the litter input:          [unitless] (default is 0.5) #CHANGE!
#'@param structural_in_lignin            A scalar indicating the fraction of structural material in the lignin:  [unitless] (default is 0.1) #CHANGE!


#       --Decomposition rate parameters
#'@param decomposition_param_RothC       A vector of 5 elements providing the decomposition rate parameters of the RothC C pools:  [1/yr] (default is: c(k.DPM = 10, k.RPM = 0.3, k.BIO = 0.66, k.HUM = 0.02, k.IOM = 0))
#'@param decomposition_param_ICBM        A vector of 3 elements providing: the 2 decomposition rate parameters of the ICBM C pools [1/yr] and the humification coefficient [unitless] (default is: c(k1=0.8,k2=0.00605,h=0.13))
#'@param decomposition_param_Century     A vector of 7 elements providing the decomposition rate parameters of the Century C pools [1/week] (default is c(STR.surface = 0.076, MET.surface = 0.28, STR.belowground = 0.094, MET.belowground = 0.35, ACT = 0.14, SLW = 0.0038, PAS = 0.00013))

#(SG)
#         Bulk density at 0-5 cm depth:          [Mg m-3]
#         Carbon:Nitrogen ratio at 0-5 cm depth: [unitless]

#(YASSO)
#         WoodySize - size of woody litter (for non-woody litter this is 0) [cm]


Call_MULTIMODEL<-function(plot_figures,simulation_length, spinup_length,
                          #computation_time_step_fwd=1, computation_time_step_spinup=1, tolti, perche ogni modello ha il suo gia predefinito
                          start_date_simulations,
                          temperature, precipitation, potential_evapotranspiration,
                          SOC_0,C_input_spinup,C_input_fwd,clay_p,silt_p,carbonat_p, soil_thickness,
                          lignin_to_nitrogen,structural_in_lignin,woodylittersize,
                          CN_Ratio, Bulk_Density, WFPS, vswc,CH4_Conc,
                          historical_land_use,
                          decomposition_param_RothC,
                          decomposition_param_ICBM,
                          decomposition_param_Century,
                          decomposition_param_Yasso07,
                          decomposition_param_AMG){
  

  
  
  #---------------------------------------------
  #CREATE INTEGRALS TIME STEPS
  #---------------------------------------------
  #Define time step spinup run [in years]
  computation_time_step_spinup = 1
  #Define time step forward run [in years]
  #computation_time_step_fwd = 1/12
  
  #Create a vector of points in time where the spinup run solution is sought
  t_spinup = seq(1,spinup_length,by=computation_time_step_spinup)
  #Create a vector of points in time where the forward solution is sought
  t_fwd = seq(1,simulation_length,by=computation_time_step_fwd)
 
   
  #CALL ROTHC
  Roth_C_fwd<-Call_RothC(plot_figures=plot_figures,
                       temperature=temperature, precipitation=precipitation, potential_evapotranspiration=potential_evapotranspiration,
                       SOC_0=SOC_0,C_input_spinup=C_input_spinup,C_input_fwd=C_input_fwd,clay_p=clay_p,soil_thickness=soil_thickness,
                       decomposition_param_RothC=decomposition_param_RothC,t_spinup=t_spinup,t_fwd=t_fwd)
  print(Roth_C_fwd)
  
  #CALL ICBM
  ICBM_C_fwd<-Call_ICBM(plot_figures=plot_figures,
                         temperature=temperature, precipitation=precipitation, potential_evapotranspiration=potential_evapotranspiration,
                         SOC_0=SOC_0,C_input_spinup=C_input_spinup,C_input_fwd=C_input_fwd,clay_p=clay_p,soil_thickness=soil_thickness,
                         decomposition_param_ICBM=decomposition_param_ICBM,t_spinup=t_spinup,t_fwd=t_fwd)
  print(ICBM_C_fwd)
  
  #CALL CENTURY
  Century_C_fwd<-Call_Century(plot_figures=plot_figures,
                              temperature=temperature, precipitation=precipitation, potential_evapotranspiration=potential_evapotranspiration,
                              SOC_0=SOC_0,C_input_spinup=C_input_spinup,C_input_fwd=C_input_fwd,clay_p=clay_p,silt_p=silt_p,
                         lignin_to_nitrogen=lignin_to_nitrogen, structural_in_lignin=structural_in_lignin,
                         decomposition_param_Century=decomposition_param_Century,
                         t_spinup=t_spinup,t_fwd=t_fwd)

  print(Century_C_fwd)
  
  #CALL YASSO07
  Yasso07_C_fwd<-Call_Yasso07(plot_figures=plot_figures,
               temperature=temperature, precipitation=precipitation, woodylittersize=woodylittersize,
               #simulation_length=simulation_length, #remove once Yasso function for temp and moist implemented
               SOC_0=SOC_0,C_input_spinup=C_input_spinup,C_input_fwd=C_input_fwd,
               decomposition_param_Yasso07=decomposition_param_Yasso07,
               t_spinup=t_spinup,t_fwd=t_fwd)

  print(Yasso07_C_fwd)
  
  #CALL AMG
  AMG_C_fwd<-Call_AMG(plot_figures=plot_figures,
                      temperature=temperature, precipitation=precipitation, potential_evapotranspiration=potential_evapotranspiration,
                      SOC_0=SOC_0,C_input_fwd=C_input_fwd,clay_p=clay_p,carbonat_p=carbonat_p,
                      decomposition_param_AMG=decomposition_param_AMG,
                      historical_land_use=historical_land_use,
                      t_fwd=t_fwd)
  
  print(AMG_C_fwd)
  
  #CALL SG
  SG_C_fwd<-Call_SG(CN_Ratio=CN_Ratio, Bulk_Density=Bulk_Density, WFPS=WFPS, vswc=vswc, Soil_Temperature=temperature, CH4_Conc=CH4_Conc)
  #Select only sought solutions
  SG_C_fwd<-SG_C_fwd[c(1:length(t_fwd)),]
  SG_C_fwd[,1]=SG_C_fwd[,1]*(1e-12*1e4*3.171e8*computation_time_step_fwd) #convert CO2: from ugC/m2/sec to tC/ha/time
  SG_C_fwd[,c(2,3)]=SG_C_fwd[,c(2,3)]*(1e-12*1e4*8760*computation_time_step_fwd) #convert CO2: from ug/m2/hour to tC/ha/time
  print(SG_C_fwd)
  
  #Plot SG fluxes
  if(computation_time_step_fwd==1/12){
    time_legend = "month"
  }else if(computation_time_step_fwd==1){
    time_legend = "yr"
  }else if(computation_time_step_fwd==1/365){
    time_legend = "day"
  }else{
    time_legend = "CHECK UNITS"
  }
  if(plot_figures==TRUE){
    legend_SG<-c("CH4", "N20")
    #plot the pools
    matplot(t_fwd, SG_C_fwd[,c(2,3)], type="l", lty=1, col=1:2,
            xlab="Time (years)", ylab=paste0("Fluxes (tC/ha/",time_legend,")"),main="SG fluxes")
    legend("topleft", legend_SG,
           lty=1, col=1:2, bty="n")
  }
  
  #PLOT MULTIMODEL
  
  if(plot_figures==TRUE){
    #plot total C stocks
    plot(t_fwd, rowSums(Roth_C_fwd[[2]]), type="l", lty=1, col="red",
            xlab="Time (years)", ylab="C stocks (MgC/ha)", ylim=c(50,70),main="Multi-model")
    lines(t_fwd,rowSums(ICBM_C_fwd[[2]]),type="l", lty=1, col="black")
    lines(t_fwd,rowSums(Century_C_fwd[[2]]),type="l", lty=1, col="green")
    lines(t_fwd,rowSums(Yasso07_C_fwd[[2]]),type="l", lty=1, col="pink")
    lines(t_fwd,rowSums(AMG_C_fwd[[2]]),type="l", lty=1, col="blue")
    legend("topleft", c("RothC", "ICBM","Century","Yasso07","AMG"),
           lty=1, col=c("red","black","green","pink","blue"),ncol=2, bty="n")
    
    #plot CO2 fluxes

    plot(t_fwd, rowSums(Roth_C_fwd[[1]]), type="l", lty=1, col="red",
         xlab="Time (years)", ylab=paste0("CO2 fluxes (MgC/ha/",time_legend,")"), ylim=c(-1,10),main="Multi-model")
    lines(t_fwd,rowSums(ICBM_C_fwd[[1]]),type="l", lty=1, col="black")
    lines(t_fwd,rowSums(Century_C_fwd[[1]]),type="l", lty=1, col="green")
    lines(t_fwd,rowSums(Yasso07_C_fwd[[1]]),type="l", lty=1, col="pink")
    lines(t_fwd,rowSums(AMG_C_fwd[[1]]),type="l", lty=1, col="blue")
    lines(t_fwd, SG_C_fwd[,1],type="l", lty=1, col="violet")
    legend("topleft", c("RothC", "ICBM","Century","Yasso07","AMG","SG"),
           lty=1, col=c("red","black","green","pink","blue","violet"), ncol=2,bty="n")
    
  }

  return(list(Roth_C_fwd,ICBM_C_fwd,Century_C_fwd))
  
}



