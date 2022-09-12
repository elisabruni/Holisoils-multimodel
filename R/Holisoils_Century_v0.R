###################################################
#Century model
###################################################
#This function launches the SoilR version of Century
#First, a spin-up run is performed
#Then, the estimated C pool fractions are used to rescale total SOC stocks with measured SOC stocks
#Finally, a forward run is perfomed
#The function provides SOC stock and CO2 fluxes for the chosen simulation lenght

# Century
#' Parton, W.J, D.S. Schimel, C.V. Cole, and D.S. Ojima. 1987.
#' Analysis of factors controlling soil organic matter levels in Great Plain
#' grasslands. Soil Science Society of America Journal 51: 1173--1179.
#' 


Call_Century<-function(plot_figures,
                     temperature, precipitation,potential_evapotranspiration,
                     SOC_0,C_input_spinup,C_input_fwd,clay_p,silt_p,
                     lignin_to_nitrogen, structural_in_lignin,
                     decomposition_param_Century,
                     t_spinup,t_fwd){
  
  
  #Convert input variables Century (monthly)
  ts_spinup_Cent<-t_spinup*12 #[months]
  ts_fwd_Cent<-t_fwd*12 #[months]
  InCent_spinup<-C_input_spinup/12 #MgC/ha/yr => MgC/ha/month
  InCent_fwd<-C_input_fwd/12 #MgC/ha/yr => MgC/ha/month
  clayCent<-clay_p/100. #% to proportion (i.e. [mass/mass])
  siltCent<-silt_p/100. #% to proportion
  
  legend_Cent<-c("Surface structural", "Surface metabolic", "Belowground structural",
    "Belowground metabolic", "Active SOM", "Slow SOM", "Passive SOM")
  
  #CHANGE if possible with real data
  #To calculate aboveground and belowground inputs
  APPT = mean(precipitation[,2])
  Pmax=-40+7.7*APPT # Max aboveground production
  Rmax=100+7.0*APPT # Max belowground production
  abvgIn=Pmax/(Pmax+Rmax) #proportion of aboveground (?)
  blgIn=Rmax/(Pmax+Rmax)  #proportion of belowground (?)

  #--
  #Rate modifiers
  fT=fT.Century1(temperature$Temp) #Temperature effects per month
  fW=fW.Century(PPT = (precipitation$Precip), PET =(potential_evapotranspiration$Potevap)) #Moisture effects per month
  
  #Rate modifier spinup
  ##############################################################################
  #MODIFY FOLLOWING LINE FOR EFFICIENCY: xi.frame_spinup_Century=mean(fT*fW)
  ##############################################################################
  xi.frame_spinup_Century=data.frame(ts_spinup_Cent,rep(mean(fT*fW),length.out=length(ts_spinup_Cent)))
  
  #Rate modifier forward
  xi.frame_Century=data.frame(ts_fwd_Cent,rep(fT*fW,length.out=length(ts_fwd_Cent)))
  
  #Run spinup
  Century_spinup = CenturyModel(t=ts_spinup_Cent,ks=decomposition_param_Century,C0=rep(0,7),surfaceIn=InCent_spinup*abvgIn, soilIn = InCent_spinup*blgIn, LN=lignin_to_nitrogen, Ls=structural_in_lignin,
                                clay=clayCent,silt=siltCent,xi=xi.frame_spinup_Century,
                                solver=deSolve.lsoda.wrapper)
  
  #Get SOC spinup
  CCentury_spinup=getC(Century_spinup)

  
  #plot the pools
  if(plot_figures==TRUE){
    matplot(ts_spinup_Cent, CCentury_spinup, type="l", lty=1, col=1:7,
            xlab="Time (months)", ylab="C stocks (MgC/ha)",main="Century spinup")
    legend("topleft", legend_Cent,
           lty=1, col=1:7, bty="n")
  }
  #####
  #2)## Relaxation
  #####
  #C after spinup in each pool
  CCentury_pools_spinup = tail(CCentury_spinup,1)
  #Total C after spinup
  CCentury_ACTtot_spinup = sum(CCentury_pools_spinup)
  
  #Proportion of C in each pool, relative to total C
  prop_ACTpools_CCentury = CCentury_pools_spinup/CCentury_ACTtot_spinup
  
  #Initial C in each pool after relaxation
  CCentury_ACTpools_relax = SOC_0*prop_ACTpools_CCentury
  
  
  #####
  #3)## Forward
  #####
  InCent=data.frame(time=ts_fwd_Cent,In=InCent_fwd)
  surfaceInCent=data.frame(time=ts_fwd_Cent,surfaceIn=InCent_fwd*abvgIn)
  soilInCent = data.frame(time=ts_fwd_Cent,soilIn=InCent_fwd*blgIn)

  Century_fwd = CenturyModel(t=ts_fwd_Cent,ks=decomposition_param_Century,C0=as.numeric(CCentury_ACTpools_relax),surfaceIn=surfaceInCent, soilIn = soilInCent, LN=lignin_to_nitrogen, Ls=structural_in_lignin,
                             clay=clayCent,silt=siltCent,xi=xi.frame_Century,
                             solver=deSolve.lsoda.wrapper)
  
  #Get SOC fwd
  CCentury_fwd=getC(Century_fwd) #MgC/ha
  
  if(plot_figures==TRUE){
    #plot the pools
    matplot(ts_fwd_Cent, CCentury_fwd, type="l", lty=1, col=1:7,
            xlab="Time (months)", ylab="C stocks (MgC/ha)",main="Century forward")
    legend("topleft", legend_Cent,
           lty=1, col=1:7, bty="n")
  }
  
  #Outputs of C
  Rt1Century=getReleaseFlux(Century_fwd) #MgC/ha/month
  
  #Cumulative Outputs of C
  Rc1Century=getAccumulatedRelease(Century_fwd)
  
  return(list(Rt1Century,CCentury_fwd))
}
