###################################################
#AMG model
###################################################
#This function launches AMG, using the OnePoolModel function in SoilR
#First, the model is initialized considering a default and constant fraction of total C that is stable, which depends on the historical land use
#Then, a forward run is perfomed
#The function provides SOC stock and CO2 fluxes for the chosen simulation lenght



Call_AMG<-function(plot_figures,
                   temperature, precipitation, potential_evapotranspiration,
                   SOC_0,C_input_fwd,clay_p,carbonat_p,
                   decomposition_param_AMG,
                   historical_land_use,
                   t_fwd){
  
  #Convert variables in AMG
  InAMG_fwd = C_input_fwd   #MgC/yr
  ts_fwd_AMG = t_fwd        #yr
  clayAMG = clay_p*10       #% to g/kg
  carbonatAMG = carbonat_p  #g/kg
  
  #CHANGE if possible with real data
  #To calculate aboveground and belowground inputs
  APPT = mean(precipitation[,2])
  Pmax=-40+7.7*APPT # Max aboveground production
  Rmax=100+7.0*APPT # Max belowground production
  abvgIn=Pmax/(Pmax+Rmax) #proportion of aboveground (?)
  blgIn=Rmax/(Pmax+Rmax)  #proportion of belowground (?)
  
  #C input are multiplied by an humification coefficient before entering the active pool
  HumifiedInAMG_fwd <- InAMG_fwd*abvgIn*as.numeric(decomposition_param_AMG["humABOVE"])+
    InAMG_fwd*blgIn*as.numeric(decomposition_param_AMG["humBELOW"])
  
  #--
  #Rate modifiers

  #Create annual dataframe for mean temperature
  temperature_annual <- aggregate(temperature["Temp"], list(format(temperature$Date,"%Y")),FUN=mean, na.rm = TRUE)
  
  #Create annual dataframe for precip
  precipitation_annual <- aggregate(precipitation["Precip"],list(format(precipitation$Date,"%Y")),FUN=sum, na.rm = TRUE)

  #Create annual dataframe for potevap
  potential_evapotranspiration_annual <- aggregate(potential_evapotranspiration["Potevap"],list(format(precipitation$Date,"%Y")),FUN=sum, na.rm = TRUE)
  

  fT=fT.AMG(temperature_annual$Temp) #Temperature effects per year
  fW=fW.AMG(precipitation_annual$Prec,potential_evapotranspiration_annual$Potevap) #Moisture effects per year
  fCl=fCl.AMG1(clayAMG)
  fCaCO3=fCaCO3.AMG1(carbonatAMG)
  
  
  #Rate modifier forward
  xi.frame_AMG=data.frame(ts_fwd_AMG,rep(fT*fW*fCl*fCaCO3,length.out=length(ts_fwd_AMG)))

  
  #####
  #1)## Initialization
  #####
  #Proportion of stable C in total C:
  if(grepl("arable", historical_land_use, ignore.case = TRUE)){ #if historical_land_use = "arable" (not case sensitive)
    stab.propr=0.65    #if historical use of the land is arable
  }else if(grepl("grassland", historical_land_use, ignore.case = TRUE) | grepl("forest", historical_land_use, ignore.case = TRUE)){ #if historical_land_use = "grassland" (not case sensitive)
    stab.propr=0.4     #if historical use of the land is grassland
  }else{
    print("What? Choose between `arable` and `grassland`")
    stop()
  }
  
  #Proportion of stable C
  CAMG_STABpool = SOC_0*stab.propr
  #Proportion of active C at the beginning of the simulation
  CAMG_ACTpool = SOC_0-CAMG_STABpool
  
  #2. no relaxation
  
  #####
  #3)## Forward
  #####


  InAMG=data.frame(time=ts_fwd_AMG,In=HumifiedInAMG_fwd)
  
  AMG_fwd=OnepModel(t=ts_fwd_AMG, k=as.numeric(decomposition_param_AMG["k0"]),
                           C0=CAMG_ACTpool,
                           In=InAMG,xi=xi.frame_AMG)
  
  #Convert variables out
  #Get SOC fwd
  CAMG_fwd_act=getC(AMG_fwd) #Active C stock
  CAMG_fwd = cbind(CAMG_fwd_act,rep(CAMG_STABpool,length(CAMG_fwd_act))) #Active and Stable C pools

  if(plot_figures==TRUE){
    #plot the pools
    matplot(ts_fwd_AMG, CAMG_fwd, type="l", lty=1, col=1:2,
            xlab="Time (years)", ylab="C stocks (MgC/ha)",main="AMG forward")
    legend("topleft", c("Active","Stable"),
           lty=1, col=1:2, bty="n")
  }
  
  #Outputs of C
  Rt1AMG_act=getReleaseFlux(AMG_fwd)
  Rt1AMG=cbind(Rt1AMG_act,rep(0,length(Rt1AMG_act))) #add stable pool fluxes (=0)
  
  #Cumulative Outputs of C
  Rc1AMG_act=getAccumulatedRelease(AMG_fwd)
  Rc1AMG = cbind(Rc1AMG_act,rep(0,length(Rc1AMG_act))) #add stable pool fluxes (=0)
  
  return(list(Rt1AMG,CAMG_fwd))
}
