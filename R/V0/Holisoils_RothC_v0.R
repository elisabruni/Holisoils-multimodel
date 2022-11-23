###################################################
#RothC model
###################################################
#This function launches the SoilR version of RothC
#First, a spin-up run is performed
#Then, the estimated C pool fractions are used to rescale total SOC stocks with measured SOC stocks
#Finally, a forward run is perfomed
#The function provides SOC stock and CO2 fluxes for the chosen simulation lenght

#RothC
#' Jenkinson, D. S., S. P. S. Andrew, J. M. Lynch, M. J. Goss, and
#' P. B. Tinker. 1990. The Turnover of Organic Carbon and Nitrogen in Soil.
#' Philosophical Transactions: Biological Sciences 329:361-368.
#' 


Call_RothC<-function(plot_figures,
                     temperature, precipitation, potential_evapotranspiration,
                     SOC_0,C_input_spinup,C_input_fwd,clay_p,soil_thickness,
                     decomposition_param_RothC,
                     t_spinup,t_fwd,
                     spinupcheck, thresholdspin){
  

  #Convert variables in RothC
  InRothC_spinup = C_input_spinup #MgC/yr
  ts_spinup_RothC = t_spinup      #yr
  
  InRothC_fwd = C_input_fwd       #MgC/yr
  ts_fwd_RothC = t_fwd            #yr
  
  #--
  #Rate modifiers
  fT=fT.RothC(temperature$Temp) #Temperature effects per month
  fW=fW.RothC(P=(precipitation$Precip), E=(potential_evapotranspiration$Potevap), #CHANGE TO HAVE NEUTRAL DATABASE (without Temp and $Precip)
              S.Thick = soil_thickness, pClay = clay_p, 
              pE = 1.0, bare = FALSE)$b #Moisture effects per month
  
  #Rate modifier spinup
  ##############################################################################
  #MODIFY FOLLOWING LINE FOR EFFICIENCY: xi.frame_spinup_RothC=mean(fT*fW)
  ##############################################################################
  xi.frame_spinup_RothC=data.frame(ts_spinup_RothC,rep(mean(fT*fW),length.out=length(ts_spinup_RothC)))
  
  #Rate modifier forward
  xi.frame_RothC=data.frame(ts_fwd_RothC,rep(fT*fW,length.out=length(ts_fwd_RothC)))
  
  #####
  #1)## Initialization
  #####
  
  FallIOM=0.049*SOC_0^(1.139) #IOM using Falloon method
  
  #Run spinup RothC
  RothC_spinup=RothCModel(t=ts_spinup_RothC,ks=decomposition_param_RothC,C0=c(DPM=0, RPM=0, BIO=0, HUM=0, IOM=FallIOM),
                          In=InRothC_spinup,clay=clay_p, xi=xi.frame_spinup_RothC)
  
  #Get SOC spinup
  CRothC_spinup=getC(RothC_spinup)
  
  #Check that steady state is reached
  INIZ=CRothC_spinup[nrow(CRothC_spinup)-spinupcheck,] #initialize SOC pool values
  for(row in 1:nrow(tail(CRothC_spinup,spinupcheck))){ 
    pools_i = tail(CRothC_spinup,spinupcheck)[row,]
    deltai = (INIZ-pools_i)/INIZ
    
    if(all(deltai<thresholdspin)){ #Check that SOC stock variation of each pool is <thresholdspin for all years
      print("spinup ok")}else{ #Otherwise stop and increase the spinup length
        print(paste(c("current delta is",deltai,collapse = " ")))
        stop("spinup length should be increased, current annual SOC variation is ")
      }
    INIZ=pools_i
  }
  
  if(plot_figures==TRUE){
    #plot the pools
    matplot(ts_spinup_RothC, CRothC_spinup, type="l", lty=1, col=1:5,
            xlab="Time (years)", ylab="C stocks (MgC/ha)",main="RothC spinup")
    legend("topleft", c("DPM", "RPM", "BIO", "HUM", "IOM"),
           lty=1, col=1:5, bty="n")
  }
  
  
  #####
  #2)## Relaxation
  #####
  #C after spinup in each pool
  #CRothC_pools_spinup = colMeans(tail(CRothC_spinup,12))
  CRothC_pools_spinup = tail(CRothC_spinup,1)
  #Average total C after spinup
  CRothC_tot_spinup = sum(CRothC_pools_spinup)
  
  #Average total active C after spinup
  CRothC_ACTtot_spinup = sum(CRothC_pools_spinup[1:4])
  
  #Proportion of C in each active pool, relative to total C
  prop_ACTpools_CRothC = CRothC_pools_spinup[1:4]/CRothC_ACTtot_spinup
  
  #Initial C in each pool after relaxation
  CRothC_ACTpools_relax = (SOC_0-FallIOM)*prop_ACTpools_CRothC
  CRothC_pools_relax = c(CRothC_ACTpools_relax,FallIOM)
  
  
  #####
  #3)## Forward
  #####
  
  InRothC=data.frame(time=ts_fwd_RothC,In=InRothC_fwd)
  
  FYM=data.frame(time=ts_fwd_RothC,FYM=rep(0,length(InRothC_fwd)))
  
  RothC_fwd=RothCModel(t=ts_fwd_RothC,ks=decomposition_param_RothC,C0=CRothC_pools_relax,
                       In=InRothC,FYM=FYM, clay=clay_p, xi=xi.frame_RothC)
  
  #Convert variables out
  #Get SOC fwd
  CRothC_fwd=getC(RothC_fwd)
  
  if(plot_figures==TRUE){
    #plot the pools
    matplot(ts_fwd_RothC, CRothC_fwd, type="l", lty=1, col=1:5,
            xlab="Time (years)", ylab="C stocks (MgC/ha)",main="RothC forward")
    legend("topleft", c("DPM", "RPM", "BIO", "HUM", "IOM"),
           lty=1, col=1:5, bty="n")
  }
  
  #Outputs of C
  Rt1RothC=getReleaseFlux(RothC_fwd) 

  #Cumulative Outputs of C
  Rc1RothC=getAccumulatedRelease(RothC_fwd)
  
  return(list(Rt1RothC,CRothC_fwd))
}


                 
