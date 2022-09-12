#----------------------------------------------------------
#AMG environmental functions affecting decomposition
#----------------------------------------------------------

fT.AMG <- function(Temp){
  #Temp = mean annual temperature
  aT<-25
  cT<-0.120 #1/Kelvin
  Tref <- 15 #Celsius
  bT <- (aT-1)*exp(cT*Tref)
  
  k_temp <- ifelse(Temp>0, aT/(1+bT*exp(-cT*Temp)), 0)
  
  return(k_temp)
}



fW.AMG <- function(Prec,Potevap){
  #Prec = cumulative annual water inputs (precipitation and irrigation water) (mm)
  #pot_evapot = potential evapotranspiration (mm)
  
  aH = 3/100
  bH = 5.247 #1/m
  
  k_precip = 1/(1+ aH*exp(-bH*(Prec-Potevap)/1000))
  
  return(k_precip)
}

fCl.AMG1<-function(Clay){
  #clay (gclay/kgsoil)
  
  aM = 2.720/1000 #g/kg
  k_clay = exp(-aM*Clay)

  return(k_clay)
}

fCl.AMG2<-function(Clay){
  #clay (gclay/kgsoil)
  
  aM =2.519/1000 #AMGv2
  k_clay = exp(-aM*Clay)
  
  return(k_clay)
}

fCaCO3.AMG1<-function(Carbonate){
  #carbonate (CaCO3) content (gCaCo3/kgsoil)
  
  cM = 1.67/1000  #g/kg

  k_carbonate = 1/(1+cM*Carbonate)
  return(k_carbonate)
}

fCaCO3.AMG2<-function(Carbonate){
  #carbonate (CaCO3) content (gCaCo3/kgsoil)
  
  cM = 1.50/1000 #AMGv2
  
  k_carbonate = 1/(1+cM*Carbonate)
  return(k_carbonate)
}

fpH.AMG2<-function(PH){
  apH = 0.112
  bpH = 8.5
  k_pH = exp(-apH*(pH-bpH)^2)
  return(k_pH)
}

fCN.AMG2<-function(CNratio){
  aCN = 0.060
  bCN = 11.
  CN_ratio_func = 0.8*exp(-aCN*(CN_ratio-bCN)^2) + 0.2
  return(k_CN)
}


