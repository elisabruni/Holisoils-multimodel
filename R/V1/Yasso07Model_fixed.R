## Corrected SoilR Yasso07Model #######################################################################
#Yasso07Model #original version

#the original SoilR version IS MISSINNG FRACTIONATION TO AWEN!!!!!
#as a result C input is not distributed to A W E N pools but all enters pool A

# FIXED SoilR Yasso07model function includes:
# 1) AWEN fractionation
# 2) decomopsition dependency on size of litter
# 3) original environmental functions used for model calibration

# by Boris Tupek boris.tupek@luke.fi
# Dec 2019

#Import SoilR library
library(SoilR)

#Yasso07 model as in Tuomi et al. 2011 ecological applications
Yasso07Modelfi <- function (t, #years
                            #decomposition rates
                            ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031, 
                                    kH = 0.0017), 
                            #transfers and feedbacks
                            pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99, 
                                   p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0, 
                                   p11 = 0.01, p12 = 0.92, pH = 0.0045),
                            # climate dependence parameters
                            beta1 = 0.096, 
                            beta2 = -0.0014, 
                            gamma = -1.21, 
                            # Woody litter size dependence parameters
                            delta1 = -1.7, 
                            delta2 = 0.86, 
                            r =  -0.306,
                            C0, #initial C , 5 element vector, e.g. c(0,0,0,0,0) 
                            In, # litter C input, data.frame(years, litter), same length as years, if AWEN faractionatoin not provided it has to be in 5 element form for each time step
                            AWEN, #5 element vector, fractionation of plant C litter input to yasso AWEN pools
                            xi = 0, # x1 != 1 will use climate data
                                    # xi = 1  will ignore climate data, no climate effect,
                            MT,# MeanTemperature
                            TA, # TemperatureAmplitude = (mothly temp. range)/2
                            PR_mm, # Precipitation_mm
                            WS, # woody size, 0 no effect  for nonwoody, 2 finewoody, 20 coarse woody
                            solver = deSolve.lsoda.wrapper, 
                            pass = FALSE){
  # AWEN fractionation # see https://en.ilmatieteenlaitos.fi/yasso-download-and-support, 
  # Yasso07 user-interface manual (pdf 450 kB) p.13-14
  t_start = min(t)
  t_end = max(t)
  
  #structural matrix
  Ap = diag(-1, 5, 5)
  Ap[1, 2] = pY[1]
  Ap[1, 3] = pY[2]
  Ap[1, 4] = pY[3]
  Ap[2, 1] = pY[4]
  Ap[2, 3] = pY[5]
  Ap[2, 4] = pY[6]
  Ap[3, 1] = pY[7]
  Ap[3, 2] = pY[8]
  Ap[3, 4] = pY[9]
  Ap[4, 1] = pY[10]
  Ap[4, 2] = pY[11]
  Ap[4, 3] = pY[12]
  Ap[5, 1:4] = pY[13]
  
  # add Woody litter size dependence to structural matrix
  # WS in cm, e.g. 0 for nonwoody, 2 for finewoody, 20 - for coarsewoody litter
  # the effect of wl size as in Eq. 3.1 in model description
  AYS.wsfun <- function(WS){
    ksY.ws = c(ksY[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),ksY[5])
    #print(ksY.ws)
    A1.ws = abs(diag(ksY.ws))
    AYS.ws = Ap %*% A1.ws
    #print(AYS.ws)
    return(AYS.ws)
  }
  AYS.ws <- AYS.wsfun(WS)
  
  #LitterInput
  if (length(In[1,]) == 6){
    LI = as.matrix(In[,2:6]) #first column for years, 2:6 for AWEN
  } else if (length(AWEN) != 5){
    stop("the AWEN litter fractionation 5 element vector must be provided if C litter input is not fractionated to AWEN")
    } else {
    #fractionate liter input to yasso07 A W E N pools
    LA =  matrix(AWEN , nrow=length(t),
                 ncol=5, byrow=TRUE) 
    LI = LA*as.vector(In[,2])  #first column for years
  }
  
  inputFluxes=function(t){
    matrix(nrow = 5, ncol = 1, LI[t,] )
  }
  inputFluxes_tm=BoundInFluxes(
    inputFluxes, 
    t_start,
    t_end)
  
  
  #environmental effect as function for matrix multiplication
  #this can be replaced by any soil TEMP-SWC function
  y07.Efun <- function(MT,TA,PR_mm){
    #MT = Temp
    #TA =TempAmplit
    PR = PR_mm/1000  # conversion from mm to meters
    
    #seasonal approximation of temperature (if annual mean)
    T1=MT+4*TA/pi*(1/sqrt(2)-1)          # Eq. 2.4 in model description
    T2=MT-4*TA/(sqrt(2)*pi)              # Eq. 2.5 in model description
    T3=MT+4*TA/pi*(1-1/sqrt(2))          # Eq. 2.6 in model description
    T4=MT+4*TA/(sqrt(2)*pi)
    TS =  exp(beta1*cbind(T1,T2,T3,T4)+beta2*(cbind(T1,T2,T3,T4)^2))*(1-exp(gamma*PR))
    TS
    apply(TS,1,mean)
  }
  
 #Add environmental effect
  xiE=function(t){
    y07.Efun(MT[t], TA[t], PR_mm[t])
  }
  
  #if xi = 1 replace climate by no effect
  if(xi == 1){
    #woody size but NO environtal effect 
    #if woody size 0 than no woody size too
    AYS.Ews_t=BoundLinDecompOp(
      function(t){AYS.ws},
      t_start,
      t_end
    )
  } else{
    #woody size and environtal effect
    AYS.Ews_t=BoundLinDecompOp(
      function(t){xiE(t)*AYS.ws},
      t_start,
      t_end
    )
  }
  
  #Yasso07 on Carlos Sierra's general model 
  Mod=GeneralModel(t, A=AYS.Ews_t, ivList= C0, inputFluxes =inputFluxes_tm,
                   solver, pass)
  return(Mod)
}

###########
# TEST
###########

# #-------
years=seq(from=1,to=10,by=1)#/365)

#Litter=data.frame(year=c(1:10),Litter=rnorm(n=10,mean=10,sd=2))
Litter=data.frame(year=c(1:10),Litter=rep(8,10))
TempData=data.frame(years,Temp=15+sin(2*pi*years)+
                      rnorm(n=length(years),mean=0,sd=1))
j=length(years) # how many years we simulate 
MeanTemperature <- TempData$Temp  
TemperatureAmplitude <- rep(diff(range(TempData[,2]))/2,j) # temperature amplitude 
Precipitation <- rep(800,j) # precipitation 800mm

MT=MeanTemperature
TA=TemperatureAmplitude
PR_mm=Precipitation
#PR_mm=Precip_year$Precip
#note conversion from mm to meters in the model's environmental function


# EXAMPLE of fixed model ##
# Modified yasso07 C. Sierra general model WITH environmental effect 
yassofix <- Yasso07Modelfi(years,
                           C0=rep(0,5), #initial carbon
                           AWEN = c(0.52,0.18,0.08,0.2,0), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                           In=Litter,#litter C input (same length as years)
                           xi = 0, # only xi = 1  will replace climate data no climate effect,
                           MT=MT,#MeanTemperature
                           TA=TA, #TemperatureAmplitude
                           PR_mm=PR_mm,#Precipitation_mm)
                           WS=2) 
getC(yassofix)
