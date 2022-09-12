## Corrected SoilR Yasso07Model #######################################################################
#install.packages("SoilR", dependencies = T)
library(SoilR)
#Yasso07Model #original version

#the original SoilR version IS MISSINNG FRACTIONATION TO AWEN!!!!!
#as a result C input is not distributed to A W E N pools but all enters pool A

# FIXED SoilR Yasso07model function includes:
# 1) AWEN fractionation
# 2) decomopsition dependency on size of litter
# 3) original environmental functions used for model calibration

# by Boris Tupek boris.tupek@luke.fi
# Dec 2019


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


yasso.matrix.fun <- function(WS,# 0,2,20 cm, when 0 ignored
                             clim, #MT, TA, PR_mm, if clim = 1 ignored
                             wetlands, # if wetlands = "y", apply 35% reduction of decomposition se Kleinen et al. (2021)
                             A.print #if A.print = "y", prints the structural matrix, "n" ignored
){
  # for below to work
  # read the structure functions to the environment from the model!
  #decomposition rates
  ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031, 
          kH = 0.0017) 
  #transfers and feedbacks
  pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99, 
         p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0, 
         p11 = 0.01, p12 = 0.92, pH = 0.0045)
  # climate dependence parameters
  beta1 = 0.096 
  beta2 = -0.0014 
  gamma = -1.21 
  # Woody litter size dependence parameters
  delta1 = -1.7 
  delta2 = 0.86 
  r =  -0.306
  
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
  
  AYS = Ap %*% abs(diag(ksY))
  
  if(wetlands=="y"){
    AYS = Ap %*% abs(diag(0.35*ksY))
    return(AYS)
  }
  
  if(A.print=="y"){
    print("default structural matrix")
    print(AYS)
  }
  
  # add Woody litter size dependence to structural matrix
  # WS in cm, e.g. 0 for nonwoody, 2 for finewoody, 20 - for coarsewoody litter
  # the effect of wl size as in Eq. 3.1 in model description
  
  AYS.wsfun <- function(WS){
    ksY.ws = c(ksY[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),ksY[5])
    A1.ws = abs(diag(ksY.ws))
    AYS.ws = Ap %*% A1.ws
    if(A.print=="y"){
      print("structural matrix modified by woody size")
      print(AYS.ws)
    }
    return(AYS.ws)
  }
  AYS.ws <- AYS.wsfun(WS)
  
  #yasso environmental function
  if(length(clim) ==3){
    
    MT = clim[1]
    TA = clim[3]
    PR_mm = clim[2]
    
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
    xiE <- y07.Efun(MT, TA, PR_mm)
    AYS.wsE <- xiE*AYS.ws
  } else{
    AYS.wsE <- 1*AYS.ws
  }
  if(A.print=="y"){
    print("returns structural matrix modified by woody size and climate")
    print(AYS.wsE)
  }
  return(AYS.wsE)
}

yasso07.soilr.fi.example <- function(){

years=seq(from=1,to=10,by=1)#/365)
Litter=data.frame(year=c(1:10),Litter=rnorm(n=10,mean=10,sd=2))
TempData=data.frame(years,Temp=15+sin(2*pi*years)+
                      rnorm(n=length(years),mean=0,sd=1))
j=length(years) # how many years we simulate 
MeanTemperature <- TempData$Temp  
TemperatureAmplitude <- rep(diff(range(TempData[,2]))/2,j) # temperature amplitude 
Precipitation <- rep(800,j) # precipitation 800mm

MT=MeanTemperature
TA=TemperatureAmplitude
PR_mm=Precipitation
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

Ct=getC(yassofix)
Rt=getReleaseFlux(yassofix) #respiration

par(mfrow=c(2,2))
#plot carbon pools
matplot( Ct, type="l", lty=1, col=1:5, xlab="Years", ylab="Carbon stocks", main ="YM07fix(envir.)")
legend("topleft", c("A","W","E","N", "H"), lty=1, col=1:5, bty="n", n = 2)
#plot respiration
matplot(years, Rt, type="l", ylab="Resiration", lty=1, col=1:2)

# Modify litter input in form of AWEN
LitterAWEN <- c(0.52,0.18,0.08,0.2,0) 
LA =  matrix(LitterAWEN , nrow=j, ncol=5, byrow=TRUE) 
LI.awen = as.data.frame(LA*as.vector(Litter[,2])) #LitterInput
names(LI.awen) = c("A","W","E","N", "H")
LI.y.awen = data.frame(years, LI.awen)
as.matrix(LI.y.awen[,2:6])

yassofix2 <- Yasso07Modelfi(years,
                           C0=rep(0,5), #initial carbon
                           AWEN = 0, # or AWEN e.g. c(0.52,0.18,0.08,0.2,0), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                           In=LI.y.awen ,#litter C input (same length as years)
                           xi = 0, # only xi = 1  will replace climate data no climate effect,
                           MT=MT,#MeanTemperature
                           TA=TA, #TemperatureAmplitude
                           PR_mm=PR_mm,#Precipitation_mm)
                           WS=2) 

Ct2=getC(yassofix2)
Rt2=getReleaseFlux(yassofix2) #respiration

#par(mfrow=c(1,2))
#plot carbon pools
matplot( Ct2, type="l", lty=1, col=1:5, xlab="Years", ylab="Carbon stocks", main ="YM07fix(envir.)")
legend("topleft", c("A","W","E","N", "H"), lty=1, col=1:5, bty="n", n = 2)
#plot respiration
matplot(years, Rt2, type="l", ylab="Resiration", lty=1, col=1:2)

Ct
Ct2
identical(Ct,Ct2)


### model diagnostics ##################

yasso.matrix.fun(WS = 2, clim = c(5, 500, 7), A.print = "y")
yasso.matrix.fun(WS = 0, clim = 1, A.print = "n")

AYS.wsfun <- function(WS){
  ksY.ws = c(ksY[1:4]*(1+delta1*WS+delta2*(WS^2))^(r),ksY[5])
  #print(ksY.ws)
  A1.ws = abs(diag(ksY.ws))
  AYS.ws = Ap %*% A1.ws
  #print(AYS.ws)
  return(AYS.ws)
}

#decomposition rates
ksY = c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031, 
        kH = 0.0017) 
#transfers and feedbacks
pY = c(p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99, 
       p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0, 
       p11 = 0.01, p12 = 0.92, pH = 0.0045)
# climate dependence parameters
beta1 = 0.096 
beta2 = -0.0014 
gamma = -1.21 
# Woody litter size dependence parameters
delta1 = -1.7 
delta2 = 0.86 
r =  -0.306
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

#calculate steady states for different litter sizes
#note: quantity litter input and  quality awens should differ for different types
u.Li <- 5 * c(0.52,0.18,0.08,0.2,0) # 4 ton c /ha litter times AWEN fractions 
plot(c(1,3),c(0,100), col ="white", xlab="litter size 0,2,20 cm",ylab ="C steady state" )
#compute steady state
for(i in 1:3){
  AYS.ws <- yasso.matrix.fun(WS = c(0,2,20)[i], clim = 1, A.print = "n")
  xss=-1*solve(AYS.ws)%*%u.Li #inverse of matrix solve(B)
  points(rep(c(1.5,2,2.5)[i],5), xss, col=1:5, pch = 16)
}
legend("topleft", c("A","W","E","N", "H"), pch =16, col=1:5, bty="n", n = 2)

#functionality ages and transit times
tau=seq(0,500)
for(i in 1:3){
  #i= 1
  AYS.ws <- AYS.wsfun(c(0,2,20)[i])
  SA=systemAge(A=AYS.ws, u=u.Li, a=tau)
  
  plot(tau, SA$systemAgeDensity, type="l")
  abline(v=SA$meanSystemAge, lwd=2) #mean
  #abline(v=SA$quantilesSystemAge[2], col=2) #medium
  matplot(tau, SA$poolAgeDensity, col = 1:5, type = "l", add = T)
  abline(v=SA$meanPoolAge, col = 1:5, lwd=2) #mean
}
legend("topright", c("A","W","E","N", "H"), pch =16, col=1:5, bty="n", n = 2)

TT=transitTime(A=AYS.ws, u=u.Li, a=tau)
names(TT)
plot(tau, TT$transitTimeDensity, type="l")
abline(v=TT$meanTransitTime)


}
cat("yasso07.soilr.fi.example() to see some examples: \n
      model runs with different litter input form, \n
      steady state estimation from the matrix and litter input \n
      model diagnostics ")



cat("\n 
\n description: Yasso07Modelfi(years, #if using monthly or daily time step modify temperature function 
\n                       C0=rep(0,5), #initial carbon
\n                       AWEN = c(0.52,0.18,0.08,0.2,0), #fractionation of litter to yasso A W E N H pools, 
\n                            #this depends on plant organ and species
\n                            #if 0 then ignored, and ignored if littter is fractionated data.frame(years, A W E N H) 
\n                       In=Litter,#if litter C input not fractionated data.frame(years, litter), AWEN must be provided
\n                       xi = 0, #0 will use climate data # xi = 1  will ignore it, no climate effect,
\n                       MT=MT,#MeanTemperature mean annual in Celsius
\n                       TA=TA, #TemperatureAmplitude range(mean monthly tmep)
\n                       PR_mm=PR_mm,#Precipitation_mm (ammual sum))
\n                       WS=2)")

cat("\n 
 \n yasso.matrix.fun returns model matrix depending on woody size and climate 
 \n  for steady state estimation, diagnostics of age and transition time
  \n #WS = 0, clim = 1, A.print = n
  \n returns default matrix ignores WS and climate, does not print
  \n #WS = 2, clim = c(5, 500, 7), A.print = y
  \n returns WS and clim ~ matrix and print matrix in the console")
    
#yasso07.soilr.fi.example()
