## SoilR Yasso20Model #######################################################################

# based on Yasso07 script by Boris Tupek boris.tupek@luke.fi
#Modified with environmental functions from Viskari et al 2022
#https://doi.org/10.5194/gmd-15-1735-2022

#Elisa Bruni ebruni93@gmail.com 2022

#Import SoilR library
library(SoilR)


#Yasso20 as in Viskari et al 2022
Yasso20Modelfi <- function (t, #years
                            #decomposition rates
                            #------------
                            #Viskari et al 2022
                            #------------
                            # #decomposition params
                            ksY = c(kA = 0.51, kW = 5.19, kE = 0.13, kN = 0.1, kH = 0.0015),
                            #transfers and feedbacks
                            # pY = c( p1 = 0.5, p2 =0, p3 = 1,
                            #         p4 = 1, p5 = 0.99, p6 = 0.,
                            #         p7 = 0, p8 = 0., p9 = 0.,
                            #         p10 = 0, p11 = 0.163, p12 = 0.,
                            #         pH = 0.0042),
                            #TESTS on Viskari parametrization
                            #Problem here is with pools A and N
                            #Decreasing for instance p4 (transfer from A to W) and p3 (trnasfer from N to A) solves the issue
                            #transfers and feedbacks
                            #Now p3 and p4 are defined as in Viskari et al., 2022 = 1-pH
                            #But we substract epsilon = 1e-16 to compensate for R calculation approximations
                            pY = c( p1 = 0.5, p2 =0., p3 = 1-0.0042-1e-15,
                                    p4 =1-0.0042-1e-15, p5 = 0.99, p6 = 0.,
                                    p7 = 0., p8 = 0., p9 = 0.,
                                    p10 = 0., p11 = 0.163, p12 = 0.,
                                    pH = 0.0042),
                            #environmental dependence params
                            beta1 = c(beta1AWE=0.158, beta1N=0.17, beta1H=0.067),
                            beta2 = c(beta2AWE=-0.002, beta2N=-0.005, beta2H=0.),
                            gamma = c(g=-1.44, gN=-2.0,gH=-6.9),
                            # Woody litter size dependence parameters
                            delta1 = -2.55,
                            delta2 = 1.24,
                            r =  -0.25,
                            #------------
                            #sample_param
                            #------------
                            # #decomposition params
                            # ksY = c(kA = 4.897147e-01, kW = 4.913873e+00, kE = 2.419735e-01, kN = 9.487642e-02,
                            #         kH = 1.302583e-03),
                            # #transfers and feedbacks
                            # # pY = c(p1 = 4.362893e-01, p2 =2.499740e-01, p3 = 9.151269e-01, p4 = 9.925823e-01,
                            # #        p5 = 8.385374e-02, p6 = 1.147678e-02, p7 = 6.083150e-04, p8 = 4.761282e-04, p9 = 6.603773e-02, p10 = 7.713417e-04,
                            # #        p11 = 1.040174e-01, p12 = 6.488076e-01, pH = 4.596472e-03),
                            # #Problem here is with pool N
                            # # Can be solved by changing p3 in sample_parameters!!
                            # #p3 needs to be set to 0.36 max to respect mass balance
                            # #https://github.com/YASSOmodel/Ryassofortran/tree/master/data
                            # pY = c(p1 = 4.362893e-01, p2 =2.499740e-01, p3 = 0.4, p4 = 9.925823e-01,
                            #        p5 = 8.385374e-02, p6 = 1.147678e-02, p7 = 6.083150e-04, p8 = 4.761282e-04, p9 = 6.603773e-02, p10 = 7.713417e-04,
                            #        p11 = 1.040174e-01, p12 = 6.488076e-01, pH = 4.596472e-03),
                            # #environmental dependence params
                            # beta1 = c(beta1AWE=9.059805e-02, beta1N=4.877247e-02, beta1H=3.518549e-02),
                            # beta2 = c(beta2AWE=-2.144096e-04, beta2N=-7.913602e-05, beta2H=-2.089906e-04),
                            # gamma = c(g=-1.808920e+00, gN=-1.172547e+00,gH=-1.253595e+01),
                            # # # Woody litter size dependence parameters
                            # delta1 = -4.389227e-01,
                            # delta2 = 1.267467e+00,
                            # r =  -2.569142e-01,
                            C0, #initial C , 5 element vector, e.g. c(0,0,0,0,0) 
                            In, # litter C input, data.frame(years, litter), same length as years, if AWEN faractionatoin not provided it has to be in 5 element form for each time step
                            AWEN, #5 element vector, fractionation of plant C litter input to yasso AWEN pools
                            xi = 0, # x1 != 1 will use climate data
                            # xi = 1  will ignore climate data, no climate effect,
                            MTm,# MeanTemperature
                            PR_mm, # Precipitation_mm
                            WS, # woody size, 0 no effect  for nonwoody, 2 finewoody, 20 coarse woody
                            solver = deSolve.lsoda.wrapper, 
                            pass = FALSE){
  
  t_start = min(t)
  t_end = max(t)
  
  # print("=====")
  # print(MTm)
  # print(PR_mm)
  # print("=====")
  
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
  
  print(" ")
  print("Ap")
  print(Ap)
  print("Sum of columns Ap")
  print(colSums(Ap))
  
  # add Woody litter size dependence to structural matrix
  # WS in cm, e.g. 0 for nonwoody, 2 for finewoody, 20 - for coarsewoody litter
  # the effect of wl size as in Eq. 3.1 in model description
  AYS.wsfun <- function(WS){
    ksY.ws = c(ksY[1:4]*(1+delta1*WS+delta2*(WS^2))^(-abs(r)),ksY[5])
    print(ksY.ws)
    A1.ws = abs(diag(ksY.ws))
    #AYS.ws = Ap %*% A1.ws
    #print(AYS.ws)
    return(A1.ws)
  }
  AYS.ws <- AYS.wsfun(WS)
  print(" ")
  print("AYS.ws")
  print(AYS.ws)
  print("Column sums AYS.ws")
  print(colSums(AYS.ws))
  
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
  y20.Efun <- function(MTmX,PR_mmX){
    #MTm = dataframe of Mean Temp monthly, as a list for multiple years
    #PR_mm = Mean Precip annual (mm)
    # print("***")
    # print(MTmX)
    # print("...")
    # print(PR_mmX)
    # print("***")
    
    func_temp_months<-function(beta1_pool,beta2_pool,monthTemp){
      #print(exp(beta1_pool*monthTemp+beta2_pool*monthTemp^2))
      return(exp(beta1_pool*monthTemp+beta2_pool*monthTemp^2))
    }
    
    
    TS_dep= sum(func_temp_months(beta1["beta1AWE"],beta2["beta2AWE"],MTmX))
    TS_depN= sum(func_temp_months(beta1["beta1N"],beta2["beta2N"],MTmX))
    TS_depH= sum(func_temp_months(beta1["beta1H"],beta2["beta2H"],MTmX))
    
    temp_dep_vec=c(TS_dep,TS_dep,TS_dep,TS_depN,TS_depH)

    #Add precipitation dependence
    ENV_dep = temp_dep_vec[1:3]*(1-exp(gamma["g"]*PR_mmX/1000))/12
    ENV_depN = temp_dep_vec[4]*(1-exp(gamma["gN"]*PR_mmX/1000))/12
    ENV_depH = temp_dep_vec[5]*(1-exp(gamma["gH"]*PR_mmX/1000))/12

    ENV_vector =  as.numeric(c(ENV_dep,ENV_depN,ENV_depH))
    # ENV_vector =  rep(1.92406,5)
    # print("ENV_vector")
    # print(ENV_vector)
    return(ENV_vector)
  }
  
  #Add environmental effect
  xiE=function(t){
    # print("=====")
    # print(y20.Efun(MTm[t,], PR_mm[t]))
    # print(MTm)
    # print(MTm[t,])
    # print(PR_mm)
    # print(PR_mm[t])
    y20.Efun(MTm[t,], PR_mm[t])
  }

  
  #if xi = 1 replace climate by no effect
  if(xi == 1){
    #woody size but NO environtal effect 
    #if woody size 0 than no woody size too
    AYS.Ews_t=BoundLinDecompOp(
      function(t){Ap %*%AYS.ws},
      t_start,
      t_end
    )
  } else{
    #woody size and environtal effect
    AYS.Ews_t=BoundLinDecompOp(
      function(t){
        if(length(xiE(t))==0){
          Ap *xiE(t)*AYS.ws
        }else{
          print("^^^^^^^^^^^^^^^^^^^^")
          print("^^^^^^^^^^^^^^^^^^^^")
          print("^^^^^^^^^^^^^^^^^^^^")
          # print(is.double(diag(xiE(t))))
          # print(is.double(AYS.ws))
          # print(is.double(Ap))
          print(diag(xiE(t)))
          print(AYS.ws)
          print(" ")
          print("Matrix xiE ")
          print(diag(xiE(t)))
          print(" ")
          print("Matrix AYS.Ews_t ")
          print(Ap %*%diag(xiE(t))%*%AYS.ws)
          print("Sum of columns AYS.Ews_t ")
          print(colSums(Ap %*%diag(xiE(t))%*%AYS.ws))
          print("Sum of rows AYS.Ews_t ")
          print(rowSums(Ap %*%diag(xiE(t))%*%AYS.ws))
          Ap %*%diag(xiE(t))%*%AYS.ws
        }
        # print("xiE(t)")
        # print(xiE(t))
        # print("xiE(t)*AYS.ws")
        # print(xiE(t)*AYS.ws)
        # print(colSums(xiE(t)*AYS.ws))
        # xiE(t)*AYS.ws
      }
      ,
      t_start,
      t_end
    )
  }
  
  #Yasso20 on Carlos Sierra's general model 
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

TempData=data.frame(seq(from=1,to=10*12,by=1),Temp=15+sin(2*pi*years)+
                      rnorm(n=length(years),mean=0,sd=1))
TempData=data.frame("Date"=seq(from = as.Date("2001-01-01"), to = as.Date("2010-12-31"), by = 'month'),"Temp"=TempData[,"Temp"])

Temp_month_split =
  do.call(rbind,split(TempData$Temp,format(TempData$Date,format="%Y")))

Precipitation <- rep(800,length(years)) # precipitation 800mm

# #Convert daily to monthly temperature
# Temp_month = aggregate(Temp_day, list(format(Temp_day$Date,"%Y-%m")),mean, na.rm = TRUE)
# Precip_year <-
#   aggregate(Precip_day["Precip"], list(format(Precip_day$Date,"%Y")),FUN=sum, na.rm = TRUE)
# PR_in=Precip_year$Precip

MTm_in=Temp_month_split
PR_in=Precipitation
  
yassofix20 <- Yasso20Modelfi(years,
                             C0=rep(0,5), #initial carbon
                             AWEN = c(0.52,0.18,0.08,0.2,0), #to separate litter to yasso AWEN pools, this depends on plant organ and species
                             In=Litter,#litter C input (same length as years)
                             xi = 0, # only xi = 1  will replace climate data no climate effect,
                             MTm=MTm_in,#MeanTemperature
                             PR_mm=PR_in,#Precipitation_mm)
                             WS=2)
getC(yassofix20)
