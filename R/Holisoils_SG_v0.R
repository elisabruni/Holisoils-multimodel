#######################################
#######################################
# Soil GHG model (SG model by Shoji Hashimoto shojih@ffpri.affrc.go.jp )
# Coded in Jan 2022
#Reference 1 (model development): Hashimoto, S. et al. (2011) Simple models for soil CO2, CH4, and N2O fluxes calibrated using a Bayesian approach and multi-site data. Ecological Modelling, 222: 1283-1292. http://dx.doi.org/10.1016/j.ecolmodel.2011.01.013
#Reference 2 (regional application): Hashimoto, S. et al. (2011) Increasing trends of soil greenhouse gas fluxes in Japanese forests from 1980 to 2009. Scientific Reports, 1, 116. https://doi.org/10.1038/srep00116
#Reference 3 (global application): Hashimoto S. (2012) A new estimation of global soil greenhouse gas fluxes using a simple data-oriented model. PLoS ONE, 7, e41962. https://doi.org/10.1371/journal.pone.0041962

#######################################
#Unit of the outputs
#CO2 efflux: ugC m-2 s-1
#CH4 uptake flux: ugC m-2 h-1
#N2O efflux: ugN m-2 h-1

#######################################
#Details of the inputs
#CN Ratio (unitless): 0-5 cm depth
#Bulk density (Mg m-3): 0-5 cm depth
#WFPS (water filled pore space, unitless) @5cm depth
#Soil temperature, (degreeC) @5cm depth

#######################################
# IF you have no WFPS, but bulk density and volumetric soil water content,
#BD2 (bulk dencity, Mg m-3)@5cm depth
#VSWC (volumetiric soil water content, %) @5cm depth
#Then, please calculate WFPS using the function, cal_WFPS()

#######################################
#############  How to run #############
#Minimum requirement
#CN Ratio, Bulk density, WFPS, Soil temperature
#If available, atmospheric CH4 concentration, ppb (e.g. 1900)

#######################################
#######################################
#######################################
# Functions
SGmodel<-function(CN_Ratio, Bulk_Density, WFPS, vswc, Soil_Temperature, CH4_Conc){
  #If WFPS not available, set it to a negative number -> it will be calculated using vswc and bulk density
  #If WFPS is available, vswc won't be used
  
  #Calculate WFPS with vswc and bulk density, if provided as a negative number:
  if(WFPS<0){
    WFPS<-as.numeric((vswc/100)/(1 - (Bulk_Density/2.65)))
  }
  
    m<-vector()
    n<-vector()
    a<-vector()
    b<-vector()
    c<-vector()
    d<-vector()
    p<-vector()

    f<-vector()
    g<-vector()
    h<-vector()
    flux<-c(0,0,0)
    names(flux)<-c("CO2_flux", "CH4_flux", "N2O_flux")

    ##############################
    #CO2 efflux
    i<-1
    m[i]<- 4.95
    n[i]<- 0.033
    a[i]<--0.10
    b[i]<- 0.50
    c[i]<- 2.69
    d[i]<- 1.52
    p[i]<- 0.067

    #CH4 uptake
    i<-2
    m[i]<- 223.2
    n[i]<- 1.76
    a[i]<--1.99
    b[i]<- 0.073
    c[i]<- 1.12
    d[i]<- 4.31
    p[i]<- 0.009

    #N2O efflux
    i<-3
    m[i]<- 48.58
    n[i]<- 0.254
    a[i]<--0.41
    b[i]<- 1.00
    c[i]<- 1.81
    d[i]<- 9.03
    p[i]<- 0.12

    ##############################
    #Atmospheric CH4 concentration
    #Ryori station, Japan average 2002-2004
    CH4_Conc_base<-1858

    if(CH4_Conc>0)
    {
        CH4_Conc_Ratio<-CH4_Conc/CH4_Conc_base
    }
    else
    {
        #Ryori station, Japan, average 2016-2020
        CH4_Conc<-1946
        CH4_Conc_Ratio<-CH4_Conc/CH4_Conc_base
    }

    ##############################
    #Calculate flux
    for (i in 1:3)
    {
        ##########
        #Soil property
        if(i==1)
        {
            f[i]<-m[i]*exp(n[i]*CN_Ratio)
        }
        if(i==2)
        {
            f[i]<-m[i]*exp(-n[i]*Bulk_Density)
        }
        if(i==3)
        {
            f[i]<-m[i]*exp(-n[i]*CN_Ratio)
        }

        ##########
        #Water
        if(WFPS<a[i])
        {
            g[i]<-0.0
        }
        else if (WFPS>c[i])
        {
            g[i]<-0.0
        }
        else
        {
            g[i]<- ((WFPS - a[i]) / (b[i] - a[i]))^ d[i]
            g[i]<-g[i]* ((WFPS - c[i]) / (b[i] - c[i]))^((d[i] * (c[i] - b[i])) / (b[i] - a[i]))
            if (g[i]<=0.0)
            {
                g[i]<-0.0
            }
        }
        ##########
        #Temperature
        h[i]=exp(p[i]*Soil_Temperature)
    }

    #Flux
    for (i in 1:3)
    {
        flux[i]<-f[i]*g[i]*h[i]

        if(i==2)
        {
            flux[i]<-flux[i]*CH4_Conc_Ratio
        }
    }
    return(flux)
}

Call_SG<-function(CN_Ratio, Bulk_Density, WFPS, vswc, Soil_Temperature, CH4_Conc){
  #Runs the SG model for multiple years
  #all variables @5cm depth
  #Soil_Temperature dataframe, first column are the dates, second column are the values of soil temperature (C)
  #CN_Ratio a scalar (unitless)
  #Bulk_Density a scalar (Mg m-3)
  #WFPS a scalar (water filled pore space, unitless)
  #vswc a scalar, (volumetric soil water content, %)
  #CH4_Conc a scalar, Atmospheric CH4 concentration base (ppb)
  temp_vec <- Soil_Temperature[,"Temp"]
  sg_vec <- data.frame(matrix(ncol = 3, nrow = length(temp_vec)))
  colnames(sg_vec)<-c("CO2_flux","CH4_flux","N2O_flux")
  j<-1
  for(i in temp_vec){
    temp_in<-i
    sg_out<-SGmodel(CN_Ratio=CN_Ratio, Bulk_Density=Bulk_Density, WFPS=WFPS, vswc=vswc, Soil_Temperature=temp_in, CH4_Conc=CH4_Conc)
    #print(sg_out)
    sg_vec[j,]<-sg_out
    j<-j+1
  }
  return(sg_vec)
}
#######################################


#######################################
#######################################
############  Model runs   ############
#######################################
#######################################
# Test runs

# With atmospheric CH4
SGmodel(15.52, 0.38, 0.32, 25.4, 11.3, 1858)
SGmodel(CN_Ratio=15.52, Bulk_Density=0.38, WFPS=0.32, vswc=25.4, Soil_Temperature=11.3, CH4_Conc=1858)

temp_vec<-c(10,11,12,11,10)
sg_vec <- data.frame(matrix(ncol = 3, nrow = length(temp_vec)))
colnames(sg_vec)<-c("CO2_flux","CH4_flux","N2O_flux")
j<-1
for(i in temp_vec){
  temp_in<-i
  sg_out<-SGmodel(CN_Ratio=15.52, Bulk_Density=0.38, WFPS=-1, vswc=25.4, Soil_Temperature=temp_in, CH4_Conc=1858)
  #print(sg_out)
  sg_vec[j,]<-sg_out
  j<-j+1
}

# Without atmospheric CH4
SGmodel(15.52, 0.38, 0.32, 25.4, 11.3, -999.0)
SGmodel(CN_Ratio=15.52, Bulk_Density=0.38, WFPS=0.32,vswc=25.4, Soil_Temperature=11.3, CH4_Conc=-999.0)

# If you have only vswc and bd data, then calculate WFS and use it.
#wfps<-cal_WFPS(vswc=25.4, bd=0.55)
SGmodel(15.52, 0.38, -999.0,25.4, 11.3, 1858)

