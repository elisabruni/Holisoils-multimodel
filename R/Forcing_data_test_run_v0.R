#Forcing
#Define input directory
loc_forc = "/PATH/TO/DATA/" #Choose the path where you stored the climate data

#---------------------------
#USER DEFINED PARAMETERS
#Decide wether to plot default figures [TRUE] or not [FALSE]
plot_figures=TRUE

#Define time step forward run [in years]
computation_time_step_fwd = 1/12

#Soil data from LUCAS (GPS: 42.83861 ; -2.39256 - nearby area under revegetated soil)
clay_site = 38 #[%]
silt_site = 50 #[%]
sand_site = 12 #[%]
ph_site = 7.35
CaCO3_site = 403 #[g/kg]
SOC_initial = 57.2129743120269 #[MgC/ha]
soil_OM_thick = 25 #[cm]

#C input data
Cinput_spinup =  3.338 # [MgC/ha/yr] From Mahmoudi et al 2021 under antropogenic disturbance
Cinput_fwd = 3.338 # [MgC/ha/yr] From Mahmoudi et al 2021 under antropogenic disturbance

#Cinput_spinup = 5.83 # [MgC/ha/yr] From Mahmoudi et al 2021 without antropogenic disturbance
#Cinput_fwd = 5.83 # [MgC/ha/yr] From Mahmoudi et al 2021 without antropogenic disturbance

#Default values
lignin_to_nitrogen_ratio=0.5 #[unitless] lignin to nitrogen ratio of the C input
structural_in_lignin_ratio=0.1 #[unitless] structural to lignin ratio of the C input
woodylittersize_scalar=2 #[cm] set to 0 for nonwoody, 2 finewoody, 20 coarse woody
historical_LU = 'grassland' #historical land use of the site. Set either to 'arable', 'grassland' or 'forest'

#All these variables should be measured @5cm depth:
CN_site=10. #[unitless] carbon to nitrogen ratio of the soil at 5cm depth
BD_site=1.3 #[Mg/m3] bulk density of the soil at 5cm depth
water_filled_pore_space=0.32 #[unitless] ratio of the pore space filled with water at 5cm depth. If value not available, set to <0
volumetric_soil_water_cont=25 #[%] volumetric soil water content at 5cm depth. If water_filled_pore_space is available, this will be ignored
CH4_conc_site=1858 #ppb atmospheric CH4 concentration

#Define at what date the simulations start
start_date_simulations_site=as.Date("2021-01-01") #Starting date of the simulations (YYYY-MM-DD)

#Decomposition rates models 
#RothC
ksRothC = c(k.DPM = 10, k.RPM = 0.3, k.BIO = 0.66, k.HUM = 0.02, 
            k.IOM = 0) #(1/yr)
#ICBM
param_ICBM = c(k1=0.8,k2=0.00605,h=0.13) #(1/yr)

#Century
#ksCent = c(STR.surface = 0.076, MET.surface = 0.28, STR.belowground = 0.094, MET.belowground = 0.35, ACT = 0.14, SLW = 0.0038, PAS = 0.00013) #(1/week)
ksCent = c(STR.surface = 1/(12*0.245), MET.surface = 1/(12*0.066), STR.belowground = 1/(12*0.245), MET.belowground = 1/(12*0.066), ACT = 1/(12*0.149), SLW = 1/(12*5.480), PAS = 1/(12*241)) #1/month

#Yasso07
paramYasso07=c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031,kH = 0.0017, #[1/yr]
               beta1=0.096, beta2=-0.0014, gamma=-1.21,
               delta1=-1.7,delta2=0.86,r=-0.306,
               p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99,
               p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0,
               p11 = 0.01, p12 = 0.92, pH = 0.0045)

#AMG
#ksAMG = 0.165 #[1/yr] #Default AMG1 is
ksAMG = c(k0=0.165,humABOVE=0.5,humBELOW=0.4)

###########################
#Read forcing climate data
###########################
#----------------
#Temperature data
#----------------
#Read daily temperature data
temp_df <- read.delim(paste0(loc_forc,"942217_TAS_2001_2005.txt"),sep = ":")[,2]
temp_df <- as.numeric(as.character((temp_df[7:length(temp_df)])))-273.15 #K to C
#Repeat because not enough data
temp_df <- c(temp_df,temp_df) #Repeat 2001-2005 twice
#Create a dataframe for daily temperature
Temp_day <- data.frame("Date"=seq(from = as.Date("2001-01-01"), to = as.Date("2010-12-31"), by = 'day'),"Temp"=temp_df)
#Convert daily temperature to monthly average temperature
Temp_month <- aggregate(Temp_day, list(format(Temp_day$Date,"%Y-%m")),mean, na.rm = TRUE)
#Add dates to dataframe
Temp_month$Date<-as.Date(paste(Temp_month$Group.1,"-01",sep=""))

#----------------
#Precipitation data
#----------------
#Read daily rain data
rain_df <- read.delim(paste0(loc_forc,"942217_PR_2001_2005.txt"),sep = ":")[,2]
rain_df <- as.numeric(as.character(rain_df[7:length(rain_df)]))# "kg m-2 s-1"
#Read daily snow data
snow_df <- read.delim(paste0(loc_forc,"942217_PRSN_2001_2005.txt"),sep = ":")[,2]
snow_df <- as.numeric(as.character(snow_df[7:length(snow_df)])) #"kg m-2 s-1"
#Calculate total precipitation and convert to mm/day
precip_df <- (rain_df+snow_df)*60*60*24 #"kg m-2 s-1" to mm/day
#Repeat because not enough data
precip_df <- c(precip_df,precip_df) #Repeat 2001-2005 twice
#Create a dataframe for daily precipitation
Precip_day <- data.frame("Date"=seq(from = as.Date("2001-01-01"), to = as.Date("2010-12-31"), by = 'day'),"Precip"=precip_df)
#Convert daily precipitation to monthly precipitation
Precip_month <- aggregate(Precip_day["Precip"], list(format(Precip_day$Date,"%Y-%m")),FUN=sum, na.rm = TRUE)
#Add dates to dataframe
Precip_month$Date<-as.Date(paste(Precip_month$Group.1,"-01",sep=""))

#This is to simulate a summer drought
#Create summer drought
Precip_month_drought <- Precip_month
Precip_month_drought$Precip <- with(Precip_month_drought, 
                                    ifelse(sub(".*-", "",Precip_month_drought["Group.1"][[1]])=="05"
                                           |sub(".*-", "",Precip_month_drought["Group.1"][[1]])=="06"
                                           |sub(".*-", "",Precip_month_drought["Group.1"][[1]])=="07", 0, Precip))

#--------------------------------
#Potential evapotranspiration
#--------------------------------
#Read monthly potevap data
potevap_df <- read.delim(paste0(loc_forc,"942217_POT_2001_2005.txt"),sep = ":")[,2]
potevap_df <- as.numeric(as.character(potevap_df[7:length(potevap_df)]))
#Convert to mm/month
potevap_df <-potevap_df*(60*60*24)*rep(c(31,28,31,30,31,30,31,31,30,31,30,31),length(potevap_df)/12)# "kg m-2 s-1" to mm/month
#Repeat because not enough data
potevap_df <- c(potevap_df,potevap_df) #Repeat 2001-2005 twice
#Create dataframe and add dates
Potevap_month <- data.frame("Date"=seq(from = as.Date("2001-01-01"), to = as.Date("2010-12-31"), by = 'month'),"Potevap"=potevap_df)
#Potevap_month$Date<-as.Date(paste(Potevap_month$Group.1,"-01",sep=""))
