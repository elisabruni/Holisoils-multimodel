#Forcing
#Define input directory
inputdir <- "/Users/ebruni/Desktop/HOLISOILS/"
loc_exp = paste0(inputdir,'DATI/Respiration_Bertrand.xlsx')
loc_forc = paste0(inputdir,"DATI/SPAIN/")


#Soil data from LUCAS (GPS: 42.83861 ; -2.39256 - nearby area under revegetated soil)
clay_site = 38 #[%]
silt_site = 50 #[%]
sand_site = 12 #[%]
ph_site = 7.35
CaCO3_site = 403 #[g/kg]
SOC_initial = 57.2129743120269 #tC/ha
soil_OM_thick = 25 #[cm]

#C input data
Cinput_spinup =  3.338 # tC/ha/yr From Mahmoudi et al 2021 under antropogenic disturbance
Cinput_fwd = 3.338 # tC/ha/yr From Mahmoudi et al 2021 under antropogenic disturbance

#Cinput_spinup = 5.83 # tC/ha/yr From Mahmoudi et al 2021 under without antropogenic disturbance
#Cinput_fwd = 5.83 # tC/ha/yr From Mahmoudi et al 2021 under without antropogenic disturbance

#Default values
lignin_to_nitrogen_ratio=0.5
structural_in_lignin_ratio=0.1
woodylittersize_scalar=2 #[0 for nonwoody, 2 finewoody, 20 coarse woody]
historical_LU = 'grassland'

#NEED TO BE @5cm depth:
CN_site=10. #ratio
BD_site=1.3 #Mg/m3
water_filled_pore_space=0.32 #ratio
volumetric_soil_water_cont=25 #%
CH4_conc_site=1858 #ppb


#Decomposition rates models 
#CHANGE FOR PARAMETERS THAT YOU READ FROM DOC
ksRothC = c(k.DPM = 10, k.RPM = 0.3, k.BIO = 0.66, k.HUM = 0.02, 
            k.IOM = 0) #(1/yr)
param_ICBM = c(k1=0.8,k2=0.00605,h=0.13) #(1/yr)
#ksCent = c(STR.surface = 0.076, MET.surface = 0.28, STR.belowground = 0.094, MET.belowground = 0.35, ACT = 0.14, SLW = 0.0038, PAS = 0.00013) #(1/week)
ksCent = c(STR.surface = 1/(12*0.245), MET.surface = 1/(12*0.066), STR.belowground = 1/(12*0.245), MET.belowground = 1/(12*0.066), ACT = 1/(12*0.149), SLW = 1/(12*5.480), PAS = 1/(12*241)) #1/month

#par_rantakari parametrization:
#Rantakari, M., Lehtonen, A., Linkosalo, T., Tuomi, M., Tamminen, P., Heikkinen, J., Liski, J., Mäkipää, R., Ilvesniemi, H. & Sievänen, R. 2012. The Yasso07 soil carbon model - Testing against repeated soil carbon inventory. 2012  Forest Ecology and Management. 286, p. 137-147. doi:10.1016/j.foreco.2012.08.041
#paramYasso07 = c(k.A=0.5172509,k.W=3.551512,k.E=0.3458914,k.N=0.2660175,k.H=2.4180325e-4, #[1/yr]
#                beta1=0.089501545,beta2=-0.0022709155,gamma=-2.935411,
#               delta1=-0.5391662,delta2=1.18574,r=-0.2632936,
#              p1=0.044852223,p2=0.0029265443,p3=0.9779027,p4=0.6373951,p5=0.3124745,
#             p6=0.018712098,p7=0.022490378,p8=0.011738963,p9=9.9046889e-4,p10=0.3361765,
#            p11=0.041966144,p12=0.089885026,p13=0.0015341907)

#par_gui parametrization:
#Tuomi, M., Rasinmäki, J., Repo, A., Vanhala, P. & Liski. J. 2011.Soil carbon model Yasso07 graphical user interface.Environmental Modeling and Software 26 (11): 1358-1362. doi:10.1016/j.envsoft.2011.05.009
#paramYasso07 = c(k.A=0.7035942673683167,k.W=5.681055545806885,k.E=0.2613542377948761,k.N=0.02810959704220295,k.H=0.0014966174494475126, #[1/yr]
#             beta1=0.09873183816671371,beta2=-0.001571640488691628,gamma=-1.2716917991638184,
#            delta1=-1.7084113359451294,delta2=0.8585553765296936,r=-0.3068014085292816,
#           p1=0.4888527989387512,p2=0.019057683646678925,p3=0.9696374535560608,p4=0.9872559905052185,p5=0.0028432635590434074,
#          p6=0.0033964612521231174,p7=1.39937037602067e-5,p8=1.7966924133361317e-5,p9=0.01218125969171524,p10=0.0027778467629104853,
#         p11=0.012695553712546825,p12=0.9713827967643738,p13=0.0042703705839812756)


#Tuomi et al., 2009
#paramYasso07 = c(k.A=0.66,k.W=4.3,k.E=0.35,k.N=0.22,k.H=3.3, #[1/yr]
#                beta1=7.6,beta2=-8.9,gamma=-1.27,
#               delta1=-0.5391662,delta2=1.18574,r=-0.2632936,
#              p1=0.32,p2=0.01,p3=0.93,p4=0.34,p5=0.0035,
#             p6=0.0035,p7=0.0015,p8=0.003,p9=0.01,p10=0.0015,
#            p11=0.03,p12=0.92,p13=0.04)

#SoiLR
#paramYasso07 = c(k.A=0.66,k.W=4.3,k.E=0.35,k.N=0.22,k.H=0.0033, #[1/yr]
#                beta1=0.089501545,beta2=-0.0022709155,gamma=-2.935411,
#               delta1=-0.5391662,delta2=1.18574,r=-0.2632936,
#              p1=0.32,p2=0.01,p3=0.93,p4=0.34,p5=0,
#             p6=0,p7=0.035,p8=0.005,p9=0.01,p10=0.0005,
#            p11=0.03,p12=0.92,p13=0.04)

#Test
#paramYasso07 = c(k.A=0.5172509,k.W=3.551512,k.E=0.3458914,k.N=0.2660175,k.H=2.4180325e-4, #[1/yr]
#               beta1=0.089501545,beta2=-0.0022709155,gamma=-2.935411,
#             delta1=-0.5391662,delta2=1.18574,r=-0.2632936,
#           p1=0.044852223,p2=0.0029265443,p3=0.9779027,p4=0.6373951,p5=0.3124745,
#         p6=0.018712098,p7=0.022490378,p8=0.011738963,p9=9.9046889e-4,p10=0.3361765,
#       p11=0.041966144,p12=0.089885026,p13=0)

#IF SOILR VERSION (BORIS PARAM)
#paramYasso07=c(k.A = 0.73, k.W = 5.8, k.E = 0.29, k.N = 0.031,k.H = 0.0017, #[1/yr]
#              beta1=0.096, beta2=-0.0014, gamma=-1.21,
#             delta1=-1.7,delta2=0.86,r=-0.306,
#            p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99,
#           p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0,
#          p11 = 0.01, p12 = 0.92, p13 = 0.0045)

#IF BORIS VERSION
paramYasso07=c(kA = 0.73, kW = 5.8, kE = 0.29, kN = 0.031,kH = 0.0017, #[1/yr]
               beta1=0.096, beta2=-0.0014, gamma=-1.21,
               delta1=-1.7,delta2=0.86,r=-0.306,
               p1 = 0.48, p2 = 0.01, p3 = 0.83, p4 = 0.99,
               p5 = 0, p6 = 0.01, p7 = 0, p8 = 0, p9 = 0.03, p10 = 0,
               p11 = 0.01, p12 = 0.92, pH = 0.0045)



#ksAMG = 0.165 #[1/yr] #Default AMG1 is
ksAMG = c(k0=0.165,humABOVE=0.5,humBELOW=0.4)

#Starting date of the simulations
start_date_simulations_site=as.Date("2021-01-01")

#Read forcing
temp_df <- read.delim(paste0(loc_forc,"942217_TAS_2001_2005.txt"),sep = ":")[,2]
temp_df <- as.numeric(temp_df[7:length(temp_df)])-273.15 #K to C
temp_df <- c(temp_df,temp_df) #Repeat 2001-2005 twice


Temp_day <- data.frame("Date"=seq(from = as.Date("2001-01-01"), to = as.Date("2010-12-31"), by = 'day'),"Temp"=temp_df)
Temp_month <- aggregate(Temp_day, list(format(Temp_day$Date,"%Y-%m")),mean, na.rm = TRUE)
Temp_month$Date<-as.Date(paste(Temp_month$Group.1,"-01",sep=""))

rain_df <- read.delim(paste0(loc_forc,"942217_PR_2001_2005.txt"),sep = ":")[,2]
rain_df <- as.numeric(rain_df[7:length(rain_df)])# "kg m-2 s-1"
snow_df <- read.delim(paste0(loc_forc,"942217_PRSN_2001_2005.txt"),sep = ":")[,2]
snow_df <- as.numeric(snow_df[7:length(snow_df)]) #"kg m-2 s-1"
precip_df <- (rain_df+snow_df)*60*60*24 #"kg m-2 s-1" to mm/day

precip_df <- c(precip_df,precip_df) #Repeat 2001-2005 twice
Precip_day <- data.frame("Date"=seq(from = as.Date("2001-01-01"), to = as.Date("2010-12-31"), by = 'day'),"Precip"=precip_df)
Precip_month <- aggregate(Precip_day["Precip"], list(format(Precip_day$Date,"%Y-%m")),FUN=sum, na.rm = TRUE)
Precip_month$Date<-as.Date(paste(Precip_month$Group.1,"-01",sep=""))

#Create summer drought
Precip_month_drought <- Precip_month
Precip_month_drought$Precip <- with(Precip_month_drought, 
                                    ifelse(sub(".*-", "",Precip_month_drought["Group.1"][[1]])=="05"
                                           |sub(".*-", "",Precip_month_drought["Group.1"][[1]])=="06"
                                           |sub(".*-", "",Precip_month_drought["Group.1"][[1]])=="07", 0, Precip))


potevap_df <- read.delim(paste0(loc_forc,"942217_POT_2001_2005.txt"),sep = ":")[,2]
potevap_df <- as.numeric(potevap_df[7:length(potevap_df)]) 
potevap_df <-potevap_df*(60*60*24)*rep(c(31,28,31,30,31,30,31,31,30,31,30,31),length(potevap_df)/12)# "kg m-2 s-1" to mm/month

potevap_df <- c(potevap_df,potevap_df) #Repeat 2001-2005 twice
Potevap_month <- data.frame("Date"=seq(from = as.Date("2001-01-01"), to = as.Date("2010-12-31"), by = 'month'),"Potevap"=potevap_df)
#Potevap_month$Date<-as.Date(paste(Potevap_month$Group.1,"-01",sep=""))
