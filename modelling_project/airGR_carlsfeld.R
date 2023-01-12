setwd('~/Documents/TUBAF/teaching/hydropedology/uebung/hydropedology_lab/modelling_project/')
library(airGRteaching)
library(lubridate)

carlsfeld <- read.csv('WBcarlsfeld1D.csv', header = TRUE, stringsAsFactors = FALSE)
carlsfeld$X <- as.POSIXct(carlsfeld$X, tz = "UTC")

ShinyGR(ObsDF = carlsfeld[1096:2253, c('X', 'Prec', 'ETo', 'BachOst', 'T')], SimPer = c("2014-10-01", "2017-10-30"))
ShinyGR(ObsDF = carlsfeld[731:1860, c('X', 'Prec', 'ETo', 'Wilzsch', 'T')], SimPer = c("2013-10-01", "2016-10-30"))
#ShinyGR(ObsDF = carlsfeld[, c('X', 'Prec', 'ETo', 'BachOst', 'T')], SimPer = c("2018-01-08", "2021-09-30"))



# Data processing for GR4J (with Q for calibration)
prep <- PrepGR(DatesR     = carlsfeld$X,
               Precip     = carlsfeld$Prec,
               PotEvap    = carlsfeld$ETo, 
               TempMean   = carlsfeld$T, 
               Qobs       = carlsfeld$BachOst,
               HydroModel = "GR4J", 
               CemaNeige  = FALSE)

# Parameter set to test
i_param_gr4j <- c(X1 = 350, X2 = 0, X3 = 90, X4 = 1.4)

# Rainfall-runoff simulation on the calibration period
i_sim_manu <- SimGR(PrepGR  = prep, 
                    Param   = i_param_gr4j,
                    WupPer  = c("2018-01-08", "2018-04-30"),
                    SimPer  = c("2018-05-01", "2021-09-30"),
                    EffCrit = "NSE",
                    verbose = TRUE)

# Calibration using NSE score 
cal_auto <- CalGR(PrepGR  = prep, 
                  CalCrit = "KGE",
                  WupPer  = c("2018-01-08", "2018-04-30"),
                  CalPer  = c("2018-05-01", "2021-09-30"))
plot(cal_auto)
