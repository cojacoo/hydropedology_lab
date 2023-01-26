setwd('~/Documents/TUBAF/teaching/hydropedology/uebung/hydropedology_lab/modelling_project/')
library(airGRteaching)
library(dplyr)
library(zoo)
library(lubridate)

carlsfeld <- read.csv('WBcarlsfeld1D.csv', header = TRUE, stringsAsFactors = FALSE)
carlsfeld$X <- as.POSIXct(carlsfeld$X, tz = 'UTC')
carlsfeld$ETo[carlsfeld$ETo<0.] <- 0.
carlsfeld$BachOst[carlsfeld$BachOst<0.0001] <- NA
carlsfeld$Wilzsch[carlsfeld$Wilzsch<0.0001] <- NA
carlsfeld$KleineBockau[carlsfeld$KleineBockau<0.0001] <- NA
carlsfeld <- carlsfeld %>% mutate(BachOst = na.approx(BachOst))
carlsfeld <- carlsfeld %>% mutate(Wilzsch = na.approx(Wilzsch))
carlsfeld <- carlsfeld %>% mutate(KleineBockau = na.approx(KleineBockau))
carlsfeld <- carlsfeld %>% mutate(T = na.approx(T))

ShinyGR(ObsDF = list("BachOst" = carlsfeld[, c('X', 'Prec', 'ETo', 'BachOst', 'T')], "Wilzsch" = carlsfeld[, c('X', 'Prec', 'ETo', 'Wilzsch', 'T')], "KleineBockau" = carlsfeld[, c('X', 'Prec', 'ETo', 'KleineBockau', 'T')]),
                     SimPer = c('2011-10-01', '2021-09-30'))



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
