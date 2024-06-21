# This .R script tests if (TQf) can be used as a reliable predictor of the size
# of post(T) contributions to raw and reaction energies
#
#
library(ggplot2)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(data.table)
rm(list = ls())

options(max.print=100000)
au2kcal = 627.5186025


# THIS WILL CONTAIN SPECIES MISSIG DATA
setwd('/Users/37348458/thermo-boot/TQ_survey')
foo_DZ <- read.csv('data_pVDZ.csv')
foo_TZ <- read.csv('data_pVTZ.csv')

#take out h as it has no correlation and screws up some fits
foo_DZ <- foo_DZ[foo_DZ$Species != "h",]
foo_TZ <- foo_TZ[foo_TZ$Species != "h",]

foo_DZ$Species




####################################################################
# GRAB JUST THE DATA WE NEED
#

raw_DZ <- data.frame(Species = foo_DZ$Species)
raw_TZ <- data.frame(Species = foo_TZ$Species)

raw_DZ[["count_val_e"]] <- foo_DZ$count_val_e
raw_DZ[["count_H"]] <- foo_DZ$count_H
raw_DZ[["count_B"]] <- foo_DZ$count_B
raw_DZ[["count_C"]] <- foo_DZ$count_C
raw_DZ[["count_N"]] <- foo_DZ$count_N
raw_DZ[["count_O"]] <- foo_DZ$count_O
raw_DZ[["count_F"]] <- foo_DZ$count_F

raw_TZ[["count_val_e"]] <- foo_TZ$count_val_e
raw_TZ[["count_H"]] <- foo_TZ$count_H
raw_TZ[["count_B"]] <- foo_TZ$count_B
raw_TZ[["count_C"]] <- foo_TZ$count_C
raw_TZ[["count_N"]] <- foo_TZ$count_N
raw_TZ[["count_O"]] <- foo_TZ$count_O
raw_TZ[["count_F"]] <- foo_TZ$count_F

raw_DZ[["Max T1(A) : DZ"]] <- foo_DZ$Max.T1_A
raw_TZ[["Max T1(A) : TZ"]] <- foo_TZ$Max.T1_A
raw_DZ[["Max T2(AB) : DZ"]] <- foo_DZ$Max.T2_AB.CCSD
raw_TZ[["Max T2(AB) : TZ"]] <- foo_TZ$Max.T2_AB.CCSD

raw_DZ[["CCSD : DZ"]] <- foo_DZ$CCSD.cc.pVDZ
raw_DZ[["CCSD : DZ per e"]] <- foo_DZ$CCSD.cc.pVDZ / foo_DZ$count_val_e
raw_TZ[["CCSD : TZ"]] <- foo_TZ$CCSD.cc.pVTZ
raw_TZ[["CCSD : TZ per e"]] <- foo_TZ$CCSD.cc.pVTZ / foo_TZ$count_val_e

raw_DZ[["CCSD(T) : DZ"]] <- foo_DZ$CCSD.T..cc.pVDZ
raw_DZ[["CCSD(T) : DZ per e"]] <- foo_DZ$CCSD.T..cc.pVDZ / foo_DZ$count_val_e
raw_TZ[["CCSD(T) : TZ"]] <- foo_TZ$CCSD.T..cc.pVTZ
raw_TZ[["CCSD(T) : TZ per e"]] <- foo_TZ$CCSD.T..cc.pVTZ / foo_TZ$count_val_e

raw_DZ[["CCSD(T-5) : DZ"]] <- foo_DZ$CCSD.T.5..cc.pVDZ
raw_DZ[["CCSD(T-5) : DZ per e"]] <- foo_DZ$CCSD.T.5..cc.pVDZ / foo_DZ$count_val_e
raw_TZ[["CCSD(T-5) : TZ"]] <- foo_TZ$CCSD.T.5..cc.pVTZ
raw_TZ[["CCSD(T-5) : TZ per e"]] <- foo_TZ$CCSD.T.5..cc.pVTZ / foo_TZ$count_val_e

raw_DZ[["CCSDT : DZ"]] <- foo_DZ$CCSDT.cc.pVDZ
raw_DZ[["CCSDT : DZ per e"]] <- foo_DZ$CCSDT.cc.pVDZ / foo_DZ$count_val_e
raw_TZ[["CCSDT : TZ"]] <- foo_TZ$CCSDT.cc.pVTZ
raw_TZ[["CCSDT : TZ per e"]] <- foo_TZ$CCSDT.cc.pVTZ / foo_TZ$count_val_e


raw_DZ[["CCSD(TQ) : DZ"]] <- foo_DZ$CCSD.TQ..cc.pVDZ
raw_DZ[["CCSD(TQ) : DZ per e"]] <- foo_DZ$CCSD.TQ..cc.pVDZ / foo_DZ$count_val_e
#raw_TZ[["CCSD(TQ) : TZ"]] <- foo_TZ$CCSD.TQ..cc.pVTZ

raw_DZ[["CCSD(TQf) : DZ"]] <- foo_DZ$CCSD.TQf..cc.pVDZ
raw_DZ[["CCSD(TQf) : DZ per e"]] <- foo_DZ$CCSD.TQf..cc.pVDZ / foo_DZ$count_val_e
#raw_TZ[["CCSD(TQf) : TZ"]] <- foo_TZ$CCSD.TQf..cc.pVTZ

raw_DZ[["CCSDT(Q)L : DZ"]] <- foo_DZ$CCSDT.Q._L.cc.pVDZ
raw_DZ[["CCSDT(Q)L : DZ per e"]] <- foo_DZ$CCSDT.Q._L.cc.pVDZ/ foo_DZ$count_val_e
#raw_TZ[["CCSDT(Q)L : TZ"]] <- foo_TZ$CCSDT.Q._L.cc.pVTZ

#Drop empty values (NaNs)
raw_DZ <- na.omit(raw_DZ)
raw_TZ <- na.omit(raw_TZ)

#Keep only species in both sets
raw_DZ <- raw_DZ[(raw_DZ$Species %in% raw_TZ$Species), ]
raw_TZ <- raw_TZ[(raw_TZ$Species %in% raw_DZ$Species), ]

#And print
raw_DZ$Species
raw_TZ$Species


####################################################################
# ADD IN THE DATA WE NEED
#
#
# T <- (T) data

raw_DZ[["CCSD(T) - CCSD : DZ"]] <- raw_DZ[["CCSD(T) : DZ"]] - raw_DZ[["CCSD : DZ"]]
raw_DZ[["CCSD(T) - CCSD : DZ per e"]] <- raw_DZ[["CCSD(T) : DZ per e"]] - raw_DZ[["CCSD : DZ per e"]]
raw_TZ[["CCSD(T) - CCSD : TZ"]] <- raw_TZ[["CCSD(T) : TZ"]] - raw_TZ[["CCSD : TZ"]]
raw_TZ[["CCSD(T) - CCSD : TZ per e"]] <- raw_TZ[["CCSD(T) : TZ per e"]] - raw_TZ[["CCSD : TZ per e"]]


raw_DZ[["CCSDT - CCSD(T) : DZ"]] <- raw_DZ[["CCSDT : DZ"]] - raw_DZ[["CCSD(T) : DZ"]]
raw_DZ[["CCSDT - CCSD(T) : DZ per e"]] <- raw_DZ[["CCSDT : DZ per e"]] - raw_DZ[["CCSD(T) : DZ per e"]]
raw_TZ[["CCSDT - CCSD(T) : TZ"]] <- raw_TZ[["CCSDT : TZ"]] - raw_TZ[["CCSD(T) : TZ"]]
raw_TZ[["CCSDT - CCSD(T) : TZ per e"]] <- raw_TZ[["CCSDT : TZ per e"]] - raw_TZ[["CCSD(T) : TZ per e"]]


raw_DZ[["CCSD(T-5) - CCSD(T) : DZ"]] <- raw_DZ[["CCSD(T-5) : DZ"]] - raw_DZ[["CCSD(T) : DZ"]]
raw_DZ[["CCSD(T-5) - CCSD(T) : DZ per e"]] <- raw_DZ[["CCSD(T-5) : DZ per e"]] - raw_DZ[["CCSD(T) : DZ per e"]]
raw_TZ[["CCSD(T-5) - CCSD(T) : TZ"]] <- raw_TZ[["CCSD(T-5) : TZ"]] - raw_TZ[["CCSD(T) : TZ"]]
raw_TZ[["CCSD(T-5) - CCSD(T) : TZ per e"]] <- raw_TZ[["CCSD(T-5) : TZ per e"]] - raw_TZ[["CCSD(T) : TZ per e"]]

####################################################################
# INITIAL TESTS
#


#Scaning for outliers
par(mfrow = c(2,3))

plot(raw_DZ[["CCSDT - CCSD(T) : DZ"]])
plot(raw_DZ[["CCSD(T-5) - CCSD(T) : DZ"]])
plot(raw_DZ[["CCSD(T) - CCSD : DZ"]])

plot(raw_TZ[["CCSDT - CCSD(T) : TZ"]])
plot(raw_TZ[["CCSD(T-5) - CCSD(T) : TZ"]])
plot(raw_TZ[["CCSD(T) - CCSD : TZ"]])

raw_DZ$Species[which.max(raw_DZ[["CCSDT - CCSD(T) : DZ"]])]
raw_DZ[["CCSDT - CCSD(T) : DZ"]]
raw_DZ[["Species"]]

par(mfrow = c(2,3))

plot(raw_DZ[["CCSDT - CCSD(T) : DZ per e"]])
plot(raw_DZ[["CCSD(T-5) - CCSD(T) : DZ per e"]])
plot(raw_DZ[["CCSD(T) - CCSD : DZ per e"]])

plot(raw_TZ[["CCSDT - CCSD(T) : TZ per e"]])
plot(raw_TZ[["CCSD(T-5) - CCSD(T) : TZ per e"]])
plot(raw_TZ[["CCSD(T) - CCSD : TZ per e"]])

raw_DZ$Species[which.max(raw_DZ[["CCSDT - CCSD(T) : DZ per e"]])]
raw_DZ$Species[which.min(raw_DZ[["CCSDT - CCSD(T) : DZ per e"]])]

raw_DZ[["CCSDT - CCSD(T) : DZ per e"]]
raw_DZ[["Species"]]



#Guessing what might be important
par(mfrow = c(2,4))
plot(raw_DZ[["CCSD(T) - CCSD : DZ"]], raw_DZ[["CCSDT - CCSD(T) : DZ"]])
plot(raw_DZ[["CCSD(T-5) - CCSD(T) : DZ"]], raw_DZ[["CCSDT - CCSD(T) : DZ"]])
plot(abs(raw_DZ[["Max T1(A) : DZ"]]), raw_DZ[["CCSDT - CCSD(T) : DZ"]])
plot(abs(raw_DZ[["Max T2(AB) : DZ"]]), raw_DZ[["CCSDT - CCSD(T) : DZ"]])


plot(raw_TZ[["CCSD(T) - CCSD : TZ"]], raw_TZ[["CCSDT - CCSD(T) : TZ"]])
plot(raw_TZ[["CCSD(T-5) - CCSD(T) : TZ"]], raw_TZ[["CCSDT - CCSD(T) : TZ"]])
plot(abs(raw_TZ[["Max T1(A) : TZ"]]), raw_TZ[["CCSDT - CCSD(T) : TZ"]])
plot(abs(raw_TZ[["Max T2(AB) : TZ"]]), raw_TZ[["CCSDT - CCSD(T) : TZ"]])

raw_DZ$Species[which.max(raw_DZ[["CCSDT - CCSD(T) : DZ"]])]
raw_DZ[["CCSDT - CCSD(T) : DZ"]]
raw_DZ[["Species"]]



#Guessing what might be important
par(mfrow = c(2,2))
plot(raw_DZ[["CCSD(T) - CCSD : DZ per e"]], raw_DZ[["CCSDT - CCSD(T) : DZ per e"]])
plot(raw_DZ[["CCSD(T-5) - CCSD(T) : DZ per e"]], raw_DZ[["CCSDT - CCSD(T) : DZ per e"]])
#plot(abs(raw_DZ[["Max T1(A) : DZ"]]), raw_DZ[["CCSDT - CCSD(T) : DZ"]])
#plot(abs(raw_DZ[["Max T2(AB) : DZ"]]), raw_DZ[["CCSDT - CCSD(T) : DZ"]])


plot(raw_TZ[["CCSD(T) - CCSD : TZ per e"]], raw_TZ[["CCSDT - CCSD(T) : TZ per e"]])
plot(raw_TZ[["CCSD(T-5) - CCSD(T) : TZ per e"]], raw_TZ[["CCSDT - CCSD(T) : TZ per e"]])
#plot(abs(raw_TZ[["Max T1(A) : TZ"]]), raw_TZ[["CCSDT - CCSD(T) : TZ"]])
#plot(abs(raw_TZ[["Max T2(AB) : TZ"]]), raw_TZ[["CCSDT - CCSD(T) : TZ"]])



stop()

####################################################################
# Models we want to test
#

# 1. (T-5) - (T)
model_E_T_pt5mt_dz <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`CCSD(T) - CCSD`, data=raw, start=c(a=0, b=1))
summary(model_E_qL_tmD)

# 2. (TQf) - D 
model_E_qL_tqfmD <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`CCSD(TQf) - CCSD`, data=raw, start=c(a=0, b=1))
summary(model_E_qL_tqfmD)

# 3. [(T) - D ] / (T)
model_E_qL_tmDdt <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`[CCSD(T) - CCSD]/CCSD(T)`, data=raw, start=c(a=0, b=-1))
summary(model_E_qL_tmDdt)

# 4. [(TQf) - D ] / (TQf)
model_E_qL_tqfmDdtqf <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`[CCSD(TQf) - CCSD]/CCSD(TQf)`, data=raw, start=c(a=0, b=-1))
summary(model_E_qL_tqfmDdtqf)

# 5. (TQf) - (T)
model_E_qL_tqfmt <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`CCSD(TQf) - CCSD(T)`, data=raw, start=c(a=0, b=-1))
summary(model_E_qL_tqfmt)

# 6. [(TQf) - (T)] / (TQf)
model_E_qL_tqfmtdtqf <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`[CCSD(TQf) - CCSD(T)]/CCSD(TQf)`, data=raw, start=c(a=0, b=-1))
summary(model_E_qL_tqfmtdtqf)

# 7. [(TQf) - (T)] / D
model_E_qL_tqfmtdD <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`[CCSD(TQf) - CCSD(T)]/CCSD`, data=raw, start=c(a=0, b=-1))
summary(model_E_qL_tqfmtdD)


#raw_non2o4 <- raw[raw$Species != "n2o4",]
#model_E_qL_tqf_t <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`CCSD(TQf) - CCSD(T)`, data=raw, start=c(a=0, b=1))
#model_E_qL_tqf_t_non2o4 <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`CCSD(TQf) - CCSD(T)`, data=raw_non2o4, start=c(a=0, b=1))

#model_E_qL_tqf_t_t5_t <- nls(`CCSDT(Q)_L - CCSD(T)` ~ a + b*`CCSD(TQf) - CCSD(T)` + c*`CCSD(T-5) - CCSD(T)`, data=raw, start=c(a=0, b=1, c=1))
#summary(model_E_qL_tqf_t_t5_t)


raw[["HLC pred (T)-D"]] <- predict(model_E_qL_tmD, newdata=raw)
raw[["HLC pred (TQf)-D"]] <- predict(model_E_qL_tqfmD, newdata=raw)
raw[["HLC pred [(T)-D]/(T)"]] <- predict(model_E_qL_tmDdt, newdata=raw)
raw[["HLC pred [(TQf)-D]/(TQf)"]] <- predict(model_E_qL_tqfmDdtqf, newdata=raw)
raw[["HLC pred (TQf)-(T)"]] <- predict(model_E_qL_tqfmt, newdata=raw)
raw[["HLC pred [(TQf)-(T)]/(TQf)"]] <- predict(model_E_qL_tqfmtdtqf, newdata=raw)
raw[["HLC pred [(TQf)-(T)]/D"]] <- predict(model_E_qL_tqfmtdD, newdata=raw)








####################################################################
# split data into reference and non-reference
#Split noref and ref data frames
ref_strs <- c("h2", "b2h6", "ch4", "h2o", "nh3", "hf")
tf <- rep(FALSE, times = nrow(raw))
for (i in 1:nrow(raw))
{
  if (any(ref_strs %in% raw$Species[i]))
  {
    tf[i] <- TRUE
  }
}
raw$is_ref <- tf


data_noref <- raw[raw$is_ref == FALSE, ]
data_ref <- raw[raw$is_ref == TRUE, ]

#remove any NaNs
#data <- na.omit(data_noref)
#ref <- na.omit(data_ref)
data <- data_noref
ref <- data_ref

#double check that the references contain H2, CH4, NH3, and H2O
for (r in ref_strs)
{
  print(r)
  if (!any(ref$Species %in% r)) stop("Reference is missing")
}

####################################################################
# Greatest Common Denominator function for a list of integers
gcd <- function(a)
{
  m <- max(a)
  for (i in ceiling(sqrt(max(a))):2 )
  {
    g <- a %% i
    if (all(g ==0)) return (i)
  }
  return(1)
}

rms <- function(x)
{
  return (sqrt(mean(x^2)))
}

####################################################################
# GENERATE REACTIONS DATA
#
# This proceeds via two steps. The first step generates the reaction names 
# and the number of references (negative indicates LHS and positive indicates RHS)
#
# Once these are generates, the rest of the columns of the dataframe are filled. 


#we want only the pairs of species that share at least B, C, N, O, or F
rxn_names = c()
num_rct = c()
num_prd = c()
num_h2 = c()
num_b2h6 = c()
num_ch4 = c()
num_nh3 = c()
num_h2o = c()
num_hf = c()
species_lhs = c()
species_rhs = c()

#Pair reactions
for (i in 1:(nrow(data)-1))
{
  for (j in (i+1):nrow(data))
  {
      #screen all reactions without common heavy atom
      if (! ( (data$count_B[i]!=0 & data$count_B[j]!=0) 
           |  (data$count_C[i]!=0 & data$count_C[j]!=0) 
           |  (data$count_N[i]!=0 & data$count_N[j]!=0) 
           |  (data$count_O[i]!=0 & data$count_O[j]!=0)
           |  (data$count_F[i]!=0 & data$count_F[j]!=0) )
      ) next
  
      species_lhs = c(species_lhs, data$Species[i])
      species_rhs = c(species_rhs, data$Species[j])
      
      #reaction balancing #Original numbers
      n_rct  = 1
      n_prd  = 1
      n_h2   = 0  
      n_b2h6 = 0
      n_ch4  = 0
      n_nh3  = 0
      n_h2o  = 0
      n_hf   = 0
      delta_H = (n_prd*data$count_H[j] - n_rct*data$count_H[i]) + n_h2*2 + n_b2h6*6 + n_ch4*4 + n_nh3*3 + n_h2o*2 + n_hf*1
      delta_B = (n_prd*data$count_B[j] - n_rct*data$count_B[i]) + n_b2h6*2
      delta_C = (n_prd*data$count_C[j] - n_rct*data$count_C[i]) + n_ch4
      delta_N = (n_prd*data$count_N[j] - n_rct*data$count_N[i]) + n_nh3
      delta_O = (n_prd*data$count_O[j] - n_rct*data$count_O[i]) + n_h2o
      delta_F = (n_prd*data$count_F[j] - n_rct*data$count_F[i]) + n_hf
      #message(sprintf("%d, %d, %s -> %s : (%#.2f R, %#.2f P) [%#.2f H, %#.2f B, %#.2f C, %#.2f N, %#.2f O, %#.2f F] ", i, j, data$Species[i], data$Species[j], n_rct, n_prd, delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
      #message(sprintf("[%#.2f H2, %#.2f B2H6, %#.2f CH4, %#.2f NH3, %#.2f H2O, %#.2f HF] ", n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf))

      
      #first Balance out B2H6
      n_b2h6 = -delta_B
      if (abs(n_b2h6) > 0)
      {
        n_rct  = 2*n_rct
        n_prd  = 2*n_prd
        n_h2   = 2*n_h2
        n_ch4  = 2*n_ch4
        n_nh3  = 2*n_nh3
        n_h2o  = 2*n_h2o
        n_hf   = 2*n_hf
        delta_H = (n_prd*data$count_H[j] - n_rct*data$count_H[i]) + n_h2*2 + n_b2h6*6 + n_ch4*4 + n_nh3*3 + n_h2o*2 + n_hf*1
        delta_B = (n_prd*data$count_B[j] - n_rct*data$count_B[i]) + n_b2h6*2
        delta_C = (n_prd*data$count_C[j] - n_rct*data$count_C[i]) + n_ch4
        delta_N = (n_prd*data$count_N[j] - n_rct*data$count_N[i]) + n_nh3
        delta_O = (n_prd*data$count_O[j] - n_rct*data$count_O[i]) + n_h2o
        delta_F = (n_prd*data$count_F[j] - n_rct*data$count_F[i]) + n_hf
      }
      #message(sprintf("%d, %d, %s -> %s : (%#.2f R, %#.2f P) [%#.2f H, %#.2f B, %#.2f C, %#.2f N, %#.2f O, %#.2f F] ", i, j, data$Species[i], data$Species[j], n_rct, n_prd, delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
      #message(sprintf("[%#.2f H2, %#.2f B2H6, %#.2f CH4, %#.2f NH3, %#.2f H2O, %#.2f HF] ", n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf))
      
      
      #Now balance C,N,O,F
      n_ch4  = -delta_C
      n_nh3  = -delta_N
      n_h2o  = -delta_O
      n_hf   = -delta_F
      delta_H = (n_prd*data$count_H[j] - n_rct*data$count_H[i]) + n_h2*2 + n_b2h6*6 + n_ch4*4 + n_nh3*3 + n_h2o*2 + n_hf*1
      delta_B = (n_prd*data$count_B[j] - n_rct*data$count_B[i]) + n_b2h6*2
      delta_C = (n_prd*data$count_C[j] - n_rct*data$count_C[i]) + n_ch4
      delta_N = (n_prd*data$count_N[j] - n_rct*data$count_N[i]) + n_nh3
      delta_O = (n_prd*data$count_O[j] - n_rct*data$count_O[i]) + n_h2o
      delta_F = (n_prd*data$count_F[j] - n_rct*data$count_F[i]) + n_hf
      #message(sprintf("%d, %d, %s -> %s : (%#.2f R, %#.2f P) [%#.2f H, %#.2f B, %#.2f C, %#.2f N, %#.2f O, %#.2f F] ", i, j, data$Species[i], data$Species[j], n_rct, n_prd, delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
      #message(sprintf("[%#.2f H2, %#.2f B2H6, %#.2f CH4, %#.2f NH3, %#.2f H2O, %#.2f HF] ", n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf))
      
      
      #Finally h2
      n_h2 = -delta_H
      if (abs(n_h2) > 0)
      {
        n_rct  = 2*n_rct
        n_prd  = 2*n_prd
        n_b2h6 = 2*n_b2h6
        n_ch4  = 2*n_ch4
        n_nh3  = 2*n_nh3
        n_h2o  = 2*n_h2o
        n_hf   = 2*n_hf
        delta_H = (n_prd*data$count_H[j] - n_rct*data$count_H[i]) + n_h2*2 + n_b2h6*6 + n_ch4*4 + n_nh3*3 + n_h2o*2 + n_hf*1
        delta_B = (n_prd*data$count_B[j] - n_rct*data$count_B[i]) + n_b2h6*2
        delta_C = (n_prd*data$count_C[j] - n_rct*data$count_C[i]) + n_ch4
        delta_N = (n_prd*data$count_N[j] - n_rct*data$count_N[i]) + n_nh3
        delta_O = (n_prd*data$count_O[j] - n_rct*data$count_O[i]) + n_h2o
        delta_F = (n_prd*data$count_F[j] - n_rct*data$count_F[i]) + n_hf
      }
     # message(sprintf("%d, %d, %s -> %s : (%#.2f R, %#.2f P) [%#.2f H, %#.2f B, %#.2f C, %#.2f N, %#.2f O, %#.2f F] ", i, j, data$Species[i], data$Species[j], n_rct, n_prd, delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
     # message(sprintf("[%#.2f H2, %#.2f B2H6, %#.2f CH4, %#.2f NH3, %#.2f H2O, %#.2f HF] ", n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf))
      
      #Find largest common denominator
      coef <- c(n_rct, n_prd, n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf)
      coef <- abs(coef)
      m <- gcd(coef)
      n_rct  <- n_rct / m
      n_prd  <- n_prd / m
      n_h2   <- n_h2 / m
      n_b2h6 <- n_b2h6 / m
      n_ch4  <- n_ch4 / m
      n_nh3  <- n_nh3 / m
      n_h2o  <- n_h2o / m
      n_hf   <- n_hf / m
   
      #Generate string for reaction name
      lhs <- sprintf("(%d)%s", n_rct, data$Species[i])
      rhs <- sprintf("(%d)%s", n_prd, data$Species[j])
      if (n_h2   > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_h2,   "h2"), sep=" ")}
      if (n_h2   < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_h2,   "h2"), sep=" ")}
      if (n_b2h6 > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_b2h6, "b2h6"), sep=" ")}
      if (n_b2h6 < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_b2h6, "b2h6"), sep=" ")}
      if (n_ch4  > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_ch4,  "ch4"), sep=" ")}
      if (n_ch4  < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_ch4,  "ch4"), sep=" ")}
      if (n_nh3  > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_nh3,  "nh3"), sep=" ")}
      if (n_nh3  < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_nh3,  "nh3"), sep=" ")}
      if (n_h2o  > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_h2o,  "h2o"), sep=" ")}
      if (n_h2o  < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_h2o,  "h2o"), sep=" ")}
      if (n_hf   > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_hf,   "hf"), sep=" ")}
      if (n_hf   < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_hf,   "hf"), sep=" ")}
      
      #message(sprintf("%d, %d, %s -> %s : [%d, %d, %d, %d, %d, %d] ", i, j, data$Species[i], data$Species[j], delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))

      rxn_names <- c(rxn_names, paste(lhs, rhs, sep=" -> "))
      num_rct = c(num_rct, n_rct)
      num_prd = c(num_prd, n_prd)
      num_h2 = c(num_h2, n_h2)
      num_b2h6 = c(num_b2h6, n_b2h6)
      num_ch4 = c(num_ch4, n_ch4)
      num_nh3 = c(num_nh3, n_nh3)
      num_h2o = c(num_h2o, n_h2o)
      num_hf = c(num_hf, n_hf)
      
  } #loop over j in initialization
}#loop over i in initialization

#ANL reactions
for (i in 1:(nrow(data)))
{
    
  species_lhs = c(species_lhs, data$Species[i])
  species_rhs = c(species_rhs, "")
    
  #reaction balancing #Original numbers
  n_rct  = 1
  n_h2   = 0  
  n_b2h6 = 0
  n_ch4  = 0
  n_nh3  = 0
  n_h2o  = 0
  n_hf   = 0
  delta_H = (-n_rct*data$count_H[i]) + n_h2*2 + n_b2h6*6 + n_ch4*4 + n_nh3*3 + n_h2o*2 + n_hf*1
  delta_B = (-n_rct*data$count_B[i]) + n_b2h6*2
  delta_C = (-n_rct*data$count_C[i]) + n_ch4
  delta_N = (-n_rct*data$count_N[i]) + n_nh3
  delta_O = (-n_rct*data$count_O[i]) + n_h2o
  delta_F = (-n_rct*data$count_F[i]) + n_hf
  #message(sprintf("%d, %d, %s -> %s : (%#.2f R, %#.2f P) [%#.2f H, %#.2f B, %#.2f C, %#.2f N, %#.2f O, %#.2f F] ", i, j, data$Species[i], data$Species[j], n_rct, n_prd, delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
  #message(sprintf("[%#.2f H2, %#.2f B2H6, %#.2f CH4, %#.2f NH3, %#.2f H2O, %#.2f HF] ", n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf))
    
    
  #first Balance out B2H6
  n_b2h6 = -delta_B
  if (abs(n_b2h6) > 0)
  {
    n_rct  = 2*n_rct
    n_h2   = 2*n_h2
    n_ch4  = 2*n_ch4
    n_nh3  = 2*n_nh3
    n_h2o  = 2*n_h2o
    n_hf   = 2*n_hf
    delta_H = (-n_rct*data$count_H[i]) + n_h2*2 + n_b2h6*6 + n_ch4*4 + n_nh3*3 + n_h2o*2 + n_hf*1
    delta_B = (-n_rct*data$count_B[i]) + n_b2h6*2
    delta_C = (-n_rct*data$count_C[i]) + n_ch4
    delta_N = (-n_rct*data$count_N[i]) + n_nh3
    delta_O = (-n_rct*data$count_O[i]) + n_h2o
    delta_F = (-n_rct*data$count_F[i]) + n_hf
  }
  #message(sprintf("%d, %d, %s -> %s : (%#.2f R, %#.2f P) [%#.2f H, %#.2f B, %#.2f C, %#.2f N, %#.2f O, %#.2f F] ", i, j, data$Species[i], data$Species[j], n_rct, n_prd, delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
  #message(sprintf("[%#.2f H2, %#.2f B2H6, %#.2f CH4, %#.2f NH3, %#.2f H2O, %#.2f HF] ", n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf))
    
    
  #Now balance C,N,O,F
  n_ch4  = -delta_C
  n_nh3  = -delta_N
  n_h2o  = -delta_O
  n_hf   = -delta_F
  delta_H = (-n_rct*data$count_H[i]) + n_h2*2 + n_b2h6*6 + n_ch4*4 + n_nh3*3 + n_h2o*2 + n_hf*1
  delta_B = (-n_rct*data$count_B[i]) + n_b2h6*2
  delta_C = (-n_rct*data$count_C[i]) + n_ch4
  delta_N = (-n_rct*data$count_N[i]) + n_nh3
  delta_O = (-n_rct*data$count_O[i]) + n_h2o
  delta_F = (-n_rct*data$count_F[i]) + n_hf
  #message(sprintf("%d, %d, %s -> %s : (%#.2f R, %#.2f P) [%#.2f H, %#.2f B, %#.2f C, %#.2f N, %#.2f O, %#.2f F] ", i, j, data$Species[i], data$Species[j], n_rct, n_prd, delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
  #message(sprintf("[%#.2f H2, %#.2f B2H6, %#.2f CH4, %#.2f NH3, %#.2f H2O, %#.2f HF] ", n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf))
    
    
  #Finally h2
  n_h2 = -delta_H
  if (abs(n_h2) > 0)
  {
    n_rct  = 2*n_rct
    n_b2h6 = 2*n_b2h6
    n_ch4  = 2*n_ch4
    n_nh3  = 2*n_nh3
    n_h2o  = 2*n_h2o
    n_hf   = 2*n_hf
    delta_H = (-n_rct*data$count_H[i]) + n_h2*2 + n_b2h6*6 + n_ch4*4 + n_nh3*3 + n_h2o*2 + n_hf*1
    delta_B = (-n_rct*data$count_B[i]) + n_b2h6*2
    delta_C = (-n_rct*data$count_C[i]) + n_ch4
    delta_N = (-n_rct*data$count_N[i]) + n_nh3
    delta_O = (-n_rct*data$count_O[i]) + n_h2o
    delta_F = (-n_rct*data$count_F[i]) + n_hf
  }
  # message(sprintf("%d, %d, %s -> %s : (%#.2f R, %#.2f P) [%#.2f H, %#.2f B, %#.2f C, %#.2f N, %#.2f O, %#.2f F] ", i, j, data$Species[i], data$Species[j], n_rct, n_prd, delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
  # message(sprintf("[%#.2f H2, %#.2f B2H6, %#.2f CH4, %#.2f NH3, %#.2f H2O, %#.2f HF] ", n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf))
    
  #Find largest common denominator
  coef <- c(n_rct, n_h2, n_b2h6, n_ch4, n_nh3, n_h2o, n_hf)
  coef <- abs(coef)
  m <- gcd(coef)
  n_rct  <- n_rct / m
  n_h2   <- n_h2 / m
  n_b2h6 <- n_b2h6 / m
  n_ch4  <- n_ch4 / m
  n_nh3  <- n_nh3 / m
  n_h2o  <- n_h2o / m
  n_hf   <- n_hf / m
    
  #Generate string for reaction name
  lhs <- sprintf("(%d)%s", n_rct, data$Species[i])
  rhs <- sprintf("")
  if (n_h2   > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_h2,   "h2"), sep=" ")}
  if (n_h2   < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_h2,   "h2"), sep=" ")}
  if (n_b2h6 > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_b2h6, "b2h6"), sep=" ")}
  if (n_b2h6 < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_b2h6, "b2h6"), sep=" ")}
  if (n_ch4  > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_ch4,  "ch4"), sep=" ")}
  if (n_ch4  < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_ch4,  "ch4"), sep=" ")}
  if (n_nh3  > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_nh3,  "nh3"), sep=" ")}
  if (n_nh3  < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_nh3,  "nh3"), sep=" ")}
  if (n_h2o  > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_h2o,  "h2o"), sep=" ")}
  if (n_h2o  < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_h2o,  "h2o"), sep=" ")}
  if (n_hf   > 0) {rhs <- paste(rhs, sprintf("+ (%d)%s", n_hf,   "hf"), sep=" ")}
  if (n_hf   < 0) {lhs <- paste(lhs, sprintf("+ (%d)%s",-n_hf,   "hf"), sep=" ")}
  #message(sprintf("%d, %d, %s -> %s : [%d, %d, %d, %d, %d, %d] ", i, j, data$Species[i], data$Species[j], delta_H, delta_B, delta_C, delta_N, delta_O, delta_F))
    
  rxn_names <- c(rxn_names, paste(lhs, rhs, sep=" -> "))
  num_rct = c(num_rct, n_rct)
  num_prd = c(num_prd, 0)
  num_h2 = c(num_h2, n_h2)
  num_b2h6 = c(num_b2h6, n_b2h6)
  num_ch4 = c(num_ch4, n_ch4)
  num_nh3 = c(num_nh3, n_nh3)
  num_h2o = c(num_h2o, n_h2o)
  num_hf = c(num_hf, n_hf)
    
}#loop over i in initialization


rxns <- data.frame(rxn_name = rxn_names,
                   species_lhs = species_lhs,
                   species_rhs = species_rhs,
                   num_rct = num_rct, 
                   num_prd = num_prd,
                   num_h2 = num_h2,
                   num_b2h6 = num_b2h6,
                   num_ch4 = num_ch4, 
                   num_nh3 = num_nh3,
                   num_h2o = num_h2o, 
                   num_hf = num_hf)

print(rxns$rxn_name)

#Now fill in columns the for the reactions, but only the numeric cols that aren't integers
H2_rowid   <- which(ref$Species %in% "h2")
B2H6_rowid <- which(ref$Species %in% "b2h6")
CH4_rowid  <- which(ref$Species %in% "ch4")
NH3_rowid  <- which(ref$Species %in% "nh3")
H2O_rowid  <- which(ref$Species %in% "h2o")
HF_rowid   <- which(ref$Species %in% "hf")
for (col in names(data))
{
  if (!(is.numeric(data[[col]]) & !is.integer(data[[col]]))) next

  new_col <- c()
  for (row in 1:nrow(rxns))
  {
    
    #two types of reactions: either with LHS and RHS species or just LHS speces
    val = 0
    rct_rowid <- which(data$Species %in% rxns$species_lhs[row])
    
    if (rxns$num_prd[row] > 0)
    {
      prd_rowid <- which(data$Species %in% rxns$species_rhs[row])
      
      val <-  rxns$num_prd[row]  * data[[col]][prd_rowid] 
            - rxns$num_rct[row]  * data[[col]][rct_rowid]
            + rxns$num_h2[row]   * ref[[col]][H2_rowid]
            + rxns$num_b2h6[row] * ref[[col]][B2H6_rowid]
            + rxns$num_ch4[row]  * ref[[col]][CH4_rowid]
            + rxns$num_nh3[row]  * ref[[col]][NH3_rowid]
            + rxns$num_h2o[row]  * ref[[col]][H2O_rowid]
            + rxns$num_hf[row]   * ref[[col]][HF_rowid]
    }
    else
    {
      val <- -rxns$num_rct[row]  * data[[col]][rct_rowid]
            + rxns$num_h2[row]   * ref[[col]][H2_rowid]
            + rxns$num_b2h6[row] * ref[[col]][B2H6_rowid]
            + rxns$num_ch4[row]  * ref[[col]][CH4_rowid]
            + rxns$num_nh3[row]  * ref[[col]][NH3_rowid]
            + rxns$num_h2o[row]  * ref[[col]][H2O_rowid]
            + rxns$num_hf[row]   * ref[[col]][HF_rowid] 
    }
    
    new_col <- c(new_col, val)
    
  }#loop over rows of rxns
  
  #Add this column to reactions
  rxns[[col]] <- new_col
  
}#loop over cols to insert

#Add number of atoms/electrons involved
col_names <- c("count_H", "count_B", "count_C", "count_N", "count_O", "count_F", "count_val_e")
for (col in col_names)
{
  new_col <- c()
  
  for (row in 1:nrow(rxns))
  {
    val = 0
    rct_rowid <- which(data$Species %in% rxns$species_lhs[row])
  
    val <- val + abs(rxns$num_rct[row] * data[[col]][rct_rowid])
    if (rxns$num_h2[row]   < 0) val <- val + abs(rxns$num_h2[row]   * ref[[col]][H2_rowid])
    if (rxns$num_b2h6[row] < 0) val <- val + abs(rxns$num_b2h6[row] * ref[[col]][B2H6_rowid])
    if (rxns$num_ch4[row]  < 0) val <- val + abs(rxns$num_ch4[row]  * ref[[col]][CH4_rowid])
    if (rxns$num_nh3[row]  < 0) val <- val + abs(rxns$num_nh3[row]  * ref[[col]][NH3_rowid])
    if (rxns$num_h2o[row]  < 0) val <- val + abs(rxns$num_h2o[row]  * ref[[col]][H2O_rowid])
    if (rxns$num_hf[row]   < 0) val <- val + abs(rxns$num_hf[row]   * ref[[col]][HF_rowid])
    
    new_col <- c(new_col, val)
  
  }#loop over rows of rxns
  
  rxns[[col]] <- new_col
}#loop over subset of cols

################################################################################
# (Q)_L - (T) predictions

par(mfrow = c(4,4))

message(sprintf("RMS of (T) - D predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` - rxns$`HLC pred (T)-D`)*au2kcal))
plot(data$`CCSD(T) - CCSD`, data$`CCSDT(Q)_L - CCSD(T)`)
points(data$`CCSD(T) - CCSD`, data$`HLC pred (T)-D`, col='red')
plot(rxns$`CCSD(T) - CCSD`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
points(rxns$`CCSD(T) - CCSD`*au2kcal, rxns$`HLC pred (T)-D`*au2kcal, col='red')

message(sprintf("RMS of (TQf) - D predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` -rxns$`HLC pred (TQf)-D`)*au2kcal))
plot(data$`CCSD(TQf) - CCSD`, data$`CCSDT(Q)_L - CCSD(T)`)
points(data$`CCSD(TQf) - CCSD`, data$`HLC pred (TQf)-D`, col='red')
plot(rxns$`CCSD(TQf) - CCSD`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
points(rxns$`CCSD(TQf) - CCSD`*au2kcal, rxns$`HLC pred (TQf)-D`*au2kcal, col='red')

message(sprintf("RMS of [(T) - D]/(T) predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` -rxns$`HLC pred [(T)-D]/(T)`)*au2kcal))
plot(data$`[CCSD(T) - CCSD]/CCSD(T)`, data$`CCSDT(Q)_L - CCSD(T)`)
points(data$`[CCSD(T) - CCSD]/CCSD(T)`, data$`HLC pred [(T)-D]/(T)`, col='red')
plot(rxns$`[CCSD(T) - CCSD]/CCSD(T)`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
points(rxns$`[CCSD(T) - CCSD]/CCSD(T)`*au2kcal, rxns$`HLC pred [(T)-D]/(T)`*au2kcal, col='red')

message(sprintf("RMS of [(TQf) - D]/(TQf) predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` -rxns$`HLC pred [(TQf)-D]/(TQf)`)*au2kcal))
plot(data$`[CCSD(TQf) - CCSD]/CCSD(TQf)`, data$`CCSDT(Q)_L - CCSD(T)`)
points(data$`[CCSD(TQf) - CCSD]/CCSD(TQf)`, data$`HLC pred [(TQf)-D]/(TQf)`, col='red')
plot(rxns$`[CCSD(TQf) - CCSD]/CCSD(TQf)`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
points(rxns$`[CCSD(TQf) - CCSD]/CCSD(TQf)`*au2kcal, rxns$`HLC pred [(TQf)-D]/(TQf)`*au2kcal, col='red')


raw[["HLC pred (TQf)-(T)"]] <- predict(model_E_qL_tqfmt, newdata=raw)
message(sprintf("RMS of [(TQf) - (T)] predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` -rxns$`HLC pred (TQf)-(T)`)*au2kcal))
plot(data$`CCSD(TQf) - CCSD(T)`, data$`CCSDT(Q)_L - CCSD(T)`)
points(data$`CCSD(TQf) - CCSD(T)`, data$`HLC pred (TQf)-(T)`, col='red')
plot(rxns$`CCSD(TQf) - CCSD(T)`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
points(rxns$`CCSD(TQf) - CCSD(T)`*au2kcal, rxns$`HLC pred (TQf)-(T)`*au2kcal, col='red')


message(sprintf("RMS of [(TQf) - (T)]/(TQf) predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` -rxns$`HLC pred [(TQf)-(T)]/(TQf)`)*au2kcal))
plot(data$`[CCSD(TQf) - CCSD(T)]/CCSD(TQf)`, data$`CCSDT(Q)_L - CCSD(T)`)
points(data$`[CCSD(TQf) - CCSD(T)]/CCSD(TQf)`, data$`HLC pred [(TQf)-(T)]/(TQf)`, col='red')
plot(rxns$`[CCSD(TQf) - CCSD(T)]/CCSD(TQf)`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
points(rxns$`[CCSD(TQf) - CCSD(T)]/CCSD(TQf)`*au2kcal, rxns$`HLC pred [(TQf)-(T)]/(TQf)`*au2kcal, col='red')

message(sprintf("RMS of [(TQf) - (T)]/D predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` -rxns$`HLC pred [(TQf)-(T)]/D`)*au2kcal))
plot(data$`[CCSD(TQf) - CCSD(T)]/CCSD`, data$`CCSDT(Q)_L - CCSD(T)`)
points(data$`[CCSD(TQf) - CCSD(T)]/CCSD`, data$`HLC pred [(TQf)-(T)]/D`, col='red')
plot(rxns$`[CCSD(TQf) - CCSD(T)]/CCSD`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
points(rxns$`[CCSD(TQf) - CCSD(T)]/CCSD`*au2kcal, rxns$`HLC pred [(TQf)-(T)]/D`*au2kcal, col='red')


#message(sprintf("RMS of (TQf) - (T) + (T-5) - (T)predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` - rxns$`HLC pred (TQf)-(T) + (T-5)-(T)`)*au2kcal))
#plot(data$`CCSD(TQf) - CCSD(T)`, data$`CCSDT(Q)_L - CCSD(T)`)
#points(data$`CCSD(TQf) - CCSD(T)`, data$`HLC pred (TQf)-(T) + (T-5)-(T)`, col='red')
#plot(rxns$`CCSD(TQf) - CCSD(T)`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
#points(rxns$`CCSD(TQf) - CCSD(T)`*au2kcal, rxns$`HLC pred (TQf)-(T)`*au2kcal, col='red')

#message(sprintf("RMS of (TQf) - (T) no n2o4 predictor is %f", rms(rxns$`CCSDT(Q)_L - CCSD(T)` - rxns$`HLC pred (TQf)-(T) no n2o4`)*au2kcal))
#plot(data$`CCSD(TQf) - CCSD(T)`, data$`CCSDT(Q)_L - CCSD(T)`)
#points(data$`CCSD(TQf) - CCSD(T)`, data$`HLC pred (TQf)-(T) no n2o4`, col='red')
#plot(rxns$`CCSD(TQf) - CCSD(T)`*au2kcal, rxns$`CCSDT(Q)_L - CCSD(T)`*au2kcal)
#points(rxns$`CCSD(TQf) - CCSD(T)`*au2kcal, rxns$`HLC pred (TQf)-(T) no n2o4`*au2kcal, col='red')





