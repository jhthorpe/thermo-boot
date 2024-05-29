# This version of the script is designed to include the information about
# Lewis structures and bonding. 
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

setwd('/Users/37348458/thermo-boot/ANL/CBS')

#select species for which ATcT has high confidence
ANL <- read.csv('ANL-CBS.csv')
ANL_clean <- subset(ANL, ANL$ATcT_unc < 0.25)

#Split noref and ref data frames
ref_strs <- c("H2", "CH4", "NH3", "H2O")
tf <- rep(FALSE, times = nrow(ANL_clean))
for (i in 1:nrow(ANL_clean))
{
  if (any(ref_strs %in% ANL_clean$Species[i]))
  {
    tf[i] <- TRUE
  }
}
ANL_clean$is_ref <- tf
ANL_clean_noref <- ANL_clean[ANL_clean$is_ref == FALSE, ]
ANL_clean_ref <- ANL_clean[ANL_clean$is_ref == TRUE, ]

#remove any NaNs
data <- na.omit(ANL_clean_noref)
ref <- na.omit(ANL_clean_ref)

#double check that the references contain H2, CH4, NH3, and H2O
for (r in ref_strs)
{
  print(r)
  if (!any(ref$Species %in% r)) stop("Reference is missing")
}


#Function that generates delta of a given column in the ANL reaction scheme
delta_anl_rxn <- function(noref, ref, col)
{
  H2_rowid <- which(ref$Species %in% "H2")
  CH4_rowid <- which(ref$Species %in% "CH4")
  NH3_rowid <- which(ref$Species %in% "NH3")
  H2O_rowid <- which(ref$Species %in% "H2O")
  return(
    ref[[col]][H2_rowid] * noref$X.H2
    + ref[[col]][CH4_rowid] * noref$X.C
    + ref[[col]][NH3_rowid] * noref$X.N
    + ref[[col]][H2O_rowid] * noref$X.O
    - noref[[col]]
  )
}



#Add in deltas for all the lewis strings
data$rxn_total <- rep(0, times = nrow(data))
lewis_strs <- c("lp_C", "lp_N", "lp_O", "r_H", "r_C", "r_N", "r_O", "s_H.H", "s_C.H", "s_N.H", "s_O.H", "s_C.C", "s_C.N", "s_C.O", "s_N.N", "s_N.O", "s_O.O", "d_C.C", "d_C.N", "d_C.O", "d_N.N", "d_N.O", "d_O.O", "t_C.C", "t_C.N", "t_C.O", "t_N.N", "t_N.O")
for (str in lewis_strs)
{
  col_name <- paste("X.", str, sep="")
  rxn_name <- paste("rxn_", str, sep="")
  
  #one electron changes
  if (grepl("r_", str, fixed=TRUE))
  {
    data[[rxn_name]] <- abs(delta_anl_rxn(data, ref, col_name))
    data$rxn_total <- data$rxn_total + data[[rxn_name]]
  }
  #two electron changes
  else
  {
    data[[rxn_name]] <- abs(delta_anl_rxn(data, ref, col_name))*2
    data$rxn_total <- data$rxn_total + data[[rxn_name]]*2
  }
}





#Rename for my sanity
data <- data %>% rename(species = Species,
                        ATcT_val = ATcT,
                        ATcT_unc = ATcT_unc,
                        CCSDt_CBS2 = RXN_t_a.Q.5.Z.au.,
                        CCSDt_CBS1 = RXN_t_a.T.Q.Z.au.,
                        CCSDt_av5z = RXN_t_a5Z,
                        CCSDt_avqz = RXN_t_aQZ,
                        CCSDt_avtz = RXN_t_aTZ
                        )



###########################################################
# HELPER FUNCTIONS

#example of a function definition
apbx <- function(x, a, b) {
  a + b*x
}

#Poly3
poly_2 <- function(x, a, b, c) {
  a + b*x + c*x^2
}

#Pade approximants
# x   -> vector of x values
# c   -> pade coefficients. 1-3 for numerator, 4 for denominator
pade_2_1 <- function(x, n0, n1, n2, d1)
{
    (n0 + n1*x + n2*x^2)/(1 + d1*x)
}

#Pade approximants
# x   -> vector of x values
# c   -> pade coefficients. 1-3 for numerator, 4 for denominator
pade_3_2 <- function(x, n0, n1, n2, n3, d1, d2)
{
  (n0 + n1*x + n2*x^2 + n3*x^3)/(1 + d1*x + d2*x^2)
}


# Heaviside function
heaviside <- function(x, a=0)
{
  ifelse(x > a, 1, 0)
}

#Testing nonlinear coverage optimization : predictor
# linear fit
#
# alpha -> current set of parameters
# x     -> possibly multidimensional prediction values
prd <- function(alpha, x)
{
  alpha[1] + alpha[2]*x
}

#Testing nonlinear coverage optimization : widthor
#linear fit, must be abs as well
#
# alpha -> current set of parameters
# x     -> possibly multidimensional prediction values
wdth <- function(alpha, x)
{
    abs(alpha[1] + alpha[2]*x)
}

#Coverage function
# returns 1 if the observed data, y, is covered by our predictor and widthor,
# and 0 otherwise. 
#
# predictor (prd) and widthor (wdth) parametric values contained in alpha, 
# user needs to know which 
# variables in alpha correspond to the prd and wdth. 
#
# alpha -> current set of parameters
# x     -> possibly multidimensional prediction values
# y     -> single observed datapoint
covr <- function(alpha, x, y)
{
    1 - heaviside(abs(y - prd(alpha[1:2], x)), a = 2*wdth(alpha[3:4], x))
}

#fractional coverage function
# returns the fraction of data values currently covered by the predictor and
# weightor
# alpha -> current set of parameters
# df_x  -> dataframe of x values, possibly multidimensional
# df_y  -> 
frc_covr <- function(alpha, df_x, df_y)
{
   val = 0.
   for (i in 1:length(df_y))
   {
     val <- val + covr(alpha, df_x[i], df_y[i]) 
   }
   val <- val / length(df_y)
   val
}

#
#weight function for optimization
# alpha -> current set of parameters
# df_x  -> dataframe of x values, possibly multidimensional
# df_y  -> 
wgt <- function(alpha, df_x, df_y)
{
    (0.95 - frc_cov(alpha, df_x, df_y))**2.   
}

#
# True CBS1 err minus model from nls
err <- function(df, md)
{
  df$CCSDt_CBS1_err - predict(md)
}

# True CBS1 err minus model from nls divided by true error
frac_err <- function(df, md)
{
  (df$CCSDt_CBS1_err - predict(md))/df$CCSDt_CBS1
}

# Multiplies by a Normal Distributed squeezer, normalized to 1
squeez <- function(x, mu, sig)
{
  dnorm(x, mu, sig)/dnorm(mu, mu, sig)
}

# Multiplies by a Normal Distributed stretcher, normalized to 1
stretch <- function(x, mu, sig)
{
  1 - squeez(x, mu, sig)
}

###########################################################
# GENERATE AUXILLARY DATA
#
# The goal here is to try and get a decent model for 
# the CBS error of {T,Q} vs {Q,5} given the distance from
# Q to {T,Q}.
#
# 

#Generate CBS error data
data$CCSDt_CBS1_err <- (data$CCSDt_CBS2 -  data$CCSDt_CBS1) 
data$CCSDt_CBS1_extdist_q <- (data$CCSDt_CBS1 -  data$CCSDt_avqz)
data$CCSDt_CBS1_extdist_t <- (data$CCSDt_CBS1 -  data$CCSDt_avtz)
data$CCSDt_CBS1_frac_err <- (data$CCSDt_CBS2 -  data$CCSDt_CBS1)/(data$CCSDt_CBS1)
data$CCSDt_CBS1_frac_extdist_q <- (data$CCSDt_CBS1 -  data$CCSDt_avqz)/(data$CCSDt_CBS1)
data$CCSDt_CBS1_frac_extdist_t <- (data$CCSDt_CBS1 -  data$CCSDt_avtz)/(data$CCSDt_CBS1)
data$CCSDt_CBS1_dist_tq <- (data$CCSDt_avqz -  data$CCSDt_avtz)
data$CCSDt_CBS1_frac_dist_tq <- (data$CCSDt_avqz -  data$CCSDt_avtz)/(data$CCSDt_CBS1)

data$CCSDt_CBS1_err_per_e <- (data$CCSDt_CBS2 -  data$CCSDt_CBS1) / data$X.all_e
data$CCSDt_CBS1_extdist_q_per_e <- (data$CCSDt_CBS1 -  data$CCSDt_avqz) / data$X.all_e
data$CCSDt_CBS1_extdist_t_per_e <- (data$CCSDt_CBS1 -  data$CCSDt_avtz) / data$X.all_e
data$CCSDt_CBS1_frac_err_per_e <- ((data$CCSDt_CBS2 -  data$CCSDt_CBS1)/(data$CCSDt_CBS1)) / data$X.all_e
data$CCSDt_CBS1_frac_extdist_q_per_e <- ((data$CCSDt_CBS1 -  data$CCSDt_avqz)/(data$CCSDt_CBS1)) / data$X.all_e
data$CCSDt_CBS1_frac_extdist_t_per_e <- ((data$CCSDt_CBS1 -  data$CCSDt_avtz)/(data$CCSDt_CBS1)) / data$X.all_e
data$CCSDt_CBS1_dist_tq_per_e <- ((data$CCSDt_avqz -  data$CCSDt_avtz)) / data$X.all_e
data$CCSDt_CBS1_frac_dist_tq_per_e <- ((data$CCSDt_avqz -  data$CCSDt_avtz)/(data$CCSDt_CBS1)) / data$X.all_e

par(mfrow = c(1,1))
hist(data$CCSDt_CBS1_err)
#hist(data$CCSDt_CBS1_err_per_e)
     



par(mfrow = c(2,3))
plot(data$CCSDt_CBS1_frac_err, data$CCSDt_CBS1_err)
plot(data$X.all_e, data$CCSDt_CBS1_err)
plot(data$rxn_total, data$CCSDt_CBS1_err)

plot(data$CCSDt_CBS1_frac_err, data$CCSDt_CBS1_err_per_e)
plot(data$X.all_e, data$CCSDt_CBS1_err_per_e)
plot(data$rxn_total, data$CCSDt_CBS1_err_per_e)




###############################################################################
# DISTRIBUTIONS
###############################################################################

mu_extdist_q <- mean(data$CCSDt_CBS1_extdist_q)
mu_extdist_t <- mean(data$CCSDt_CBS1_extdist_t)
mu_CBS1 <- mean(data$CCSDt_CBS1)
mu_q <- mean(data$CCSDt_avqz)
mu_t <- mean(data$CCSDt_avtz)
mu_frac_extdist_q <- mean(data$CCSDt_CBS1_frac_extdist_q)
mu_frac_extdist_t <- mean(data$CCSDt_CBS1_frac_extdist_t)

sig_extdist_q <- sd(data$CCSDt_CBS1_extdist_q)
sig_extdist_t <- sd(data$CCSDt_CBS1_extdist_t)
sig_CBS1 <- sd(data$CCSDt_CBS1)
sig_q <- sd(data$CCSDt_avqz)
sig_t <- sd(data$CCSDt_avtz)
sig_frac_extdist_q <- sd(data$CCSDt_CBS1_frac_extdist_q)
sig_frac_extdist_t <- sd(data$CCSDt_CBS1_frac_extdist_t)


###############################################################################
# STEP 1: MEAN ERROR
###############################################################################
mu_CBS1_err <- mean(data$CCSDt_CBS1_err)
sd_CBS1_err <- sd(data$CCSDt_CBS1_err - mu_CBS1_err)
amax_CBS1_err <- max(abs(data$CCSDt_CBS1_err - mu_CBS1_err))

mu_CBS1_err_per_e <- mean((data$CCSDt_CBS1_err- mu_CBS1_err)/ data$X.all_e)
sd_CBS1_err_per_e <- sd((data$CCSDt_CBS1_err - mu_CBS1_err)/ data$X.all_e)
amax_CBS1_err_per_e <- max(abs((data$CCSDt_CBS1_err - mu_CBS1_err)/ data$X.all_e))

mean(abs(data$CCSDt_CBS1_err - mu_CBS1_err))
print(sd_CBS1_err)
print(amax_CBS1_err)

print(mean(abs(data$CCSDt_CBS1_err - mu_CBS1_err)/data$X.all_e))
print(sd_CBS1_err_per_e)
print(amax_CBS1_err_per_e)

mean(abs(data$CCSDt_CBS1_err - mu_CBS1_err))
sqrt(mean((data$CCSDt_CBS1_err - mu_CBS1_err)^2))
max(abs(data$CCSDt_CBS1_err - mu_CBS1_err))

mean(abs((data$CCSDt_CBS1_err - mu_CBS1_err) / data$X.all_e))
sqrt(mean(((data$CCSDt_CBS1_err - mu_CBS1_err) / data$X.all_e)^2)) 
max(abs((data$CCSDt_CBS1_err - mu_CBS1_err) / data$X.all_e)) 

par(mfrow = c(2,1))
hist(data$CCSDt_CBS1_err - mu_CBS1_err)
hist((data$CCSDt_CBS1_err - mu_CBS1_err) / data$X.all_e)




print("END OF STEP 1: MEAN")

par(mfrow = c(3,2))
plot(data$X.all_e, data$CCSDt_CBS1_err - mu_CBS1_err)
plot(data$X.all_e, (data$CCSDt_CBS1_err - mu_CBS1_err)/data$X.all_e)

#plot(data$rxn_total, data$CCSDt_CBS1_err - mu_CBS1_err)
#plot(data$rxn_total, (data$CCSDt_CBS1_err - mu_CBS1_err)/data$X.all_e)

plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err - mu_CBS1_err)
plot(data$CCSDt_CBS1, (data$CCSDt_CBS1_err - mu_CBS1_err)/data$X.all_e)

plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err - mu_CBS1_err)
plot(data$CCSDt_CBS1_extdist_q, (data$CCSDt_CBS1_err - mu_CBS1_err)/data$X.all_e)



###########################################################################
# STEP 2 : Fit number of atoms
# WTF why is this better....
###########################################################################

#constant plus num e plus num e pair
func <- function(df, a, a_H=0, a_C=0, a_N=0, a_O=0)
{
   return(
     a 
     + (a_H)*df$X.H + (a_C)*df$X.C*4 + (a_N)*df$X.N*5 + (a_O)*df$X.O*6
   )
   
}

model <- nls(data$CCSDt_CBS1_err ~ func(data, a, a_H, a_C, a_N, a_O), data=data, start=c(a=-0.51,a_H=0.1, a_C=0.1, a_N=0.1, a_O=0.1))
summary(model)
data$model <- predict(model)


#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(err(data, model) / data$X.all_e))
sqrt(mean((err(data, model) / data$X.all_e)^2)) 
sd(err(data, model) / data$X.all_e) 
max(abs(err(data, model) / data$X.all_e)) 


par(mfrow = c(3,2))

plot(data$X.all_e, data$CCSDt_CBS1_err)
points(data$X.all_e, data$model, type='p', col='red')
plot(data$X.all_e, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1, data$model, type='p', col='red')
plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1_extdist_q, data$model, type='p', col='red')
plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err - data$model)


print("END OF 2")

###########################################################################
# STEP 2 : Fit number of atoms
# WTF why is this better....
###########################################################################

#constant plus num e plus num e pair
func <- function(df, a,a_C=0, a_N=0, a_O=0)
{
  return(
    a 
    + (a_C)*df$X.C*4*3/2 + (a_N)*df$X.N*5*4/2 + (a_O)*df$X.O*6*5/2
  )
  
}

model <- nls(data$CCSDt_CBS1_err ~ func(data, a, a_C, a_N, a_O), data=data, start=c(a=-0.51, a_C=0.1, a_N=0.1, a_O=0.1))
summary(model)
data$model <- predict(model)


#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(err(data, model) / data$X.all_e))
sqrt(mean((err(data, model) / data$X.all_e)^2)) 
sd(err(data, model) / data$X.all_e) 
max(abs(err(data, model) / data$X.all_e)) 


par(mfrow = c(3,2))

plot(data$X.all_e, data$CCSDt_CBS1_err)
points(data$X.all_e, data$model, type='p', col='red')
plot(data$X.all_e, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1, data$model, type='p', col='red')
plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1_extdist_q, data$model, type='p', col='red')
plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err - data$model)


print("END OF 2: JUST ATOMS")

###########################################################################
# STEP 3 : lone-pair plus radical e
###########################################################################

#constant plus num e plus num e pair
func <- function(df, a, lp_C=0, lp_N=0, lp_O=0, r_H=0, r_C=0, r_N=0, r_O=0, s_HH=0, s_CH=0, s_NH=0, s_OH=0, s_CC=0, s_CN=0, s_CO=0, s_NN=0, s_NO=0, s_OO=0, d_CC=0, d_CN=0, d_CO=0, d_NN=0, d_NO=0, d_OO=0, t_CC=0, t_CN=0, t_CO=0, t_NN=0, t_NO=0)
{
  return(
    a 
    + lp_C*df$X.lp_C + lp_N*df$X.lp_N + lp_O*df$X.lp_O
    + r_H*df$X.r_H + r_C*df$X.r_C + r_N*df$X.r_N + r_O*df$X.r_O 
    
  #  a 
  #  + lp_C*df$rxn_lp_C + lp_N*df$rxn_lp_N + lp_O*df$rxn_lp_O
  #  + r_H*df$rxn_r_H + r_C*df$rxn_r_C + r_N*df$rxn_r_N + r_O*df$rxn_r_O 
  )
}

model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O), data=data, start=c(a=-0.51, lp_C=0.1, lp_N=0.1, lp_O=0.1, r_H=0.1, r_C=0.1, r_N=0.1, r_O=0.1))
summary(model)
data$model <- predict(model)

#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(err(data, model) / data$X.all_e))
sqrt(mean((err(data, model) / data$X.all_e)^2)) 
sd(err(data, model) / data$X.all_e) 
max(abs(err(data, model) / data$X.all_e)) 


par(mfrow = c(3,2))

plot(data$X.all_e, data$CCSDt_CBS1_err)
points(data$X.all_e, data$model, type='p', col='red')
plot(data$X.all_e, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1, data$model, type='p', col='red')
plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1_extdist_q, data$model, type='p', col='red')
plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err - data$model)

###########################################################################
# STEP 3 : lone-pair/radical e/singles
###########################################################################

print(data$rxn_s_H.H)

print(data$rxn_s_C.H)
print(data$rxn_s_N.H)
print(data$rxn_s_O.H)
#constant plus num e plus num e pair
func <- function(df, a, lp_C=0, lp_N=0, lp_O=0, r_H=0, r_C=0, r_N=0, r_O=0, s_HH=0, s_CH=0, s_NH=0, s_OH=0, s_CC=0, s_CN=0, s_CO=0, s_NN=0, s_NO=0, s_OO=0, d_CC=0, d_CN=0, d_CO=0, d_NN=0, d_NO=0, d_OO=0, t_CC=0, t_CN=0, t_CO=0, t_NN=0, t_NO=0)
{
  return(
    #a 
    #+ lp_C*df$X.lp_C + lp_N*df$X.lp_N + lp_O*df$X.lp_O
    #+ r_H*df$X.r_H + r_C*df$X.r_C + r_N*df$X.r_N + r_O*df$X.r_O 
    #+ s_HH*df$X.s_H.H + s_CH*df$X.s_C.H + s_NH*df$X.s_N.H + s_OH*df$X.s_O.H
    #+ s_CC*df$X.s_C.C + s_CN*df$X.s_C.N + s_CO*df$X.s_C.O 
    #+ s_NN*df$X.s_N.N + s_NO*df$X.s_N.O 
    #+ s_OO*df$X.s_O.O
    #
      a 
      + lp_C*df$rxn_lp_C + lp_N*df$rxn_lp_N + lp_O*df$rxn_lp_O
      + r_H*df$rxn_r_H + r_C*df$rxn_r_C + r_N*df$rxn_r_N + r_O*df$rxn_r_O 
      + s_HH*df$rxn_s_H.H + s_CH*df$rxn_s_C.H + s_NH*df$rxn_s_N.H + s_OH*df$rxn_s_O.H
      + s_CC*df$rxn_s_C.C + s_CN*df$rxn_s_C.N + s_CO*df$rxn_s_C.O 
      + s_NN*df$rxn_s_N.N + s_NO*df$rxn_s_N.O 
      + s_OO*df$rxn_s_O.O
      
  )
}
#This one works, but seems that full set of X-H contains redundent info
model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, 0), data=data, start=c(a=0.174355, lp_C=-0.073996 , lp_N=-0.024296, lp_O=-0.038120, r_H=-0.347418, r_C=-0.051776, r_N=-0.074486, r_O=-0.087346, s_HH=0.185226, s_CH=-0.075798, s_NH=-0.070630))

#Adding in C-X now
model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, 0, s_CC, s_CN, s_CO), data=data, start=c(a=0.088055, lp_C=-0.018452 , lp_N=-0.024296, lp_O=-0.038120, r_H=-0.265217, r_C=-0.051776, r_N=-0.074486, r_O=-0.087346, s_HH=0.189326, s_CH=-0.083977, s_NH=-0.063486, s_CC=0.046611, s_CN=0.033847, s_CO=0.045068))

#Adding in N-X now
model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, 0, s_CC, s_CN, s_CO, s_NN, s_NO), data=data, start=c(a=0.101722, lp_C=-0.018452 , lp_N=-0.050854, lp_O=-0.046037, r_H=-0.265217, r_C=-0.017631, r_N=-0.010248, r_O=-0.072045, s_HH=0.189326, s_CH=-0.083977, s_NH=-0.063486, s_CC=0.049026, s_CN=0.033847, s_CO=0.052968, s_NN=0.035128, s_NO=0.085973))

#Adding in O-X now, and replaced H-H with O-H
model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, 0, s_CH, s_NH, s_OH, s_CC, s_CN, s_CO, s_NN, s_NO, s_OO), data=data, start=c(a=0.064089, lp_C=-0.024744 , lp_N=-0.045833, lp_O=-0.003477, r_H=-0.183140, r_C=-0.010248, r_N=-0.022725, r_O=-0.068088, s_CH=-0.050479, s_NH=-0.036269, s_OH=0.065607, s_CC=0.047337, s_CN=0.038151, s_CO=0.077671, s_NN=0.042359, s_NO=0.120844, s_OO=0.165714))



#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH), data=data, start=c(a=0.214607, lp_C=-0.073996 , lp_N=-0.024296, lp_O=-0.038120, r_H=-0.347418, r_C=-0.051776, r_N=-0.074486, r_O=-0.087346, s_HH=0.185226, s_CH=-0.075798, s_NH=-0.070630, s_OH=0.1))
#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH), data=data, start=c(a=0.387398, lp_C=-0.156017, lp_N=-0.004407, lp_O=0.124360, r_H=-0.375234, r_C=-0.067189, r_N=-0.153371, r_O=-0.088513, s_HH=0.001, s_CH=0.001, s_NH=0.001, s_OH=0.001))
#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH, s_CC, s_CN, s_CO, s_NN, s_NO, s_OO), data=data, start=c(a=-0.51, lp_C=0.1, lp_N=0.1, lp_O=0.1, r_H=0.1, r_C=0.1, r_N=0.1, r_O=0.1, s_HH=0.1, s_CH=0.1, s_NH=0.1, s_OH=0.1, s_CC=0.1, s_CN=0.1, s_CO=0.1, s_NN=0.1, s_NO=0.1, s_OO=0.1))
summary(model)
data$model <- predict(model)

#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(err(data, model) / data$X.all_e))
sqrt(mean((err(data, model) / data$X.all_e)^2)) 
sd(err(data, model) / data$X.all_e) 
max(abs(err(data, model) / data$X.all_e)) 


par(mfrow = c(3,2))

plot(data$X.all_e, data$CCSDt_CBS1_err)
points(data$X.all_e, data$model, type='p', col='red')
plot(data$X.all_e, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1, data$model, type='p', col='red')
plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1_extdist_q, data$model, type='p', col='red')
plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err - data$model)

###########################################################################
# STEP 4 : lone-pair/radical e/singles/doubles
###########################################################################

print(data$rxn_s_H.H)

print(data$rxn_s_C.H)
print(data$rxn_s_N.H)
print(data$rxn_s_O.H)
#constant plus num e plus num e pair
func <- function(df, a, lp_C=0, lp_N=0, lp_O=0, r_H=0, r_C=0, r_N=0, r_O=0, s_HH=0, s_CH=0, s_NH=0, s_OH=0, s_CC=0, s_CN=0, s_CO=0, s_NN=0, s_NO=0, s_OO=0, d_CC=0, d_CN=0, d_CO=0, d_NN=0, d_NO=0, d_OO=0, t_CC=0, t_CN=0, t_CO=0, t_NN=0, t_NO=0)
{
  return(
    #a 
    #+ lp_C*df$X.lp_C + lp_N*df$X.lp_N + lp_O*df$X.lp_O
    #+ r_H*df$X.r_H + r_C*df$X.r_C + r_N*df$X.r_N + r_O*df$X.r_O 
    #+ s_HH*df$X.s_H.H + s_CH*df$X.s_C.H + s_NH*df$X.s_N.H + s_OH*df$X.s_O.H
    #+ s_CC*df$X.s_C.C + s_CN*df$X.s_C.N + s_CO*df$X.s_C.O 
    #+ s_NN*df$X.s_N.N + s_NO*df$X.s_N.O 
    #+ s_OO*df$X.s_O.O
    #
    a 
    + lp_C*df$rxn_lp_C + lp_N*df$rxn_lp_N + lp_O*df$rxn_lp_O
    + r_H*df$rxn_r_H + r_C*df$rxn_r_C + r_N*df$rxn_r_N + r_O*df$rxn_r_O 
    + s_HH*df$rxn_s_H.H + s_CH*df$rxn_s_C.H + s_NH*df$rxn_s_N.H + s_OH*df$rxn_s_O.H
    + s_CC*df$rxn_s_C.C + s_CN*df$rxn_s_C.N + s_CO*df$rxn_s_C.O 
    + s_NN*df$rxn_s_N.N + s_NO*df$rxn_s_N.O 
    + s_OO*df$rxn_s_O.O
    + d_CC*df$rxn_d_C.C + d_CN*df$rxn_d_C.N + d_CO*df$rxn_d_C.O 
    + d_NN*df$rxn_d_N.N + d_NO*df$rxn_d_N.O 
    + d_OO*df$rxn_d_O.O
  )
}

#Trying all at once
model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, 0, s_CH, s_NH, s_OH, s_CC, s_CN, s_CO, s_NN, s_NO, s_OO, d_CC, d_CN, d_CO, d_NN, d_NO, d_OO), data=data, start=c(a=0.064089, lp_C=-0.024744 , lp_N=-0.045833, lp_O=-0.003477, r_H=-0.183140, r_C=-0.010248, r_N=-0.022725, r_O=-0.068088, s_CH=-0.050479, s_NH=-0.036269, s_OH=0.065607, s_CC=0.047337, s_CN=0.038151, s_CO=0.077671, s_NN=0.042359, s_NO=0.120844, s_OO=0.165714, d_CC=0.01, d_CN=0.01, d_CO=0.01, d_NN=0.01, d_NO=0.01, d_OO=0.01))



#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH), data=data, start=c(a=0.214607, lp_C=-0.073996 , lp_N=-0.024296, lp_O=-0.038120, r_H=-0.347418, r_C=-0.051776, r_N=-0.074486, r_O=-0.087346, s_HH=0.185226, s_CH=-0.075798, s_NH=-0.070630, s_OH=0.1))
#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH), data=data, start=c(a=0.387398, lp_C=-0.156017, lp_N=-0.004407, lp_O=0.124360, r_H=-0.375234, r_C=-0.067189, r_N=-0.153371, r_O=-0.088513, s_HH=0.001, s_CH=0.001, s_NH=0.001, s_OH=0.001))
#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH, s_CC, s_CN, s_CO, s_NN, s_NO, s_OO), data=data, start=c(a=-0.51, lp_C=0.1, lp_N=0.1, lp_O=0.1, r_H=0.1, r_C=0.1, r_N=0.1, r_O=0.1, s_HH=0.1, s_CH=0.1, s_NH=0.1, s_OH=0.1, s_CC=0.1, s_CN=0.1, s_CO=0.1, s_NN=0.1, s_NO=0.1, s_OO=0.1))
summary(model)
data$model <- predict(model)

#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(err(data, model) / data$X.all_e))
sqrt(mean((err(data, model) / data$X.all_e)^2)) 
sd(err(data, model) / data$X.all_e) 
max(abs(err(data, model) / data$X.all_e)) 


par(mfrow = c(3,2))

plot(data$X.all_e, data$CCSDt_CBS1_err)
points(data$X.all_e, data$model, type='p', col='red')
plot(data$X.all_e, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1, data$model, type='p', col='red')
plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1_extdist_q, data$model, type='p', col='red')
plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err - data$model)

###########################################################################
# STEP 4 : lone-pair/radical e/singles/doubles/triples
###########################################################################

print(data$rxn_s_H.H)

print(data$rxn_s_C.H)
print(data$rxn_s_N.H)
print(data$rxn_s_O.H)
#constant plus num e plus num e pair
func <- function(df, a, lp_C=0, lp_N=0, lp_O=0, r_H=0, r_C=0, r_N=0, r_O=0, s_HH=0, s_CH=0, s_NH=0, s_OH=0, s_CC=0, s_CN=0, s_CO=0, s_NN=0, s_NO=0, s_OO=0, d_CC=0, d_CN=0, d_CO=0, d_NN=0, d_NO=0, d_OO=0, t_CC=0, t_CN=0, t_CO=0, t_NN=0, t_NO=0)
{
  return(
    #a 
    #+ lp_C*df$X.lp_C + lp_N*df$X.lp_N + lp_O*df$X.lp_O
    #+ r_H*df$X.r_H + r_C*df$X.r_C + r_N*df$X.r_N + r_O*df$X.r_O 
    #+ s_HH*df$X.s_H.H + s_CH*df$X.s_C.H + s_NH*df$X.s_N.H + s_OH*df$X.s_O.H
    #+ s_CC*df$X.s_C.C + s_CN*df$X.s_C.N + s_CO*df$X.s_C.O 
    #+ s_NN*df$X.s_N.N + s_NO*df$X.s_N.O 
    #+ s_OO*df$X.s_O.O
    #
    a 
    + lp_C*df$rxn_lp_C + lp_N*df$rxn_lp_N + lp_O*df$rxn_lp_O
    + r_H*df$rxn_r_H + r_C*df$rxn_r_C + r_N*df$rxn_r_N + r_O*df$rxn_r_O 
    + s_HH*df$rxn_s_H.H + s_CH*df$rxn_s_C.H + s_NH*df$rxn_s_N.H + s_OH*df$rxn_s_O.H
    + s_CC*df$rxn_s_C.C + s_CN*df$rxn_s_C.N + s_CO*df$rxn_s_C.O 
    + s_NN*df$rxn_s_N.N + s_NO*df$rxn_s_N.O 
    + s_OO*df$rxn_s_O.O
    + d_CC*df$rxn_d_C.C + d_CN*df$rxn_d_C.N + d_CO*df$rxn_d_C.O 
    + d_NN*df$rxn_d_N.N + d_NO*df$rxn_d_N.O 
    + d_OO*df$rxn_d_O.O
    + t_CC*df$rxn_t_C.C + t_CN*df$rxn_t_C.N + t_CO*df$rxn_t_C.O 
    + t_NN*df$rxn_t_N.N + t_NO*df$rxn_t_N.O 
  )
}

#Trying all at once
model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, 0, s_CH, s_NH, s_OH, s_CC, s_CN, s_CO, s_NN, s_NO, s_OO, d_CC, d_CN, d_CO, d_NN, d_NO, d_OO, t_CC, t_CN, t_CO, t_NN, t_NO), data=data, start=c(a=0.064089, lp_C=-0.024744 , lp_N=-0.045833, lp_O=-0.003477, r_H=-0.183140, r_C=-0.010248, r_N=-0.022725, r_O=-0.068088, s_CH=-0.050479, s_NH=-0.036269, s_OH=0.065607, s_CC=0.047337, s_CN=0.038151, s_CO=0.077671, s_NN=0.042359, s_NO=0.120844, s_OO=0.165714, d_CC=0.01, d_CN=0.01, d_CO=0.01, d_NN=0.01, d_NO=0.01, d_OO=0.01, t_CC=0.01, t_CN=0.01, t_CO=0.01, t_NN=0.01, t_NO=0.01))



#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH), data=data, start=c(a=0.214607, lp_C=-0.073996 , lp_N=-0.024296, lp_O=-0.038120, r_H=-0.347418, r_C=-0.051776, r_N=-0.074486, r_O=-0.087346, s_HH=0.185226, s_CH=-0.075798, s_NH=-0.070630, s_OH=0.1))
#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH), data=data, start=c(a=0.387398, lp_C=-0.156017, lp_N=-0.004407, lp_O=0.124360, r_H=-0.375234, r_C=-0.067189, r_N=-0.153371, r_O=-0.088513, s_HH=0.001, s_CH=0.001, s_NH=0.001, s_OH=0.001))
#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, lp_C, lp_N, lp_O, r_H, r_C, r_N, r_O, s_HH, s_CH, s_NH, s_OH, s_CC, s_CN, s_CO, s_NN, s_NO, s_OO), data=data, start=c(a=-0.51, lp_C=0.1, lp_N=0.1, lp_O=0.1, r_H=0.1, r_C=0.1, r_N=0.1, r_O=0.1, s_HH=0.1, s_CH=0.1, s_NH=0.1, s_OH=0.1, s_CC=0.1, s_CN=0.1, s_CO=0.1, s_NN=0.1, s_NO=0.1, s_OO=0.1))
summary(model)
data$model <- predict(model)

#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(err(data, model) / data$X.all_e))
sqrt(mean((err(data, model) / data$X.all_e)^2)) 
sd(err(data, model) / data$X.all_e) 
max(abs(err(data, model) / data$X.all_e)) 


par(mfrow = c(3,2))

plot(data$X.all_e, data$CCSDt_CBS1_err)
points(data$X.all_e, data$model, type='p', col='red')
plot(data$X.all_e, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1, data$model, type='p', col='red')
plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err - data$model)

plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err)
points(data$CCSDt_CBS1_extdist_q, data$model, type='p', col='red')
plot(data$CCSDt_CBS1_extdist_q, data$CCSDt_CBS1_err - data$model)

