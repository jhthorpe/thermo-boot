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

setwd('/Users/37348458/thermo-boot/TQ_survey')
raw <- read.csv('data_pVDZ.csv')

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
# GCD function for a list of integers
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

####################################################################
# GENERATE PAIR REACTIONS
#
# This proceeds via two steps. The first step generates the reaction names 
# and the number of references (negative indicates LHS and positive indicates RHS)
#
# Once these are generates, the rest of the columns of the dataframe are filled. 

rr <- which(data$Species %in% "n2o")
print(data$Species[rr])
print(data$count_H[rr])
print(data$count_B[rr])
print(data$count_C[rr])
print(data$count_N[rr])
print(data$count_O[rr])
print(data$count_F[rr])

rr <- which(data$Species %in% "hnco")
print(data$Species[rr])
print(data$count_H[rr])
print(data$count_B[rr])
print(data$count_C[rr])
print(data$count_N[rr])
print(data$count_O[rr])
print(data$count_F[rr])

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




print(rxn_names)

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

stop()


####################################################################
# RAW ENERGIES
# We specifically care about...
# (Q)_L - (T)
# (Q)_L - T
# T - (T)
#
# And additionally basis set dependence. 



stop()
# HERE HERE HERE
####################################################################

#Function that generates delta of a given column in the ANL reaction scheme
delta_anl_rxn <- function(noref, ref, col)
{
  H2_rowid <- which(ref$Species %in% "h2")
  CH4_rowid <- which(ref$Species %in% "ch4")
  NH3_rowid <- which(ref$Species %in% "nh3")
  H2O_rowid <- which(ref$Species %in% "h2o")
  HF_rowid <- which(ref$Species %in% "hf")
  B2H6_rowid <- which(ref$Species %in% "b2h6")
  return(
    ref[[col]][H2_rowid] * noref$X.H2
    + ref[[col]][CH4_rowid] * noref$X.C
    + ref[[col]][NH3_rowid] * noref$X.N
    + ref[[col]][H2O_rowid] * noref$X.O
    + ref[[col]][HF_rowid] * noref$X.F
    + ref[[col]][B2H6_rowid] * noref$X.B
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
func <- function(df, a, a_C=0, a_N=0, a_O=0)
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
# STEP 3 : Fit number of atoms + 
###########################################################################

#constant plus num e plus num e pair
func <- function(df, a, a_H=0, a_C=0, a_N=0, a_O=0, a_X=0, a_Y, b_C=0, b_N=0, b_O=0, b_dst=0)
{
  X = (df$X.all_e - (df$X.C*6 + df$X.N*7 + df$X.O*8))
  X = (df$X.val_e - (df$X.C*4 + df$X.N*5 + df$X.O*6))
  return(
    a 
    #+ df$CCSDt_CBS1*((a_H)*df$X.H + (a_C)*df$X.C*4*3/2 + (a_N)*df$X.N*5*4/2 + (a_O)*df$X.O*6*5/2)
    + (a_C)*df$X.C*4*3/2 + (a_N)*df$X.N*5*4/2 + (a_O)*df$X.O*6*5/2 
    + a_X*df$CCSDt_CBS1*((a_H)*df$X.H + (a_C)*df$X.C*4*3/2 + (a_N)*df$X.N*5*4/2 + (a_O)*df$X.O*6*5/2)
    #+ a_X*df$CCSDt_CBS1*((a_H)*df$X.H + (a_C)*df$X.C*4 + (a_N)*df$X.N*5 + (a_O)*df$X.O*6)
  )
  #+ a_X*df$CCSDt_CBS1*((a_C)*df$X.C*4*3/2 + (a_N)*df$X.N*5*4/2 + (a_O)*df$X.O*6*5/2)
  #+ a_X*df$CCSDt_CBS1_extdist_q*((a_C)*df$X.C*4*3/2 + (a_N)*df$X.N*5*4/2 + (a_O)*df$X.O*6*5/2)
  #+ (a_X)*X*(X-1)/2 #not much help
  #+ df$CCSDt_CBS1_extdist_q*((b_C)*df$X.C*4 + (b_N)*df$X.N*5 + (b_O)*df$X.O*6)
  #+ a_X*df$CCSDt_CBS1*((a_C)*df$X.C*4 + (a_N)*df$X.N*5 + (a_O)*df$X.O*6)
  #+ a_X*df$CCSDt_CBS1_extdist_q*((a_H)*df$X.H + (a_C)*df$X.C*4 + (a_N)*df$X.N*5 + (a_O)*df$X.O*6)
}

#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, a_C, a_N, a_O, b_C, b_N, b_O), data=data, start=c(a=-0.1443717, a_C=0.029030, a_N=0.023552, a_O=0.026599, b_C=0.1, b_N=0.1, b_O=0.1))
model <- nls(data$CCSDt_CBS1_err ~ func(data, a, a_H, a_C, a_N, a_O, a_X), data=data, start=c(a=-0.1443717, a_H=0.1, a_C=0.029030, a_N=0.023552, a_O=0.026599, a_X=0.01))
#model <- nls(data$CCSDt_CBS1_err ~ func(data, a, a_C, a_N, a_O), data=data, start=c(a=-0.1443717, a_C=0.029030, a_N=0.023552, a_O=0.026599))

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


print("END OF 3: ATOMS + EXTRA VALANCE")




