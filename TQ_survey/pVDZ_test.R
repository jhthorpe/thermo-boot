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






