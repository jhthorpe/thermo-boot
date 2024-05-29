# Change this for your particular working directory
library(ggplot2)
library(tidyr)
library(tibble)
library(hrbrthemes)
library(dplyr)
library(RColorBrewer)
rm(list = ls())

options(max.print=100000)

setwd('/Users/37348458/thermo-boot/ANL')

#select species for which ATcT has high confidence
ANL <- read.csv('stable_TZ.csv')
ANL_clean <- subset(ANL, ANL$ATcT.UNC < 0.25)

#throw out reference species
# NOTE : THIS REQUIRES THAT ISREF IS ONE OF THE COLUMN NAMES IN YOUR .CSV
tf <- rep(FALSE, times = nrow(ANL_clean))
for (i in 1:nrow(ANL_clean))
{
   if (ANL_clean$ISREF[i] == 1)
  {
    tf[i] <- TRUE
  }
}
ANL_clean$is_ref <- tf
ANL_clean_noref <- ANL_clean[ANL_clean$is_ref == FALSE, ]



#data is the frame used for everything after this point
#Grab a subset for tests
data <- data.frame(ANL_clean_noref$SPECIES,
                  ANL_clean_noref$ATcT.VAL,
                  ANL_clean_noref$ATcT.UNC,
                  ANL_clean_noref$HoF.ruCCSD.T..CBS2,
                  ANL_clean_noref$HoF.ruCCSD.T..CBS1,
                  ANL_clean_noref$HoF.ruCCSD.T..a5z,
                  ANL_clean_noref$HoF.ruCCSD.T..aqz,
                  ANL_clean_noref$HoF.ruCCSD.T..atz,
                  ANL_clean_noref$HoF.E0.ccsd.t..tz,
                  NTOTAL = ANL_clean_noref$X..H + ANL_clean_noref$X..C + ANL_clean_noref$X..N + ANL_clean_noref$X..O
                  )


#throw out nans
data <- na.omit(data)


#Rename for my sanity
data <- data %>% rename(species = ANL_clean_noref.SPECIES,
                        ATcT_val = ANL_clean_noref.ATcT.VAL,
                        ATcT_unc = ANL_clean_noref.ATcT.UNC,
                        CCSDt_CBS2 = ANL_clean_noref.HoF.ruCCSD.T..CBS2,
                        CCSDt_CBS1 = ANL_clean_noref.HoF.ruCCSD.T..CBS1,
                        CCSDt_av5z = ANL_clean_noref.HoF.ruCCSD.T..a5z,
                        CCSDt_avqz = ANL_clean_noref.HoF.ruCCSD.T..aqz,
                        CCSDt_avtz = ANL_clean_noref.HoF.ruCCSD.T..atz,
                        ZPE = ANL_clean_noref.HoF.E0.ccsd.t..tz)



#hist(data$CCSDt_CBS1_frac_err)
#hist(data$CCSDt_CBS1_frac_extdist)

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


###########################################################
#TESTING ZONE
#
# The goal here is to try and get a decent model for 
# the CBS error of {T,Q} vs {Q,5} given the distance from
# Q to {T,Q}.
#
# 

#Generate said data
data$CCSDt_CBS1_err <- (data$CCSDt_CBS2 -  data$CCSDt_CBS1) 
data$CCSDt_CBS1_extdist_q <- (data$CCSDt_CBS1 -  data$CCSDt_avqz)
data$CCSDt_CBS1_extdist_t <- (data$CCSDt_CBS1 -  data$CCSDt_avtz)
data$CCSDt_CBS1_frac_err <- (data$CCSDt_CBS2 -  data$CCSDt_CBS1)/(data$CCSDt_CBS1)
data$CCSDt_CBS1_frac_extdist_q <- (data$CCSDt_CBS1 -  data$CCSDt_avqz)/(data$CCSDt_CBS1)
data$CCSDt_CBS1_frac_extdist_t <- (data$CCSDt_CBS1 -  data$CCSDt_avtz)/(data$CCSDt_CBS1)
data$CCSDt_CBS1_dist_tq <- (data$CCSDt_avqz -  data$CCSDt_avtz)
data$CCSDt_CBS1_frac_dist_tq <- (data$CCSDt_avqz -  data$CCSDt_avtz)/(data$CCSDt_CBS1)

#looks like hono is an outlier
print(data$species[which.max(data$CCSDt_CBS1_frac_err)])
data <- data[data$species != "ONOH;t",]

# Sort by x values for pretty plots
data <- data[order(data$CCSDt_CBS1_frac_extdist_q), ]


###########################################################################
# ERR of ERR functions
err <- function(df, md)
{
  df$CCSDt_CBS1_err - predict(md)
}

frac_err <- function(df, md)
{
  (df$CCSDt_CBS1_err - predict(md))/df$CCSDt_CBS1
}


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

squeez <- function(x, mu, sig)
{
  #dnorm(x, mu, sig)*dnorm(0, 0, 1)/dnorm(mu, mu, sig)
  dnorm(x, mu, sig)/dnorm(mu, mu, sig)
}

stretch <- function(x, mu, sig)
{
  1 - squeez(x, mu, sig)
}

plot(data$CCSDt_CBS1_extdist_q, squeez(data$CCSDt_CBS1_extdist_q, mu_extdist_q, sig_extdist_q))
plot(data$CCSDt_CBS1_frac_extdist_q, stretch(data$CCSDt_CBS1_frac_extdist_q, mu_frac_extdist_q, sig_frac_extdist_q))


###########################################################################
# STEP 1
#
# Simple fit to CBS total
#
# err(X) = a + b*CBS(X)

func <- function(df, a=0, b=0)
{
  a +  (b )*df$CCSDt_CBS1
}

model <- nls(data$CCSDt_CBS1_err ~ func(data, mya, myb), data=data, start=c(mya=-0.51, myb=0.002))
summary(model)


#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(frac_err(data, model)))
sqrt(mean((frac_err(data, model))^2))
max(abs(frac_err(data, model)))


par(mfrow = c(2,2))
plot(data$CCSDt_CBS1_err)
lines(predict(model), col='red', type='p')

plot(data$CCSDt_CBS1_frac_err)
lines(predict(model)/data$CCSDt_CBS1, col='red', type='p')

plot(err(data, model))

plot(frac_err(data, model))




print("END OF 1")




###########################################################################
# STEP 2
#
# Add in linear relationship with extrapolation distance
# and quadratic dependence of extrapolation distance with CBS
#
# err(X) = a + (b + d*QDST(X))*CBS(X) + (c)*QDST(X) 

func <- function(df, a=0, b=0, c=0, d=0)
{
  a +  (b + d*df$CCSDt_CBS1_extdist_q)*df$CCSDt_CBS1 + (c)*df$CCSDt_CBS1_extdist_q  
}

par(mfrow = c(3,2))
plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err)
plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err*2)


model <- nls(data$CCSDt_CBS1_err ~ func(data, mya, myb, myc, myd), data=data, start=c(mya=-0.51, myb=0.002, myc=-0.335, myd=0.003))
summary(model)

#par(mfrow = c(3,2))
#plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err)
#plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err * stretch(data$CCSDt_CBS1_frac_extdist_q, mu_frac_extdist_q, sig_frac_extdist_q))



par(mfrow = c(2,2))
plot(data$CCSDt_CBS1_err)
lines(predict(model), col='red', type='p')

plot(data$CCSDt_CBS1_frac_err)
lines(predict(model)/data$CCSDt_CBS1, col='red', type='p')

plot(err(data, model))

plot(frac_err(data, model))




#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(frac_err(data, model)))
sqrt(mean((frac_err(data, model))^2))
max(abs(frac_err(data, model)))



print("END OF 2")


print("INVESTIGATING DEPENDENCES")
par(mfrow = c(3,2))
plot(data$CCSDt_CBS1_frac_extdist_q, err(data, model))
plot(data$CCSDt_CBS1_frac_extdist_q, frac_err(data, model))
plot(data$CCSDt_CBS1_extdist_q, err(data, model))
plot(data$CCSDt_CBS1_extdist_q, frac_err(data, model))
plot(data$CCSDt_CBS1, err(data, model))
plot(data$CCSDt_CBS1, frac_err(data, model))


###########################################################################
# STEP 3
#
# Quadratic relationship with CBS
#
# err(X) = a + (b + d*QDST(X) + )*CBS(X) + (c )*QDST(X) 

func <- function(df, a=0, b=0, c=0, d=0, e=0)
{
  a +  (b + d*df$CCSDt_CBS1_extdist_q + e*df$CCSDt_CBS1)*df$CCSDt_CBS1 + (c)*df$CCSDt_CBS1_extdist_q 
}

par(mfrow = c(3,2))
plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err)
plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err*2)


model <- nls(data$CCSDt_CBS1_err ~ func(data, mya, myb, myc, myd, mye), data=data, start=c(mya=-0.51, myb=0.002, myc=-0.335, myd=0.003, mye=0.1))
summary(model)

#par(mfrow = c(3,2))
#plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err)
#plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err * stretch(data$CCSDt_CBS1_frac_extdist_q, mu_frac_extdist_q, sig_frac_extdist_q))



par(mfrow = c(2,2))
plot(data$CCSDt_CBS1_err)
lines(predict(model), col='red', type='p')

plot(data$CCSDt_CBS1_frac_err)
lines(predict(model)/data$CCSDt_CBS1, col='red', type='p')

plot(err(data, model))

plot(frac_err(data, model))




#Stats
mean(abs(err(data, model)))
sqrt(mean((err(data, model))^2))
max(abs(err(data, model)))

mean(abs(frac_err(data, model)))
sqrt(mean((frac_err(data, model))^2))
max(abs(frac_err(data, model)))




#Trying my hand at contors
data$model_err <- err(data, model)
data$model_frac_err <- frac_err(data, model)
#data %>% ggplot(aes(x=CCSDt_CBS1_frac_extdist_q, y=CCSDt_CBS1_extdist_q, fill=CCSDt_CBS1_err)) + geom_point(aes(colour=CCSDt_CBS1_err)) + theme_ipsum()
#data %>% ggplot(aes(x=CCSDt_CBS1, y=CCSDt_CBS1_extdist_q, colour=CCSDt_CBS1_err))  + geom_point() + theme_ipsum() + scale_color_viridis_c(option="magma")
#data %>% ggplot(aes(x=CCSDt_CBS1_frac_extdist_q, y=CCSDt_CBS1_extdist_q, colour=model_err))  + geom_point() + theme_ipsum() + scale_color_viridis_c(option="magma")
#data %>% ggplot(aes(x=CCSDt_CBS1_frac_extdist_q, y=CCSDt_CBS1_extdist_q, colour=model_err))  + geom_point()  + scico::scale_color_scico(palette = "vik")




print("END OF 3")

###########################################################################
# STEP 4
#
# modelling the fractional error
#
# err(X) = a + (b + d*QDST(X) + e*CBS(X) + f*CBS(X)^2)*CBS(X) + (c )*QDST(X) 



func_frac <- function(df, a=0, b=0, c=0, d=0)
{
  (a + b*df$CCSDt_CBS1_frac_extdist_q)*stretch(df$CCSDt_CBS1_frac_extdist_q, mu_frac_extdist_q, sig_frac_extdist_q) + c*df$CCSDt_CBS1_frac_extdist_q*df$CCSDt_CBS1*stretch(df$CCSDt_CBS1, mu_CBS1, sig_CBS1)
}

func_err <- function(df, fmd, a=0, b=0, c=0, d=0, e=0)
{
  a + df$CCSDt_CBS1*(b*predict(fmd) + c) 
}


model_frac <- nls(data$CCSDt_CBS1_frac_err ~ func_frac(data, mya, myb, myc), data=data, start=c(mya=0., myb=-1.71, myc=0.01))
model_err <- nls(data$CCSDt_CBS1_err ~ func_err(data, model_frac, mya, myb, myc), data=data, start=c(mya=-0.51, myb=0.33, myc=0.01))
summary(model_frac)
summary(model_err)

par(mfrow = c(3,2))

plot(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err)
points(data$CCSDt_CBS1_frac_extdist_q, data$CCSDt_CBS1_frac_err - frac_err(data, model_err), type='p', col='red')
plot(data$CCSDt_CBS1_frac_extdist_q, frac_err(data, model_err))

plot(data$CCSDt_CBS1_frac_extdist_q/data$CCSDt_CBS1^2, data$CCSDt_CBS1_frac_err)
points(data$CCSDt_CBS1_frac_extdist_q/data$CCSDt_CBS1^2, data$CCSDt_CBS1_frac_err - frac_err(data, model_err), type='p', col='red')
plot(data$CCSDt_CBS1_frac_extdist_q/data$CCSDt_CBS1^2, frac_err(data, model_err))


plot(data$CCSDt_CBS1, data$CCSDt_CBS1_frac_err)
points(data$CCSDt_CBS1, data$CCSDt_CBS1_frac_err - frac_err(data, model_err), type='p', col='red')
plot(data$CCSDt_CBS1, frac_err(data, model_err))





#Stats
mean(abs(err(data, model_err)))
sqrt(mean((err(data, model_err))^2))
max(abs(err(data, model_err)))

mean(abs(frac_err(data, model_err)))
sqrt(mean((frac_err(data, model_err))^2))
max(abs(frac_err(data, model_err)))


print("END OF 4")


