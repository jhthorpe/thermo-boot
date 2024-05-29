# Change this for your particular working directory
library ("dplyr")

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
                        ZPE = ANL_clean_noref.HoF.E0.ccsd.t..tz)



#hist(data$CCSDt_CBS1_frac_err)
#hist(data$CCSDt_CBS1_frac_expdist)

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
data$CCSDt_CBS1_expdist <- (data$CCSDt_CBS1 -  data$CCSDt_avqz)
data$CCSDt_CBS1_frac_err <- (data$CCSDt_CBS2 -  data$CCSDt_CBS1)/(data$CCSDt_CBS1)
data$CCSDt_CBS1_frac_expdist <- (data$CCSDt_CBS1 -  data$CCSDt_avqz)/(data$CCSDt_CBS1)

#looks like hono is an outlier
print(data$species[which.max(data$CCSDt_CBS1_frac_err)])
data <- data[data$species != "ONOH;t",]

# Sort by x values for pretty plots
data <- data[order(data$CCSDt_CBS1_frac_expdist), ]


#What if we try weighting so that the middle isn't so important?
mu = mean(data$CCSDt_CBS1_frac_expdist)
sig = sd(data$CCSDt_CBS1_frac_expdist)
mval = dnorm(mu, mu, sig)

#Fits
frac_err <- nls(CCSDt_CBS1_frac_err ~ apbx(CCSDt_CBS1_frac_expdist, mya, myb), data=data, start=c(mya=0.1, myb=-1.0))
frac_err_2 <- nls(CCSDt_CBS1_frac_err ~ apbx(CCSDt_CBS1_frac_expdist, mya, myb)*(1 - dnorm(CCSDt_CBS1_frac_expdist, mu, sig)/mval)*predict(frac_err), data=data, start=c(mya=0.1, myb=-1.0))
frac_err_3 <- nls(CCSDt_CBS1_frac_err ~ apbx(CCSDt_CBS1_frac_expdist, mya, myb)*(1 - dnorm(CCSDt_CBS1_frac_expdist, mu, sig)/mval), data=data, start=c(mya=0.1, myb=-1.0))

raw_err_4 <- nls(CCSDt_CBS1_err ~ apbx(CCSDt_CBS1_expdist, mya, myb), data=data, start=c(mya=0.1, myb=-1.0))


print(sqrt(mean(data$CCSDt_CBS1_frac_expdist - predict(frac_err))^2))
print(sqrt(mean(data$CCSDt_CBS1_frac_expdist - predict(frac_err_2))^2))
print(sqrt(mean(data$CCSDt_CBS1_frac_expdist - predict(frac_err_3))^2))

#Plotting errors of the error prediction
par(mfrow = c(3,1))
plot(data$CCSDt_CBS1_frac_expdist,  data$CCSDt_CBS1_frac_err)
lines(data$CCSDt_CBS1_frac_expdist, predict(frac_err), col='red')
lines(data$CCSDt_CBS1_frac_expdist, predict(frac_err_2), col='green')
lines(data$CCSDt_CBS1_frac_expdist, predict(frac_err_3), col='blue')


#plots of errors from our linear / linear-Bayesian predictors
plot(data$CCSDt_CBS1_frac_expdist, data$CCSDt_CBS1_frac_err - predict(frac_err))
points(data$CCSDt_CBS1_frac_expdist, data$CCSDt_CBS1_frac_err - predict(frac_err_2), col='red')
points(data$CCSDt_CBS1_frac_expdist, data$CCSDt_CBS1_frac_err - predict(frac_err_3), col='blue')


#Plots of results
data$CCSDt_CBS1_pred_err <- predict(frac_err) * data$CCSDt_CBS1
data$CCSDt_CBS1_pred_err_2 <- predict(frac_err_2) * data$CCSDt_CBS1
data$CCSDt_CBS1_pred_err_3 <- predict(frac_err_3) * data$CCSDt_CBS1
data$CCSDt_CBS1_pred_err_4 <- predict(raw_err_4) 

#data <- data[order(NTOTAL),]

plot(data$CCSDt_CBS1_err)
lines(data$CCSDt_CBS1_pred_err, col = 'red', type = 'p')
lines(data$CCSDt_CBS1_pred_err_2, col = 'green', type = 'p')
lines(data$CCSDt_CBS1_pred_err_3, col = 'blue', type = 'p')
lines(data$CCSDt_CBS1_pred_err_4, col = 'purple', type = 'p')











###########################################################################
# NEW 
#
# Looks like the vast majority of the error is centered around a mean value. Let's include that info...
#
# So let's model the error of molecule X as:
# 
# err(X) = a + f(X)*CBS(X),
#
# where a is the mean value of the error across the species. Then...
#
# err(X) / CBS(X) = a/CBS(X) + f(X)
# 
mean_err <- mean(data$CCSDt_CBS1_err)

mu = mean(data$CCSDt_CBS1_frac_expdist)
sig = sd(data$CCSDt_CBS1_frac_expdist)
mval = dnorm(mu, mu, sig)

#Adjusted fractional error fit
frac_err_adj <- nls((data$CCSDt_CBS1_err / data$CCSDt_CBS1) - (mean_err/data$CCSDt_CBS1) ~ apbx(CCSDt_CBS1_frac_expdist, mya, myb), data=data, start=c(mya=0.1, myb=-1.0))
frac_err_adj_2 <- nls((data$CCSDt_CBS1_err / data$CCSDt_CBS1) - (mean_err/data$CCSDt_CBS1) ~ apbx(CCSDt_CBS1_frac_expdist, mya, myb)*(1 - dnorm(CCSDt_CBS1_frac_expdist, mu, sig)/mval ), data=data, start=c(mya=0.1, myb=-1.0))
frac_err_adj_3 <- nls((data$CCSDt_CBS1_err / data$CCSDt_CBS1) - (mean_err/data$CCSDt_CBS1) ~ apbx(CCSDt_CBS1_frac_expdist, mya, myb) + apbx(CCSDt_CBS1_frac_expdist, mya, myb)*(1 - dnorm(CCSDt_CBS1_frac_expdist, mu, sig)/mval ), data=data, start=c(mya=0.1, myb=-1.0))
frac_err_adj_4 <- nls((data$CCSDt_CBS1_err / data$CCSDt_CBS1) - (mean_err/data$CCSDt_CBS1) ~ apbx(CCSDt_CBS1_frac_expdist, mya, myb)*(1 - dnorm(CCSDt_CBS1_frac_expdist, mu, sig)/mval )*predict(frac_err_adj), data=data, start=c(mya=0.1, myb=-1.0))


#Plots of the fractional error vs fractional extrapolation distance
plot(data$CCSDt_CBS1_frac_expdist, (data$CCSDt_CBS1_err / data$CCSDt_CBS1) - (mean_err/data$CCSDt_CBS1))
points(data$CCSDt_CBS1_frac_expdist, predict(frac_err_adj), type='l', col='red')
points(data$CCSDt_CBS1_frac_expdist, predict(frac_err_adj_2), type='l', col='blue')
points(data$CCSDt_CBS1_frac_expdist, predict(frac_err_adj_3), type='l', col='green')
points(data$CCSDt_CBS1_frac_expdist, predict(frac_err_adj_4), type='l', col='purple')



#Plot of remaining fractional error after extrapolation distance
#plot(data$CCSDt_CBS1_expdist, data$CCSDt_CBS1_err)
#plot(data$CCSDt_CBS1_expdist, data$CCSDt_CBS1_err)
plot(data$CCSDt_CBS1_expdist, data$CCSDt_CBS1_err - (mean_err + predict(frac_err_adj)*data$CCSDt_CBS1))
err_vs_expdist <- nls(data$CCSDt_CBS1_err - (mean_err + predict(frac_err_adj)*data$CCSDt_CBS1) ~ apbx(CCSDt_CBS1_expdist, mya, myb), data=data, start=c(mya=0.1, myb=0.5))



plot(data$CCSDt_CBS1_err) 
lines(rep(mean_err, length(data$CCSDt_CBS1)), type='l', col='black')
#lines(mean_err + predict(frac_err_adj)*data$CCSDt_CBS1, type='p', col='red')
lines(mean_err + predict(frac_err_adj)*data$CCSDt_CBS1_expdist, type='p', col='red')
lines(mean_err + predict(frac_err_adj_2)*data$CCSDt_CBS1, type='p', col='blue')
lines(mean_err + predict(frac_err_adj_3)*data$CCSDt_CBS1, type='p', col='green')
lines(mean_err + predict(frac_err_adj_4)*data$CCSDt_CBS1, type='p', col='purple')


print(mean_err)
sqrt(mean((data$CCSDt_CBS1_err - mean_err)^2))
sqrt(mean((data$CCSDt_CBS1_err - (mean_err + predict(frac_err_adj)*data$CCSDt_CBS1))^2))
sqrt(mean((data$CCSDt_CBS1_err - (mean_err + predict(frac_err_adj_2)*data$CCSDt_CBS1))^2))
sqrt(mean((data$CCSDt_CBS1_err - (mean_err + predict(frac_err_adj_3)*data$CCSDt_CBS1))^2))
sqrt(mean((data$CCSDt_CBS1_err - (mean_err + predict(frac_err_adj_4)*data$CCSDt_CBS1))^2))

#expdist_adj <- nls(data$CCSDt_CBS1_err ~ mean_err + mya*data$CCSDt_CBS1_expdist, data=data, start=c(mya=-1.5))
#plot(data$CCSDt_CBS1_err) 
#lines(predict(expdist_adj), type='p', col='red')
#lines(mean_err - 0.5*mean_err/data$CCSDt_CBS1_expdist^0.1, type='p', col='red')

###########################################################################
# NEW NEW
#
# This is getting messy, so let's sort all this out
#
# So let's model the error of molecule X as:
# 
# 
#
# 

# MODEL_1
# X is assumed to be the following 
# err(X) = a + f(X)*DST(X) + g(X)*CBS(X)
#
# model 0 : f(X) = 0, g(X) = 0
# model 1 : f(X) = a1, g(X) = a2
# model 2 : f(X) = a1, g(X) = a2 + b2*FRCDST 
#
# Variables:
# df -> dataframe 
fnc_1 <- function(df, a=0, b=0, c=0, d=0, e=0)
{
  a + (b)*df$CCSDt_CBS1_expdist + (c)*df$CCSDt_CBS1 + (d)*df$CCSDt_CBS1_frac_expdist*df$CCSDt_CBS1_frac_expdist
}

fnc_2 <- function(df, a=0, b=0, c=0, d=0, e=0, f=0)
{
  a + (b + c*df$CCSDt_CBS1_expdist + d*df$CCSDt_CBS1)*df$CCSDt_CBS1_expdist + (e + f*df$CCSDt_CBS1)*df$CCSDt_CBS1
}

#adding beysian for the fractional stuff
fnc_3 <- function(df, a=0, b=0, c=0, d=0, e=0, f=0)
{
  a + (b + c*df$CCSDt_CBS1_expdist + d*df$CCSDt_CBS1)*df$CCSDt_CBS1_expdist + (e + f*df$CCSDt_CBS1)*df$CCSDt_CBS1
}


#histograms
#extrapolation distances
mu_dist = mean(data$CCSDt_CBS1_expdist)
sig_dist = sd(data$CCSDt_CBS1_expdist)
mval_dist = dnorm(mu_dist, mu_dist, sig_dist)

#fractional extrapolation distances
mu_frac_dist = mean(data$CCSDt_CBS1_frac_expdist)
sig_frac_dist = sd(data$CCSDt_CBS1_frac_expdist)
mval_frac_dist = dnorm(mu_frac_dist, mu_frac_dist, sig_frac_dist)

#CBS values
mu_CBS = mean(data$CCSDt_CBS1)
sig_CBS = sd(data$CCSDt_CBS1)
mval_CBS = dnorm(mu_CBS, mu_CBS, sig_CBS)

#max error values
max_err = max(abs(data$CCSDt_CBS1_err))
max_frac_err = max(abs(data$CCSDt_CBS1_frac_err))

model_0 <- nls(data$CCSDt_CBS1_err ~ fnc_1(data, mya, 0, 0, 0, 0), data=data, start=c(mya=-0.51))
model_1 <- nls(data$CCSDt_CBS1_err ~ fnc_1(data, mya, myb, myc, 0, 0), data=data, start=c(mya=-0.51, myb=-0.5, myc=0.1))
model_2 <- nls(data$CCSDt_CBS1_err ~ fnc_1(data, mya, myb, myc, myd, 0), data=data, start=c(mya=-0.51, myb=-0.5, myc=0.1, myd = -0.1))
model_3 <- nls(data$CCSDt_CBS1_err ~ fnc_2(data, mya, myb, myc, myd, mye, myf), data=data, start=c(mya=-0.51, myb=0.5, myc=0.1, myd=0.1, mye=0.1, myf=0.1))
model_4 <- nls(data$CCSDt_CBS1_err ~ fnc_2(data, mya, myb, myc, mye), data=data, start=c(mya=-0.51, myb=0.5, myc=0.1, mye=0.1))


#model_2 <- nls(data$CCSDt_CBS1_err ~ fnc_1(data, mya, myb, myc), data=data, start=c(mya=-0.51, myb=-0.5, myc=0.1))
#model_3 <- nls(data$CCSDt_CBS1_err ~ fnc_1(data, mya, myb, myc, myd), data=data, start=c(mya=-0.51, myb=-0.5, myc=0.1, myd=0.))


par(mfrow = c(4,2))


#plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err)
#plot(data$CCSDt_CBS1, data$CCSDt_CBS1_frac_err)
#plot(data$CCSDt_CBS1, data$CCSDt_CBS1_frac_err)

#PREDICITONS HISTOGRAM
plot(data$CCSDt_CBS1_err)

points(predict(model_0), col='black', type='l')
#points(predict(model_1), col='red', type='p')
#points(predict(model_2), col='blue', type='p')
#points(predict(model_3), col='green', type='p')
#points(predict(model_4), col='purple', type='p')


plot(data$CCSDt_CBS1_frac_err )
#points(predict(model_0) / data$CCSDt_CBS1, col='black', type='l')
#points(predict(model_1) / data$CCSDt_CBS1, col='red', type='p')
#points(predict(model_2) / data$CCSDt_CBS1, col='blue', type='p')
#points(predict(model_3) / data$CCSDt_CBS1, col ='green', type='p')
#points(predict(model_4) / data$CCSDt_CBS1, col ='purple', type='p')


#"layer 1" 
plot(data$CCSDt_CBS1_expdist, data$CCSDt_CBS1_err)
plot(data$CCSDt_CBS1_expdist, data$CCSDt_CBS1_frac_err)

#CBS vs err
plot(data$CCSDt_CBS1, data$CCSDt_CBS1_err)
plot(data$CCSDt_CBS1, data$CCSDt_CBS1_frac_err)

#total val vs frac err
plot(data$CCSDt_CBS1_frac_expdist, data$CCSDt_CBS1_err)
plot(data$CCSDt_CBS1_frac_expdist, data$CCSDt_CBS1_frac_err)

mean(abs(data$CCSDt_CBS1_err - predict(model_0)))
mean(abs(data$CCSDt_CBS1_err - predict(model_1)))
mean(abs(data$CCSDt_CBS1_err - predict(model_2)))
mean(abs(data$CCSDt_CBS1_err - predict(model_3)))
mean(abs(data$CCSDt_CBS1_err - predict(model_4)))

sqrt(mean((data$CCSDt_CBS1_err - predict(model_0))^2))
sqrt(mean((data$CCSDt_CBS1_err - predict(model_1))^2))
sqrt(mean((data$CCSDt_CBS1_err - predict(model_2))^2))
sqrt(mean((data$CCSDt_CBS1_err - predict(model_3))^2))
sqrt(mean((data$CCSDt_CBS1_err - predict(model_4))^2))


max(abs(data$CCSDt_CBS1_err - predict(model_0)))
max(abs(data$CCSDt_CBS1_err - predict(model_1)))
max(abs(data$CCSDt_CBS1_err - predict(model_2)))
max(abs(data$CCSDt_CBS1_err - predict(model_3)))
max(abs(data$CCSDt_CBS1_err - predict(model_4)))


mean(abs(data$CCSDt_CBS1_err))
sqrt(mean((data$CCSDt_CBS1_err)^2))
max(abs(data$CCSDt_CBS1_err))

mean(abs(data$CCSDt_CBS1_frac_err))
sqrt(mean((data$CCSDt_CBS1_frac_err)^2))
max(abs(data$CCSDt_CBS1_frac_err))


# JHT TESTING TESTING
#
# STEP 1 : remove the mean

fnc_1 <- function(df, a)
{
  a*(df$CCSDt_CBS1_err/df$CCSDt_CBS1_err)
}

par(mfrow = c(4,2))
model_1 <- nls(data$CCSDt_CBS1_err ~ fnc_1(data, mya), data=data, start=c(mya=-0.51))

err <- function(df, md)
{
  df$CCSDt_CBS1_err - predict(md)
}
frac_err <- function(df, md)
{
  (df$CCSDt_CBS1_err - predict(md))/df$CCSDt_CBS1
}

plot(err(data, model_1))
plot(frac_err(data, model_1))

#"layer 1" 
plot(data$CCSDt_CBS1_expdist, err(data, model_1))
plot(data$CCSDt_CBS1_expdist, frac_err(data, model_1))

#CBS vs err
plot(data$CCSDt_CBS1, err(data, model_1))
plot(data$CCSDt_CBS1, frac_err(data, model_1))

#total val vs frac err
plot(data$CCSDt_CBS1_frac_expdist, err(data, model_1))
plot(data$CCSDt_CBS1_frac_expdist, frac_err(data, model_1))

mean(abs(err(data, model_1)))
sqrt(mean((err(data, model_1))^2))
max(abs(err(data, model_1)))

mean(abs(frac_err(data, model_1)))
sqrt(mean((frac_err(data, model_1))^2))
max(abs(frac_err(data, model_1)))



# JHT TESTING TESTING
#
# STEP 2 : CBS linear model

fnc_2 <- function(df, a=0, b=0)
{
  a + b*df$CCSDt_CBS1
}

par(mfrow = c(4,2))
model_2 <- nls(data$CCSDt_CBS1_err ~ predict(model_1) + fnc_2(data, mya, myb), data=data, start=c(mya=0.1, myb=0.1))
summary(model_2)


plot(err(data, model_2))
plot(frac_err(data, model_2))

#"layer 1" 
plot(data$CCSDt_CBS1_expdist, err(data, model_2))
plot(data$CCSDt_CBS1_expdist, frac_err(data, model_2))

#CBS vs err
plot(data$CCSDt_CBS1, err(data, model_2))
plot(data$CCSDt_CBS1, frac_err(data, model_2))

#total val vs frac err
plot(data$CCSDt_CBS1_frac_expdist, err(data, model_2))
plot(data$CCSDt_CBS1_frac_expdist, frac_err(data, model_2))

mean(abs(err(data, model_2)))
sqrt(mean((err(data, model_2))^2))
max(abs(err(data, model_2)))

mean(abs(frac_err(data, model_2)))
sqrt(mean((frac_err(data, model_2))^2))
max(abs(frac_err(data, model_2)))

print("END OF 2")

# JHT TESTING TESTING
#
# STEP 3 : linear in exp dist
# 

fnc_3 <- function(df, a=0, b=0)
{
  a + b*df$CCSDt_CBS1_expdist
}

par(mfrow = c(4,2))
model_3 <- nls(data$CCSDt_CBS1_err ~ predict(model_2) + fnc_3(data, mya, myb), data=data, start=c(mya=0.1, myb=0.1))
summary(model_3)

plot(err(data, model_3))
plot(frac_err(data, model_3))

#"layer 1" 
plot(data$CCSDt_CBS1_expdist, err(data, model_3))
plot(data$CCSDt_CBS1_expdist, frac_err(data, model_3))

#CBS vs err
plot(data$CCSDt_CBS1, err(data, model_3))
plot(data$CCSDt_CBS1, frac_err(data, model_3))

#total val vs frac err
plot(data$CCSDt_CBS1_frac_expdist, err(data, model_3))
plot(data$CCSDt_CBS1_frac_expdist, frac_err(data, model_3))

mean(abs(err(data, model_3)))
sqrt(mean((err(data, model_3))^2))
max(abs(err(data, model_3)))

mean(abs(frac_err(data, model_3)))
sqrt(mean((frac_err(data, model_3))^2))
max(abs(frac_err(data, model_3)))

print("END OF 3")

# JHT TESTING TESTING
#
# STEP 4 : Linear in coupling btween expdist and CB1
#
# Add quadratic in expdist + frac_expdist
# 

plot(err(data, model_3))
plot(frac_err(data, model_3))

#"layer 1" 
plot(data$CCSDt_CBS1_expdist*data$CCSDt_CBS1_frac_expdist, err(data, model_3))
plot(data$CCSDt_CBS1_expdist*data$CCSDt_CBS1_frac_expdist, frac_err(data, model_3))

#CBS vs err
plot(data$CCSDt_CBS1_expdist*data$CCSDt_CBS1, err(data, model_3))
plot(data$CCSDt_CBS1_expdist*data$CCSDt_CBS1, frac_err(data, model_3))

fnc_4 <- function(df, a=0, b=0)
{
  (a + b*(df$CCSDt_CBS1_frac_expdist*df$CCSDt_CBS1_expdist)^2)
}

par(mfrow = c(4,2))
model_4 <- nls(data$CCSDt_CBS1_err ~ predict(model_3) + fnc_4(data, mya, myb), data=data, start=c(mya=0.1, myb=0.1))
summary(model_4)

#"layer 1" 
#plot(data$CCSDt_CBS1_expdist*data$CCSDt_CBS1_frac_expdist, err(data, model_4))
#plot(data$CCSDt_CBS1_expdist*data$CCSDt_CBS1_frac_expdist, frac_err(data, model_4))

#CBS vs err
#plot(data$CCSDt_CBS1_expdist*data$CCSDt_CBS1, err(data, model_4))
#plot(data$CCSDt_CBS1_expdist*data$CCSDt_CBS1, frac_err(data, model_4))


plot(err(data, model_4))
plot(frac_err(data, model_4))

#"layer 1" 
plot(data$CCSDt_CBS1_expdist, err(data, model_4))
plot(data$CCSDt_CBS1_expdist, frac_err(data, model_4))

#CBS vs err
plot(data$CCSDt_CBS1, err(data, model_4))
plot(data$CCSDt_CBS1, frac_err(data, model_4))

#total val vs frac err
plot(data$CCSDt_CBS1_frac_expdist, err(data, model_4))
plot(data$CCSDt_CBS1_frac_expdist, frac_err(data, model_4))


mean(abs(err(data, model_4)))
sqrt(mean((err(data, model_4))^2))
max(abs(err(data, model_4)))

mean(abs(frac_err(data, model_4)))
sqrt(mean((frac_err(data, model_4))^2))
max(abs(frac_err(data, model_4)))

print("END OF 4")



stop()

#mu_frac_expdist = mean(data$CCSDt_CBS1_frac_expdist)
#sig_frac_expdist = sd(data$CCSDt_CBS1_frac_expdist)
#mval_frac_expdist = dnorm(mu_frac_expdist, mu_frac_expdist, sig_frac_expdist)

#fnc_4 <- function(df, a=0, b=0)
#{
#  (a + b*df$CCSDt_CBS1_frac_expdist)*(1-dnorm(df$CCSDt_CBS1_frac_expdist, mu_frac_expdist, sig_frac_expdist)/mval_frac_expdist)
#}










