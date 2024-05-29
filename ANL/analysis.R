# Change this for your particular working directory
library ("dplyr")

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
                  ANL_clean_noref$HoF.E0.ccsd.t..tz
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

print(data)


data$CCSDt_CBS1_frac_err <- (data$CCSDt_CBS2 -  data$CCSDt_CBS1)/(data$CCSDt_CBS2)
data$CCSDt_CBS1_frac_expdist <- (data$CCSDt_CBS1 -  data$CCSDt_avqz)/(data$CCSDt_CBS1)


print(data$species[which.max(data$CCSDt_CBS1_frac_err)])


plot(data$CCSDt_CBS1_frac_expdist,  data$CCSDt_CBS1_frac_err)


            



