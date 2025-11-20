##test diversity relationship with biotic damage on live tree and biotic damage on dead tree
library(data.table);library(tidyverse);library(glmmTMB);library(DHARMa)
library(performance);library(spdep)

Livedata <- fread('Biotic damage on live trees test.csv') 

MAGM <- glmmTMB(AGMRsd ~ FinBioDmgsc+ logm.SDIs + logCVDBH_sc,
                #dispformula = ~FinBioDmgsc+logm.SDIs + logCVDBH_sc, 
                family= tweedie(), data= Livedata)
summary(MAGM)
sim <- simulateResiduals(fittedModel = MAGM)
plot(sim) 
testZeroInflation(sim)# p = 0.496
testDispersion(sim)# p = 0.136
check_collinearity(MAGM) # VIF < 2
spres <- resid(MAGM)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(MAGM, time = Livedata$time)# p = 0.98

explained_var <- var(fitted(MAGM))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2 # 0.8199



MANM <- glmmTMB(ANMRcsd ~ FinBioDmgsc+ logm.SDIs + logCVDBH_sc,
                dispformula = ~logm.SDIs + logCVDBH_sc, 
                family= tweedie(), data= Livedata)
summary(MANM)

#AIC       BIC    logLik -2*log(L)  df.resid 
#2375.0    2415.9   -1179.5    2359.0      1211 


#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  0.717625   0.066776  10.747  < 2e-16 ***
#  FinBioDmgsc  0.134857   0.017422   7.741 9.90e-15 ***
#  logm.SDIs    0.002621   0.021210   0.124    0.902    
#logCVDBH_sc -0.114458   0.016628  -6.883 5.85e-12 ***


#Dispersion model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -1.19210    0.15608  -7.638 2.21e-14 ***
#  logm.SDIs   -0.51117    0.04066 -12.571  < 2e-16 ***
#  logCVDBH_sc  0.01553    0.03901   0.398    0.691    
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
sim <- simulateResiduals(fittedModel = MANM)
plot(sim) 
testZeroInflation(sim)# p = 1
testDispersion(sim)# p = 0.232
check_collinearity(MANM) # VIF < 2
#Term  VIF   VIF 95% CI adj. VIF Tolerance Tolerance 95% CI
#FinBioDmgsc 1.17 [1.11, 1.26]     1.08      0.85     [0.79, 0.90]
#logm.SDIs 1.13 [1.08, 1.22]     1.06      0.88     [0.82, 0.93]
#logCVDBH_sc 1.16 [1.10, 1.25]     1.08      0.87     [0.80, 0.91]
spres <- resid(MAGM)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(MANM, time = Livedata$time)# p = 0.9031

explained_var <- var(fitted(MANM))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2 # 0.519066


Dmg <- glmmTMB(FinBioDmgsd ~ logm.SDIs + logCVDBH_sc+logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                 MATave_sc + CMIave_sc+ logSA_sc ,
               dispformula = ~logFDq0s+
                 MATave_sc + CMIave_sc+ logSA_sc,
               family= tweedie, data= Livedata)

summary(Dmg)
sim <- simulateResiduals(fittedModel = Dmg)
plot(sim) 
testZeroInflation(sim)
testDispersion(sim)
check_collinearity(Dmg)
spres <- resid(Dmg)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(Dmg, time = Livedata$time)

explained_var <- var(fitted(Dmg))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2
