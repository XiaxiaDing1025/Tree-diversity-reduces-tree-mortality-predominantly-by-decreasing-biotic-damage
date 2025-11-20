library(see); library(tidyr);library(ggplot2)
library(lme4);library(data.table);library(lmerTest)
library(boot);library(car);library(Matrix);library(TMB)
library(grid);library(gridExtra);library(maps);library(carData)
library(PerformanceAnalytics);library(fBasics);library(car);
library(reshape2);library(vegan);library(dplyr);library(MuMIn);  
library(tidyverse);library(plyr);library(multcomp);library(cAIC4);
library(ggrepel);library(ggfortify);library(mgcv);library(corrplot)
library(broom);library(scales);library(glmmTMB);library(DHARMa)
library(here);library(optimx);library(ggpmisc)

rm(list=ls())
setwd("E:/XXD_LU_PhD/Dive_mortality/Mechanism of NBC/data/Proc_231204")
Ind_BC<-fread("BC0_BA_Bio231205.csv") %>%unique
setnames(Ind_BC, "meas_yr", "FinY")

Ind_BC$Species[Ind_BC$Species=='Salix'] <- 'Salix sp'
Ind_BC$Species[Ind_BC$Species=='Alnus  viridis'] <- 'Alnus viridis'
Ind_BC$Species[Ind_BC$Species == 'Prunus'] <- 'Prunus sp'
Ind_BC$Species <- gsub(" ", "_", Ind_BC$Species)
Ind_BC[, .N, by= Species][order(N),]
Ind_BC[,.N, by= Status] [order(Status),]
D1 <- Ind_BC
D1<-D1[Status=="G"|Status=="I",] 
D1[, .N, by= Species][order(Species),]
# Step 1: Calculate BA per tree
D1[, BA := pi * (DBH / 200)^2]

# Step 2: Summarize to species-plot level
species_BA <- D1[, .(SBA = sum(BA)), by = .(PLOT_ID, FinY, Species)]
species_BA[, TBA := sum(SBA), by = .(PLOT_ID, FinY)]
species_BA[, RBA := (SBA / TBA) * 100]

# Step 3: Calculate Pinus dominance
Pinus <- species_BA[grepl("Pinus", Species), .(Pinus_RBA = sum(RBA)), by = .(PLOT_ID, FinY)]
Pinus <- Pinus[Pinus_RBA >= 70, ]
summary(Pinus)
uniqueN(Pinus$PLOT_ID)
uniqueN(D1$PLOT_ID)



##----------------------------------------------------------------------------##
setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/To_Nature_202502/R")
DF <- fread('NEW_DF_of_Chapter2_Mortality_with_ATA_ACMIA_XXD.csv')
DFA <- DF[!PLOT_ID %in% Pinus$PLOT_ID]
DFA <- DFA[, c( 'FinY', 'AGMRsd', 'ANMRcsd', 'DeadBioDmg_sc', 'logm.SDIs', 'logCVDBH_sc', 
               'logFDq0s', 'loghillPDq0_sc','logPC1s' , 'logPC2s','MATave_sc', 'CMIave_sc', 'logSA_sc')]
fwrite(DFA,'Data of mortality test excluding Pinus dominnat.csv' )
DFA <- fread('Data of mortality test excluding Pinus dominnat.csv')

hist(DFA$AGM_R, breaks = 200)
hist(DFA$logAGMRsd, breaks = 200)
hist(DFA$DeadBioDmg_Bio, breaks = 200)
##********Starting from full model, at q = 0**********************************##
agmq0 <- glmmTMB(AGMRsd ~ logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                   MATave_sc + CMIave_sc + logSA_sc,
                 dispformula = ~ MATave_sc  + poly(logSA_sc,2)  ,
                 family= tweedie(), data= DFA )
summary(agmq0)
library(performance)
check_collinearity(agmq0)
simtwmort <- simulateResiduals(fittedModel = agmq0)
plot(simtwmort)
testZeroInflation(simtwmort)
testDispersion(simtwmort)



MAGM <- glmmTMB(AGMRsd ~ DeadBioDmg_Biocsc+ logm.SDIs + logCVDBH_sc,
                dispformula = ~logm.SDIs + logCVDBH_sc, 
                family= tweedie(), data= DFA)
summary(MAGM)
sim <- simulateResiduals(fittedModel = MAGM)
plot(sim) 
testZeroInflation(sim)
testDispersion(sim)
check_collinearity(MAGM)
performance::r2_nakagawa(MAGM)




Dmg <- glmmTMB(DeadBioDmgsd ~ logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                 MATave_sc + CMIave_sc+ logSA_sc ,
               dispformula = ~logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                 MATave_sc + CMIave_sc+ logSA_sc,
               family= tweedie, data= DFA)

summary(Dmg)




