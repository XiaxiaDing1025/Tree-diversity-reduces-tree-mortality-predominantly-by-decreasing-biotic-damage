##*************Mortality--insect and pathogen damage---diversity**************##
##----------------------------------------------------------------------------##
#Tree diversity reduces tree mortality predominantly by decreasing biotic damage#
##                             Xiaxia Ding & Han Chen                         ##
##----------------------------------------------------------------------------##

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
##----------------------------------------------------------------------------##


rm(list = ls())
oop <- options(na.action = "na.fail")
options(digits=3, scipen=8) 
apatheme=theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/To_Nature_202502/R")
DFA <- fread('BC_Data_Regression_Model_Fig_1_2_3_250131.csv')
##********Starting from full model, at q = 0**********************************##
agmq0 <- glmmTMB(AGMRsd ~ logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                   MATave_sc + CMIave_sc + logSA_sc,
                 dispformula = ~ loghillPDq0_sc+MATave_sc  + poly(logSA_sc,2)  ,
                 family= tweedie(), data= DFA )

agmq <- lm(AGMRsd ~ logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                   MATave_sc + CMIave_sc + logSA_sc,
           data= DFA )

summary(agmq0)

library(performance)
check_collinearity(agmq0)
simtwmort <- simulateResiduals(fittedModel = agmq0)
plot(simtwmort)
testZeroInflation(simtwmort)
testDispersion(simtwmort)
spres <- resid(agmq0)
#### spatial data can not be made to the public as per our data use agreement with the data providers#####
# bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
# bc_moran 
DHARMa:: testTemporalAutocorrelation(agmq0, time = DFA$time)

explained_var <- var(fitted(agmq0))
total_var <- explained_var + var(simtwmort$scaledResiduals)
r2 <- explained_var / total_var
r2 


##****Starting from full model, at q = 1, shows final mode, biomass mortality****##**
agmq1 <- glmmTMB(AGMRsd ~ logFDq1s + loghillPDq1_sc+ logPC1s + logPC2s + 
                   MATave_sc + CMIave_sc  + logSA_sc,
                 dispformula = ~ loghillPDq1_sc+MATave_sc +poly(logSA_sc,2) ,
                 family= tweedie(), data= DFA )

summary(agmq1)
check_collinearity(agmq1)
simtwmort <- simulateResiduals(fittedModel = agmq1)
plot(simtwmort)
testZeroInflation(simtwmort)
testDispersion(simtwmort)
spres <- resid(agmq1)
# bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
# bc_moran 
DHARMa:: testTemporalAutocorrelation(agmq1, time = DFA$time)#1

# Calculate residual \( R^2 \)
explained_var <- var(fitted(agmq1))
total_var <- explained_var + var(simtwmort$scaledResiduals)
r2 <- explained_var / total_var
r2 


##****Starting from full model, at q = 2, shows final mode, biomass mortality**##** 
agmq2 <- glmmTMB(AGMRsd ~ logFDq2s + loghillPDq2_sc+ logPC1s + logPC2s + 
                   MATave_sc + CMIave_sc + logSA_sc ,
                 dispformula = ~  loghillPDq2_sc+MATave_sc + poly(logSA_sc,2) ,
                 family= tweedie, data= DFA )

summary(agmq2)
simtwmort <- simulateResiduals(fittedModel = agmq2)
plot(simtwmort)
testZeroInflation(simtwmort)
testDispersion(simtwmort)
check_collinearity(agmq2)
spres <- resid(agmq2)
# bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
# bc_moran 
DHARMa:: testTemporalAutocorrelation(agmq2, time = DFA$time) 
explained_var <- var(fitted(agmq2))
total_var <- explained_var + var(simtwmort$scaledResiduals)
r2 <- explained_var / total_var
r2 

AIC(agmq0, agmq1, agmq2)
saveRDS(agmq0, "agmq0.rds") ## avoid time comsuming re-run model
saveRDS(agmq1, "agmq1.rds") 
saveRDS(agmq2, "agmq2.rds") 

agmq0 <- readRDS("agmq0.rds")


##**********Starting from full model, at q = 0, stem mortality****************####*

anmq0 <- glmmTMB(ANMRcsd ~ logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                   MATave_sc+ CMIave_sc + logSA_sc ,
                 dispformula = ~ loghillPDq0_sc+logPC1s+ MATave_sc+CMIave_sc + logSA_sc  ,
                 family= tweedie(), data= DFA )


summary(anmq0)
check_collinearity(anmq0)
simtwmort <- simulateResiduals(fittedModel = anmq0)
plot(simtwmort)
testZeroInflation(simtwmort)
testDispersion(simtwmort)
check_collinearity(anmq1)
spres <- resid(anmq0)
# bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
# bc_moran 
DHARMa:: testTemporalAutocorrelation(anmq0, time = DFA$time)
explained_var <- var(fitted(anmq0))
total_var <- explained_var + var(simtwmort$scaledResiduals)
r2 <- explained_var / total_var
r2

##**********Starting from full model, at q = 1, stem mortality****************####*

anmq1 <- glmmTMB(ANMRcsd ~ logFDq1s + loghillPDq1_sc+ logPC1s + logPC2s + 
                   MATave_sc+ CMIave_sc + logSA_sc ,
                 dispformula = ~ loghillPDq1_sc+logPC1s+ MATave_sc+CMIave_sc + logSA_sc  ,
                 family= tweedie(), data= DFA )


summary(anmq1)
check_collinearity(anmq1)
simtwmort <- simulateResiduals(fittedModel = anmq1)
plot(simtwmort)
testZeroInflation(simtwmort)
testDispersion(simtwmort)
check_collinearity(anmq1)
spres <- resid(anmq1)
# bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
# bc_moran 
DHARMa:: testTemporalAutocorrelation(anmq1, time = DFA$time)
explained_var <- var(fitted(anmq1))
total_var <- explained_var + var(simtwmort$scaledResiduals)
r2 <- explained_var / total_var
r2


##**********Starting from full model, at q = 2, stem mortality****************####*

anmq2 <- glmmTMB(ANMRcsd ~ logFDq2s + loghillPDq2_sc+ logPC1s + logPC2s + 
                   MATave_sc+ CMIave_sc + logSA_sc ,
                 dispformula = ~ loghillPDq2_sc+logPC1s+ MATave_sc+CMIave_sc + logSA_sc  ,
                 family= tweedie(), data= DFA )


summary(anmq2)
check_collinearity(anmq2)
simtwmort <- simulateResiduals(fittedModel = anmq1)
plot(simtwmort)
testZeroInflation(simtwmort)
testDispersion(simtwmort)
check_collinearity(anmq1)
spres <- resid(anmq2)
# bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
# bc_moran 
DHARMa:: testTemporalAutocorrelation(anmq2, time = DFA$time)
explained_var <- var(fitted(anmq1))
total_var <- explained_var + var(simtwmort$scaledResiduals)
r2 <- explained_var / total_var
r2

AIC(anmq0, anmq1, anmq2)
saveRDS(anmq0, "anmq0.rds") ## avoid re-run model 
saveRDS(anmq1, "agmq1.rds") 
saveRDS(anmq2, "agmq2.rds") 

##********Fig 1 ****************************************####
library(sjPlot)
library(ggeffects)
library(base)
DFA <- fread('BC_Data_Regression_Model_Fig_1_2_3_250131.csv')
agmq0 <- readRDS("agmq0.rds")

## Biomass panel
## 6 rows of values were excluded from the plot (only each point), not the regression
## Please refer to the bivariate plot
predmortfd <- ggpredict(agmq0,  terms = c('logFDq0s [all]'))
ppmortfd <- ggplot()+
  geom_point(data = DFA, aes(logFDq0s, AGMRsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 6))+
  geom_line(data= predmortfd, aes(x, predicted, color= group))+
  geom_ribbon(data= predmortfd, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y = expression(Biomass~mortality~("%"~yr^-1)))+
  labs(x= ' ')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.48'", parse = TRUE, size = 2.8) +
  annotate("text",-Inf, Inf, hjust = -0.7, vjust = 2.8,
          label="~italic(r) == '0.07'",parse=TRUE,size=3)+
  annotate("text",-Inf, Inf, hjust = -0.65, vjust = 3.6,
           label="~italic(p) == 0.15",parse=TRUE,size=3)
   
ppmortfd

#### PD0 ###################
predmortpd <- ggpredict(agmq0,  terms = c('loghillPDq0_sc [all]'))
ppmortpd <- ggplot()+
  geom_point(data = DFA, aes(loghillPDq0_sc, AGMRsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 6))+
  geom_line(data= predmortpd, aes(x, predicted, color= group))+
  geom_ribbon(data= predmortpd, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= ' ')+
  labs(x= ' ')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.48'", parse = TRUE, size = 2.8) +
  annotate('text', -Inf, Inf, hjust = -0.7, vjust = 2.8, 
          label="~italic(r) == '-0.11'", parse= TRUE, size= 3)+
  annotate("text",-Inf, Inf, hjust = -0.65, vjust = 3.6,
           label="~italic(p) == 0.03",parse=TRUE,size=3)
ppmortpd

# DFA$pred_BM_PD<-resid(agmq0)+ mean(DFA$AGMRsd) + (-0.11408) *DFA$loghillPDq0_sc
# ggplot(data=DFA, mapping=aes(x=loghillPDq0_sc, y=pred_BM_PD))+geom_point()+geom_smooth()

## Stem panel   
## 4 rows of values were excluded (here we say the points), not influence the regression
predanmfd <- ggpredict(anmq0,  terms = c('logFDq0s [all]'))
ppanmfd <- ggplot()+
  geom_point(data = DFA, aes(logFDq0s, ANMRcsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 6))+
  geom_line(data= predanmfd, aes(x, predicted, color= group))+
  geom_ribbon(data= predanmfd, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y = expression(Stem~mortality~("%"~yr^-1)))+
  labs(x= 'Hill functional diversity')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.50'", parse = TRUE, size = 2.8) +
  annotate("text",-Inf, Inf, hjust = -0.7, vjust = 2.8,
           label="~italic(r) == '0.01'",parse=TRUE,size=3)+
  annotate("text",-Inf, Inf, hjust = -0.65, vjust = 3.6,
           label="~italic(p) == 0.64",parse=TRUE,size=3)
ppanmfd


predanmpd <- ggpredict(anmq0,  terms = c('loghillPDq0_sc [all]'))
ppanmpd <- ggplot()+
  geom_point(data = DFA, aes(loghillPDq0_sc, ANMRcsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 6))+
  geom_line(data= predanmpd, aes(x, predicted, color= group))+
  geom_ribbon(data= predanmpd, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= ' ')+
  labs(x= 'Hill phylogenetic diversity')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.50'", parse = TRUE, size = 2.8) +
  annotate('text', -Inf, Inf, hjust = -0.7, vjust = 2.8, 
           label="~italic(r) == '-0.07'", parse= TRUE, size= 3)+
  annotate("text",-Inf, Inf, hjust = -0.65, vjust = 3.6,
           label="~italic(p) == 0.02",parse=TRUE,size=3)
ppanmpd


setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/NC/NC_Code")
tiff(filename = '2025_Biomass_Stem_Mortality_Diversity_NC.tiff', width = 1800, height = 1600, res = 300)
cowplot::plot_grid(ppmortfd, ppmortpd, ppanmfd, ppanmpd,
                   align = 'hv', labels = 'auto', label_size = 12,
                   ncol = 2)
dev.off()

### percent change for total PD effect on mortality 
## Total effect on biomass mortality = 0.17*0.71 
## Total effect on stem mortality = 0.17*0.18 

agmq0 <- readRDS("agmq0.rds")
anmq0 <- readRDS("anmq0.rds")

BM_minPD<-mean(DFA$AGMRsd)-(0.15*.71)*min(DFA$loghillPDq0_sc)
BM_maxPD<-mean(DFA$AGMRsd)-(0.15*.71)*max(DFA$loghillPDq0_sc)
## if mean is used like our previous PNAS paper
percent_BM1<-(BM_maxPD-BM_minPD)/((BM_minPD+BM_maxPD)/2)*100 ##reduced by 65%


SM_minPD<-mean(DFA$ANMRcsd)-(0.15*.18)*min(DFA$loghillPDq0_sc)
SM_maxPD<-mean(DFA$ANMRcsd)-(0.15*.18)*max(DFA$loghillPDq0_sc)
percent_SM1<-(SM_maxPD-SM_minPD)/((SM_minPD+SM_maxPD)/2)*100 ##reduced by 8%


## assuming total effects = (0.12*-0.16*0.71)+(0.09*-0.1*0.71)
BM_minFD<-mean(DFA$AGMRsd)-0.02*min(DFA$logFDq0s)
BM_maxFD<-mean(DFA$AGMRsd)-0.02*max(DFA$logFDq0s)
## if mean is used like our previous PNAS paper
percent_BM1<-(BM_maxFD-BM_minFD)/((BM_minFD+BM_maxFD)/2)*100 ##reduced by 65%


SM_minFD<-mean(DFA$ANMRcsd)-0.005*min(DFA$logFDq0s)
SM_maxFD<-mean(DFA$ANMRcsd)-0.005*max(DFA$logFDq0s)
percent_SM1<-(SM_maxPD-SM_minPD)/((SM_minPD+SM_maxPD)/2)*100 ##reduced by 8%

##****************Mortality relationship with taxonomic diversity*************## 
TD <- glmmTMB(AGMRsd~ loghillDivq0_sc+ logPC1s + logPC2s + 
                MATave_sc + CMIave_sc + logSA_sc ,
              dispformula = ~ logPC1s +  
                MATave_sc + poly(logSA_sc,2) ,
              family= tweedie, data= DFA)
summary(TD)
sim <- simulateResiduals(fittedModel = TD)
plot(sim)
testZeroInflation(sim)
testDispersion(sim)
res <- resid(TD)
bc_moran<-moran.test(res,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 

explained_var <- var(fitted(TD))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2



ANMTD <- glmmTMB(ANMRcsd~ loghillDivq0_sc+ logPC1s + logPC2s + 
                   MATave_sc + CMIave_sc  + logSA_sc ,
                 dispformula = ~ loghillDivq0_sc+  logPC1s + logPC2s + 
                   MATave_sc + CMIave_sc  + poly(logSA_sc,2) ,
                 family= tweedie, data= DFA)
summary(ANMTD)
sim <- simulateResiduals(fittedModel = ANMTD)
plot(sim)
testZeroInflation(sim)#1
testDispersion(sim)#0.8
res <- resid(ANMTD)
bc_moran<-moran.test(res,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 

explained_var <- var(fitted(ANMTD))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2 #0.5








