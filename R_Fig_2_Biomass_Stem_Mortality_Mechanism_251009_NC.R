##*************Mortality--insect and pathogen damage---diversity**************##
##----------------------------------------------------------------------------##
#Tree diversity reduces tree mortality predominately by decreasing biotic damage#
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
library(sjPlot)
library(ggeffects)
library(MASS)

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


MAGM <- glmmTMB(AGMRsd ~ DeadBioDmg_sc+ logm.SDIs + logCVDBH_sc,
                dispformula = ~logm.SDIs + logCVDBH_sc, 
                family= tweedie(), data= DFA)
summary(MAGM)

sim <- simulateResiduals(fittedModel = MAGM)
plot(sim) 
testZeroInflation(sim)
testDispersion(sim)
check_collinearity(MAGM)
performance::r2_nakagawa(MAGM)
spres <- resid(MAGM)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(MAGM, time = DFA$time)

explained_var <- var(fitted(MAGM))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2

MANM <- glmmTMB(ANMRcsd ~ DeadBioDmg_sc+ logm.SDIs+ logCVDBH_sc,        
                dispformula = ~DeadBioDmg_sc+ logm.SDIs+ logCVDBH_sc,
                family= tweedie(), 
          ##   control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")), ## take forever
                data= DFA)
summary(MANM)
sim <- simulateResiduals(fittedModel = MANM)
plot(sim) 
testZeroInflation(sim)
testDispersion(sim)
check_collinearity(MANM)
spres <- resid(MANM)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(MANM, time = DFA$time)

explained_var <- var(fitted(MAGM))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2


## **** Fig relating mortality and mechanisms *********###

predf4a <- ggpredict(MAGM,  terms = c('DeadBioDmg_sc [all]'))
f4a <- ggplot()+
  geom_point(data = DFA, aes(DeadBioDmg_sc, AGMRsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  geom_line(data= predf4a, aes(x, predicted, color= group))+
  geom_ribbon(data= predf4a, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  ylim(c(0, 6))+
  apatheme+
  labs(y = expression(atop("Biomass mortality", "(% "~yr^-1~")")))+
  labs(x= expression(' '))+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.89'", parse = TRUE, size = 2.8) +
  annotate('text', -Inf, Inf, hjust = -0.7, vjust = 2.8, 
           label="~italic(r) == '0.71'", parse= TRUE, size= 3)+
  annotate("text",-Inf, Inf, hjust = -0.58, vjust = 3.6,
           label="~italic(p) < 0.001",parse=TRUE,size=3)
f4a



predf4b <- ggpredict(MAGM,  terms = c('logm.SDIs [all]'))
f4b <- ggplot()+
  geom_point(data = DFA, aes(logm.SDIs, AGMRsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 6))+
  geom_line(data= predf4b, aes(x, predicted, color= group))+
  geom_ribbon(data= predf4b, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= ' ')+
  labs(x= ' ')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.89'", parse = TRUE, size = 2.8) +
  annotate("text",-Inf, Inf, hjust = -0.7, vjust = 2.8,
           label="~italic(r) == '-0.07'",parse=TRUE,size=3)+
  annotate("text",-Inf, Inf, hjust = -0.7, vjust = 3.6,
           label="~italic(p) == 0.08",parse=TRUE,size=3)
f4b


predf4c <- ggpredict(MAGM,  terms = c('logCVDBH_sc [all]'))
f4c <- ggplot()+
  geom_point(data = DFA, aes(logCVDBH_sc, AGMRsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 6))+
  geom_line(data= predf4c, aes(x, predicted, color= group))+
  geom_ribbon(data= predf4c, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= ' ')+
  labs(x= ' ')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.89'", parse = TRUE, size = 2.8) +
  annotate('text', -Inf, Inf, hjust = -0.7, vjust = 2.8, 
           label="~italic(r) == '0.03'", parse= TRUE, size= 3)+
  annotate("text",-Inf, Inf, hjust = -0.65, vjust = 3.6,
           label="~italic(p) == 0.39",parse=TRUE,size=3)
f4c



predf4d <- ggpredict(MANM,  terms = c('DeadBioDmg_sc [all]'))
f4d <- ggplot()+
  geom_point(data = DFA, aes(DeadBioDmg_sc, ANMRcsd),color= '#b5dcf0', size= 1, pch= 1)+ 
  geom_line(data= predf4d, aes(x, predicted, color= group))+
  geom_ribbon(data= predf4d, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  ylim(c(0, 6))+
  apatheme+
  labs(y = expression(atop("Stem mortality", "(% "~yr^-1~")")))+
  labs(x= expression(Biotic~damage))+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.81'", parse = TRUE, size = 2.8) +
  annotate('text', -Inf, Inf, hjust = -0.7, vjust = 2.8, 
           label="~italic(r) == '0.18'", parse= TRUE, size=3)+
  annotate("text",-Inf, Inf, hjust = -0.55, vjust = 3.6,
           label="~italic(p) < 0.001",parse=TRUE,size=3)
f4d



predf4e <- ggpredict(MANM,  terms = c('logm.SDIs [all]'))
f4e <- ggplot()+
  geom_point(data = DFA, aes(logm.SDIs, ANMRcsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 6))+
  geom_line(data= predf4e, aes(x, predicted, color= group))+
  geom_ribbon(data= predf4e, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= ' ')+
  labs(x= 'Stand density index')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.81'", parse = TRUE, size = 2.8) +
  annotate("text",-Inf, Inf, hjust = -0.7, vjust = 2.8,
           label="~italic(r) == '0.02'",parse=TRUE,size=3)+
  annotate("text",-Inf, Inf, hjust = -0.65, vjust = 3.6,
           label="~italic(p) ==0.47",parse=TRUE,size=3)
f4e



predf4f <- ggpredict(MANM,  terms = c('logCVDBH_sc [all]'))
f4f <- ggplot()+
  geom_point(data = DFA, aes(logCVDBH_sc, ANMRcsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 6))+
  geom_line(data= predf4f, aes(x, predicted, color= group))+
  geom_ribbon(data= predf4f, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= ' ')+
  labs(x= 'Size inequality')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.81'", parse = TRUE, size = 2.8) +
  annotate('text', -Inf, Inf, hjust = -0.65, vjust = 2.8, 
           label="~italic(r) == '-0.12'", parse= TRUE, size= 3)+
  annotate("text",-Inf, Inf, hjust = -0.6, vjust = 3.6,
           label="~italic(p) < 0.001",parse=TRUE,size=3)
f4f


setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/NC/NC_Code")
tiff(filename = '2025_Mechanism_mortality_NC.tiff', width = 2200, height = 1400, res = 300)
cowplot::plot_grid(f4a,  f4c, f4b,f4d,  f4f,f4e, 
                   align = 'hv', labels = 'auto', label_size = 12,
                   ncol = 3)
dev.off()

cor.test(DFA$DeadBioDmgsd_c, DFA$logCVDBH_sc)
cor.test(DFA$DeadBioDmgsd_c, DFA$logm.SDIsd)


## percent change 

100*(mean(DFA$AGMRsd)+0.7072*max(DFA$DeadBioDmgsd_c)-mean(DFA$AGMRsd)+0.7072*min(DFA$DeadBioDmgsd_c))/(mean(DFA$AGMRsd)+0.7072*min(DFA$DeadBioDmgsd_c))
## 1360% for biomass mortality change

100*(mean(DFA$ANMRcsd)+0.1805*max(DFA$DeadBioDmgsd_c)-(mean(DFA$ANMRcsd)+0.1805*min(DFA$DeadBioDmgsd_c)))/(mean(DFA$ANMRcsd)+0.1805*min(DFA$DeadBioDmgsd_c))
## 61% for stem mortality change

