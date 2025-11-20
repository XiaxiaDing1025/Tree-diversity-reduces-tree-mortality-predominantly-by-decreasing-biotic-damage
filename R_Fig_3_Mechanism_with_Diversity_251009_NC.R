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

##--------------------Test mechanism relationship with diversity--------------##
##*******************************Biotic damage********************************##

Dmg <- glmmTMB(DeadBioDmgsd ~ logm.SDIs+ logCVDBH_sc+logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                 MATave_sc + CMIave_sc+ logSA_sc ,
               dispformula = ~logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
                 MATave_sc + CMIave_sc+ logSA_sc,
               family= tweedie, data= DFA)

summary(Dmg)
sim <- simulateResiduals(fittedModel = Dmg)
plot(sim) 
testZeroInflation(sim)
testDispersion(sim)
check_collinearity(Dmg)
spres <- resid(Dmg)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(Dmg, time = DFA$time)

explained_var <- var(fitted(Dmg))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2



##*****************************Stand density index****************************##
SDI <- lmer(logm.SDIsd ~ logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
              MATave_sc + CMIave_sc+ logSA_sc+ 
              (1|y_m_x_m03),
            # control = lmerControl(optimizer = "nloptwrap",calc.derivs = FALSE), 
            data= DFA)
summary(SDI)
hist(resid(SDI))
check_collinearity(SDI)
spres <- resid(SDI)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(SDI, time = DFA$time)
MuMIn::r.squaredGLMM(SDI)



##**************************Size inequality**********************************##
CV <- lmer(logCVDBHsd ~ logFDq0s + loghillPDq0_sc+ logPC1s + logPC2s + 
             MATave_sc + CMIave_sc+ logSA_sc+(1|y_m_x_m03),DFA)
            
             
summary(CV)
hist(resid(CV))
spres <- resid(CV)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(CV, time = DFA$time)
MuMIn::r.squaredGLMM(CV)



#*********************plot mechanism with Hill diversity***********************#
predf3a <- ggpredict(Dmg,  terms = c('logFDq0s [all]'))
f3a <- ggplot()+
  geom_point(data = DFA, aes(logFDq0s, DeadBioDmgsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0,5))+
  geom_line(data= predf3a, aes(x, predicted, color= group))+
  geom_ribbon(data= predf3a, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y = expression(Biotic~damage~("%")))+
  labs(x= ' ')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.68'", parse = TRUE, size = 2.8) +
  annotate("text", -Inf, Inf, hjust = -0.7, vjust = 2.8,
           label = "~italic(r) == '0.06'", parse = TRUE, size = 3) +
  annotate("text", -Inf, Inf, hjust = -0.65, vjust = 3.6,
           label = "~italic(p) == 0.41", parse = TRUE, size = 3)
f3a



predf3b <- ggpredict(Dmg,  terms = c('loghillPDq0_sc [all]'))
f3b <- ggplot()+
  geom_point(data = DFA, aes(loghillPDq0_sc, DeadBioDmgsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(0, 5))+
  geom_line(data= predf3b, aes(x, predicted, color= group))+
  geom_ribbon(data= predf3b, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= " ")+
  labs(x= ' ')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.68'", parse = TRUE, size = 2.8) +
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 2.8,
           label = "~italic(r) == '-0.15'", parse = TRUE, size = 3) +
  annotate("text", -Inf, Inf, hjust = -0.55, vjust = 3.6,
           label = "~italic(p) == 0.003", parse = TRUE, size = 3)
f3b



predf3c <- ggpredict(SDI,  terms = c('logFDq0s [all]'))
f3c <- ggplot()+
  geom_point(data = DFA, aes(logFDq0s, logm.SDIsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(7.5, 14))+
  geom_line(data= predf3c, aes(x, predicted, color= group))+
  geom_ribbon(data= predf3c, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y = "Stand density index\n(m³ ha⁻¹)")+
  labs(x= 'Hill functional diversity')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.79'", parse = TRUE, size = 2.8) +
  annotate("text", -Inf, Inf, hjust = -0.7, vjust = 2.8,
           label = "~italic(r) == '0.09'", parse = TRUE, size = 3) +
  annotate("text", -Inf, Inf, hjust = -0.55, vjust = 3.6,
           label = "~italic(p) == 0.014", parse = TRUE, size = 3)
f3c



predf3d <- ggpredict(SDI,  terms = c('loghillPDq0_sc [all]'))
f3d <- ggplot()+
  geom_point(data = DFA, aes(loghillPDq0_sc, logm.SDIsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  ylim(c(7.5, 14))+
  geom_line(data= predf3d, aes(x, predicted, color= group))+
  geom_ribbon(data= predf3d, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= ' ')+
  labs(x= 'Hill phylogenetic diversity')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.79'", parse = TRUE, size = 2.8) +
  annotate("text", -Inf, Inf, hjust = -0.7, vjust = 2.8,
           label = "~italic(r) == '0.06'", parse = TRUE, size = 3) +
  annotate("text", -Inf, Inf, hjust = -0.65, vjust = 3.6,
           label = "~italic(p) == 0.13", parse = TRUE, size = 3)

f3d



predf3e <- ggpredict(CV,  terms = c('logFDq0s [all]'))
f3e <- ggplot()+
  geom_point(data = DFA, aes(logFDq0s, logCVDBHsd), color= '#b5dcf0', size= 1, pch= 1)+ 
  geom_line(data= predf3e, aes(x, predicted, color= group))+
  geom_ribbon(data= predf3e, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= 'Size inequality ')+
  labs(x= '')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.79'", parse = TRUE, size = 2.8) +
  annotate("text", -Inf, Inf, hjust = -0.7, vjust = 2.8,
           label = "~italic(r) == '0.12'", parse = TRUE, size = 3) +
  annotate("text", -Inf, Inf, hjust = -0.55, vjust = 3.6,
           label="~italic(p) < 0.001", parse = TRUE, size = 3)

f3e



predf3f <- ggpredict(CV,  terms = c('loghillPDq0_sc [all]'))
f3f <- ggplot()+
  geom_point(data = DFA, aes(loghillPDq1_sc, logCVDBHsd), color= '#b5dcf0', size= 1,pch= 1)+ 
  geom_line(data= predf3f, aes(x, predicted, color= group))+
  geom_ribbon(data= predf3f, aes(x= x, ymin= conf.low, ymax= conf.high, fill= group), alpha= 0.2)+
  labs(y= ' ')+
  labs(x= ' ')+apatheme+
  theme(legend.position = 'none')+
  annotate("text", -Inf, Inf, hjust = -0.6, vjust = 1,
           label = "~italic(R)^2 == '0.79'", parse = TRUE, size = 2.8) +
  annotate("text", -Inf, Inf, hjust = -0.7, vjust = 2.8,
           label = "~italic(r) == '0.04'", parse = TRUE, size = 3) +
  annotate("text", -Inf, Inf, hjust = -0.65, vjust = 3.6,
           label = "~italic(p) == 0.24", parse = TRUE, size = 3)

f3f


setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/NC/NC_Code")
tiff(filename = '2025_Mechanism_Diversity_NC.tiff', width = 1550, height = 2100, res = 300)
cowplot::plot_grid(f3a, f3b,  f3e, f3f,f3c, f3d, label_size = 12,
                   align = 'hv', labels = 'auto',  ncol = 2)
dev.off()



















