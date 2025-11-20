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
library(piecewiseSEM)

rm(list = ls())
oop <- options(na.action = "na.fail")
options(digits=3, scipen=8) 
apatheme=theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/To_Nature_202502/R")

DF<-fread("NEW_DF_of_Chapter2_Mortality_with_ATA_ACMIA_XXD.csv")
DFA <- DF[, c( 'PLOT_ID','AGMR', 'ANMR','loghillDivq0' ,'logFDq0', 'loghillPDq0', 'logPC1', 'logPC2', 
               'DeadBioDmg',  'logm.SDI', 'logCVDBH', 'MATave', 'CMIave', 'logSA', 'y_m_x_m03')]

fwrite(DFA, 'BC_Data_SEM_250131.csv')



DFA <- fread('BC_Data_SEM_250131.csv')


SEM2 <- psem( 
  lm(AGMR ~ DeadBioDmg + logCVDBH + logPC2 + CMIave +logFDq0 ,  data = DFA),
  lm(DeadBioDmg ~ logm.SDI+logCVDBH+loghillPDq0+  logPC1 + logPC2 + MATave + CMIave + logSA,  data = DFA), 
  lmer(logm.SDI ~ logFDq0 + CMIave + logSA + (1|y_m_x_m03), data = DFA),
  lmer(logCVDBH ~ logFDq0 +  logPC1 + logPC2 + MATave + CMIave + logSA+ (1|y_m_x_m03), data = DFA),
  lm(logFDq0 ~ CMIave +MATave, DFA),
  lm(loghillPDq0 ~ MATave+ logSA, DFA),
  
  logCVDBH %~~% logm.SDI,
  logCVDBH %~~% DeadBioDmg,
  logm.SDI %~~% DeadBioDmg,
  logFDq0 %~~% loghillPDq0, 
  logPC1 %~~% logFDq0, 
  logPC1 %~~% loghillPDq0
)

summary(SEM2)  # direct effect of FD0 on BM shall be omitted as this is not based on glmmTMB;
###     other paths shall rely on glmmTMB models if discrepancies occur

COEF_TSSEM <- coefs(SEM2) 
fwrite(COEF_TSSEM, '2025_Biomass mortalty sem.csv')
resid <- residuals(SEM2)
summary(resid)




SEMB <- psem(
  lm(ANMR ~   logPC2 + DeadBioDmg + logm.SDI + logCVDBH, DFA),
  lm(DeadBioDmg ~  logm.SDI+logCVDBH+loghillPDq0+ logPC1 + logPC2 +MATave + CMIave+ logSA , DFA), 
  lmer(logm.SDI ~ logFDq0 + CMIave+ logSA+ (1|y_m_x_m03), DFA),
  lmer(logCVDBH ~ logFDq0 + logPC1 + logPC2 + (1|y_m_x_m03)+MATave + CMIave + logSA , DFA),
  lm(logFDq0 ~ CMIave+MATave, DFA),
  lm(loghillPDq0 ~ MATave+ logSA, DFA),
  
  logCVDBH %~~% logm.SDI,
  logCVDBH %~~% DeadBioDmg,
  logm.SDI %~~% DeadBioDmg,
  logFDq0 %~~% loghillPDq0, 
  logPC1 %~~% logFDq0, 
  logPC1 %~~% loghillPDq0
 )

summary(SEMB)

COEF_TSSEM <- coefs(SEMB) 
fwrite(COEF_TSSEM, '2025_Stem mortalty sem.csv')


##================ direct and indirect path of component sem==================##
setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/NC/NC_Code")
my_theme <-   theme_classic(base_size = 12)  +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, color="black", face="bold.italic"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0, face = "italic"), 
    axis.title.x = element_text(color="black", size=12),
    axis.title.y = element_text(color="black", size=12))

df_BMR <- data.frame(
  path = rep(c("Direct", "Indirect"), times = 10),
  variable = rep(c("BD", "SI", 'SDI','FDq0', "PDq0","CWMpc1",'CWMpc2',"MAT", "CMI", "Standage"), 
                  each = 2),
  value = c(
   0.71, 0, 
   0, -0.11, 
   0, -0.07, 
   0, -0.02, 
   0, -0.11, 
   0, 0.12, 
   -0.07, 0.09, 
   -0.32, -0.23, 
   -0.28, -0.21, 
   -0.23, -0.23
  ))

df_BMR$variable <- factor(df_BMR$variable, 
                          levels = c("BD", "SI", 'SDI','FDq0', "PDq0","CWMpc1",'CWMpc2',"MAT", "CMI", "Standage"),
                          labels = c(expression(Biotic~damage),expression(Size~inequality), 
                                     expression(Stand~density~index), expression(FDq0),expression(PDq0),
                                     expression(CWM[pc1]), expression(CWM[pc2]), 
                                     expression(MAT),expression(CMI),  expression(Stand~age)))

df_BMR$value <- as.numeric(df_BMR$value)
df_BMR$group <- 'Biomass~mortality~rate'
BarBMR <- ggplot(df_BMR,aes(x=variable,y=value,fill=path))+
  geom_bar(stat="identity", width=0.5)+
  facet_grid(~ factor(group, levels= c('Biomass~mortality~rate')), labeller = label_parsed)+
  scale_x_discrete(labels = ggplot2:::parse_safe,name=c(""),limits=rev)+
  ylab(expression(paste(Standardized~effect~size~"("~beta~")")))+
  theme(panel.spacing.y=unit(1, "lines"))+
  geom_hline(yintercept=0)+
  my_theme+ coord_flip()+theme(legend.position =  c(0.8, 0.2),legend.text = element_text(size= 10),
                               legend.title = element_text(size= 10))+theme(plot.margin=unit(c(1,1,1,0),"lines"))

BarBMR
ggsave("2025_Biomass mortlaity_Direct and indirect path effects_NC.png", width = 11, height = 10, units = "cm")




df_SMR <- data.frame(
  path = rep(c("Direct", "Indirect"), times = 10),
  variable = rep(c("BD", "SI", 'SDI','FDq0', "PDq0","CWMpc1",'CWMpc2',"MAT", "CMI", "Standage"), 
                 each = 2),
  value = c(
    0.18, 0, 
    0, -0.15, 
    0, -0.02, 
    0, -0.02, 
    0, -0.03, 
    0, 0.02, 
    0.07, 0.05, 
    -0.05, -0.09, 
    -0.07, -0.04, 
    -0.17, -0.12
  ))

df_SMR$variable <- factor(df_SMR$variable, 
                          levels = c("BD", "SI", 'SDI','FDq0', "PDq0","CWMpc1",'CWMpc2',"MAT", "CMI", "Standage"),
                          labels = c(expression(Biotic~damage),expression(Size~inequality), 
                                     expression(Stand~density~index), expression(FDq0),expression(PDq0),
                                     expression(CWM[pc1]), expression(CWM[pc2]), 
                                     expression(MAT),expression(CMI),  expression(Stand~age)))

df_SMR$value <- as.numeric(df_SMR$value)
df_SMR$group <- 'Stem~mortality~rate'
barSMR <- ggplot(df_SMR,aes(x=variable,y=value,fill=path))+
  geom_bar(stat="identity", width=0.5)+
  facet_grid(~ factor(group, levels= c('Stem~mortality~rate')), labeller = label_parsed)+
  scale_x_discrete(labels = ggplot2:::parse_safe,name=c(""),limits=rev)+
  ylab(expression(paste(Standardized~effect~size~"("~beta~")")))+
  theme(panel.spacing.y=unit(1, "lines"))+
  geom_hline(yintercept=0)+
  my_theme+ coord_flip()+theme(legend.position =  c(0.8, 0.7),legend.text = element_text(size= 10),
                               legend.title = element_text(size= 10))+theme(plot.margin=unit(c(1,1,1,0),"lines"))

barSMR
ggsave("2025_Stem mortlaity_Direct and indirect path effects_NC.png", width = 11, height = 10, units = "cm")




TDSEMBio <- psem( 
  lm(AGMR ~ loghillDivq0 + logPC2 + DeadBioDmg + logCVDBH + CMIave,  data = DFA),
  lm(DeadBioDmg ~logm.SDI+logCVDBH+ loghillDivq0 + logPC1 + logPC2 + MATave + CMIave + logSA,  data = DFA), 
  lmer(logm.SDI ~ loghillDivq0 + CMIave + logSA+ (1|y_m_x_m03), data = DFA),
  lmer(logCVDBH ~ loghillDivq0 + logPC1 + logPC2 + MATave + CMIave + logSA+(1|y_m_x_m03), data = DFA),
  lm(loghillDivq0 ~ MATave+ CMIave +logSA, DFA),
  
  MATave %~~% CMIave,
  logCVDBH %~~% logm.SDI,
  logCVDBH %~~% DeadBioDmg,
  logPC1 %~~% loghillDivq0
)
summary(TDSEMBio)
COEF_TSSEM <- coefs(TDSEMBio) 
fwrite(COEF_TSSEM, '2025_TD_Biomass mortalty sem.csv')


TDSEMStem <- psem( 
  lm(ANMR ~  logPC2 + DeadBioDmg+ logm.SDI + logCVDBH,  data = DFA),
  lm(DeadBioDmg ~ logm.SDI+logCVDBH+loghillDivq0 + logPC1 + logPC2 + MATave + CMIave + logSA,  data = DFA), 
  lmer(logm.SDI ~ loghillDivq0 + CMIave + logSA+ (1|y_m_x_m03), data = DFA),
  lmer(logCVDBH ~ loghillDivq0 + logPC1 + logPC2 + MATave + CMIave + logSA+(1|y_m_x_m03), data = DFA),
  lm(loghillDivq0 ~ MATave+ CMIave+logSA, DFA),
  
  MATave %~~% CMIave,
  logCVDBH %~~% logm.SDI,
  logCVDBH %~~% DeadBioDmg,
  logPC1 %~~% loghillDivq0
)

summary(TDSEMStem)
COEF_TSSEM <- coefs(TDSEMStem) 
fwrite(COEF_TSSEM, '2025_TD_Stem mortalty sem.csv')




TDDmg <- glmmTMB(DeadBioDmgsd ~logm.SDIs+logCVDBH_sc+ loghillDivq0_sc+  logPC1s + logPC2s + 
                   MATave_sc + CMIave_sc+ logSA_sc ,
                 dispformula = ~loghillDivq0_sc+ logPC2s + 
                   MATave_sc + CMIave_sc,
                 family= tweedie, data= DFA)
summary(TDDmg)
sim <- simulateResiduals(fittedModel = TDDmg)
plot(sim) 
testZeroInflation(sim)
testDispersion(sim)
check_collinearity(Dmg)
spres <- resid(Dmg)
bc_moran<-moran.test(spres,listw=nb2listw(nbsr1,zero.policy=TRUE),zero.policy=TRUE)
bc_moran 
DHARMa:: testTemporalAutocorrelation(Dmg, time = DFA$time)

explained_var <- var(fitted(TDDmg))
total_var <- explained_var + var(sim$scaledResiduals)
r2 <- explained_var / total_var
r2

