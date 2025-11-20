##*************Mortality--insect and pathogen damage---diversity**************##
##----------------------------------------------------------------------------##
#Tree diversity reduces tree mortality predominately by decreasing biotic damage# 
##                             Xiaxia Ding & Han Chen                         ##
##----------------------------------------------------------------------------##

setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/To_Nature_202502/R")
DFA <- fread('BC_Data_Fig_S3_S4_S5_S10_S11_250131.csv')
#------------plot mortality rate against PS and length
# Custom function to plot with dynamic p-value labeling
number_ticks <- function(n) {function(limits) pretty(limits, n)}
plot_custom_ggplot2 <- function(df, x_col, y_col) {
  # Fit the linear model to get p-value
  model <- lm(as.formula(paste(y_col, "~", x_col)), data = df)
  p_value <- summary(model)$coefficients[2, 4]  # Get the p-value for the slope
  # Create a label based on the p-value
  p_label <- if (p_value < 0.001) {
    "P < 0.001"
  } else {
    paste("P =", round(p_value, 3))
  }
  ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_point(size = 1, color = '#a8c4ae', alpha= 0.5, shape= 1) +
    geom_smooth(method = 'lm', se = TRUE, col = 'black', size= 0.5) +
    # Adding the p-value label directly
    annotate("text", x = max(df[[x_col]], na.rm = TRUE) * 0.5, 
             y = max(df[[y_col]], na.rm = TRUE) * 0.9, 
             label = p_label, size = 3, hjust = 0) +
    scale_y_continuous(breaks=number_ticks(4))+
    scale_x_continuous(breaks=number_ticks(4))+
    theme(panel.spacing.y=unit(1, "lines"))+
    theme(
      strip.text.y  = element_blank(),
      strip.text.x =  element_text(size = 10),
      strip.background = element_blank(), 
      #strip.placement = "in",
      plot.background = element_blank(),
      plot.margin=unit(c(1, 1, 1, 1),"lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = NA,fill=NA),
      panel.background = element_blank(),
      plot.title = element_text(size = 10),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      text = element_text(size = 10), 
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.text.x = element_blank(),   # Remove x-axis text
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.line.x = element_blank())
}

#--------observed dependent variables
#--------by plot size

anmr_ps <- plot_custom_ggplot2(DFA,  'PS', 'ANMR')+
  xlab("")+ ylab(expression('Stem mortality (% '*yr^-1*')'))+ylim(c(0, 28))

agmr_ps <- plot_custom_ggplot2(DFA,  'PS','AGMR')+
  xlab(expression('')) + 
  ylab(expression('Biomass mortality (% '*yr^-1*')'))+ylim(c(0, 50))

dmg_ps <- plot_custom_ggplot2(DFA,  'PS', 'DeadBioDmg')+
  xlab(expression('Plot size (m'^2*')')) + 
  ylab(expression('Biotic damage (%)'))+
  scale_x_continuous(breaks = seq(min(DFA$PS), max(DFA$PS), by = 1000)) + 
  theme(
    axis.title.x = element_text(size = 10),  
    axis.line.x  = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),  
    axis.text.x = element_text(size = 10),   
    axis.ticks.length = unit(0.25, "cm") )+ylim(c(0, 100))  



#------by census length
anmr_en <- plot_custom_ggplot2(DFA,  'Length', 'ANMR')+
  xlab("")+ ylab(expression('Stem mortality (% '*yr^-1*')'))+ylim(c(0, 28))

agmr_en <- plot_custom_ggplot2(DFA,  'Length','AGMR')+ylim(c(0, 50))+
  xlab(expression('')) + 
  ylab(expression('Biomass mortality (% '*yr^-1*')'))+ylim(c(0, 50))

dmg_en <- plot_custom_ggplot2(DFA,  'Length', 'DeadBioDmg')+
  xlab("Census length (years)")+ ylab(expression('Biotic damage (%)'))+
  scale_x_continuous(breaks = seq(min(DFA$Length), max(DFA$Length), by = 10)) + 
  theme(
    axis.title.x = element_text(size = 10),  
    axis.line.x  = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),  
    axis.text.x = element_text(size = 10),   
    axis.ticks.length = unit(0.25, "cm") )+ylim(c(0, 100))  

library(cowplot);library(ggpubr);library(gridExtra);library(grid);library(lemon)
tiff("SF_plot of mortality against ps and length.tiff",height=4000,width=3500,res=600,compression="lzw")

plot_grid(
  plot_grid(anmr_ps, agmr_ps, dmg_ps,  ncol = 1, labels = c('a'), label_x = 0, label_y = 1, align = 'v'),
  plot_grid(anmr_en, agmr_en, dmg_en,  ncol = 1, labels = c('b'), label_x = 0, label_y = 1, align = 'v'),
  ncol = 2, align = 'v'
)
dev.off()



#--------corrected dependent variables

#--------by plot size
anmrc_ps <- plot_custom_ggplot2(DFA,  'PS', 'ANMRc')+
  xlab(expression('Plot size (m'^2*')')) +
  ylab(expression('Stem mortality (% '*yr^-1*')'))+
  scale_x_continuous(breaks = seq(min(DFA$PS), max(DFA$PS), by = 1000)) + 
  theme(
    axis.title.x = element_text(size = 10),  
    axis.line.x  = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),  
    axis.text.x = element_text(size = 10),   
    axis.ticks.length = unit(0.25, "cm") ) +ylim(c(0, 28)) 


#------by census length
anmrc_en <- plot_custom_ggplot2(DFA,  'Length', 'ANMRc')+
  xlab("Census length (years)")+ 
  ylab(expression('Stem mortality (% '*yr^-1*')'))+
  scale_x_continuous(breaks = seq(min(DFA$Length), max(DFA$Length), by = 10))+
  theme(
    axis.title.x = element_text(size = 10),  
    axis.line.x  = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),  
    axis.text.x = element_text(size = 10),   
    axis.ticks.length = unit(0.25, "cm") )+ylim(c(0, 28))


tiff("SF_plot of corrected mortality against ps and length.tiff",height=1600,width=3300,res=600,compression="lzw")
plot_grid(
  plot_grid(anmrc_ps, ncol = 1, labels = c('a'), label_x = 0, label_y = 1, align = 'v'),
  plot_grid(anmrc_en, ncol = 1, labels = c('b'), label_x = 0, label_y = 1, align = 'v'),
  ncol = 2, align = 'v'
)
dev.off()


#-----------------------------plot correlation
plot_correlation <- function(x, y) {
  # Calculate correlation
  corr_value <- round(cor(DFA[[x]], DFA[[y]]), 3)
  ggplot(DFA, aes_string(x = x, y = y)) +
    geom_point(size = 1, color = '#a8c4ae', alpha= 0.5, shape= 1) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", size = 0.5) +
    labs( x = NULL, y = NULL) +
    annotate("text", x = Inf, y = Inf, label = paste("r =", corr_value), 
             hjust = 2, vjust = 1.5, color = "black", size = 3) +
    scale_y_continuous(breaks=number_ticks(4))+
    scale_x_continuous(breaks=number_ticks(4))+
    theme(panel.spacing.y=unit(1, "lines"))+
    theme(
      strip.text.y  = element_blank(),
      strip.text.x =  element_text(size = 9),
      strip.background = element_blank(), 
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = NA,fill=NA),
      panel.background = element_blank(),
      plot.title = element_text(size = 9),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 9),
      text = element_text(size = 9), 
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9))
}

plot_ANMR <- plot_correlation("ANMR", "ANMRc")+
  xlab(expression('Observed stem mortality (% '*yr^-1*')'))+
  ylab(expression('Corrected stem mortality (% '*yr^-1*')'))

ggsave("SF_correlation_plot.png", plot = plot_ANMR, width = 3, height = 3, dpi = 300)
dev.off()


##****************************************************************************##
library(ggpmisc)
my.formula <- y ~ x 
number_ticks <- function(n) {function(limits) pretty(limits, n)}
plot_custom_ggplot2 <- function(df, mapping) {
  ggplot(data = df, mapping = mapping) +
    geom_point(pch=1,size=1, color= '#a8c4ae', alpha= 0.5)+
    geom_smooth(method = 'lm',linewidth= 0.5, 
                col ='black',fill= '#ff7f50',  se= TRUE,formula = y~x)+
    stat_poly_eq(aes(label =  paste(after_stat(adj.rr.label), 
                                    after_stat(p.value.label), sep = "*\", \"*")),
                 formula = my.formula,label.x = 0.1, rr.digits =3,p.digits = 3)+
    scale_y_continuous(breaks=number_ticks(4))+
    scale_x_continuous(breaks=number_ticks(4))+
    theme(panel.spacing.y=unit(1, "lines"))+
    theme(
      strip.text.y  = element_blank(),
      strip.text.x =  element_text(size = 12),
      strip.background = element_blank(), 
      strip.placement = "in",
      plot.background = element_blank(),
      plot.margin=unit(c(1,1,1,1),"lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = NA,fill=NA),
      panel.background = element_blank(),
      plot.title = element_text(size = 12),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      text = element_text(size = 12), 
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12))
}


##-------------Bivariate plot of Hill Div, FD, PD at q =0,1,2 with------------##
##--------------------------------Biomass mortality rate----------------------##
library(reshape2)
AttriDF <- DFA[, c( 'AGMRsd', 
                    'loghillDivq0_sc', 'logFDq0s', 'loghillPDq0_sc', 'DeadBioDmg_sc',
                    'loghillDivq1_sc', 'logFDq1s', 'loghillPDq1_sc', 'logCVDBH_sc', 
                    'loghillDivq2_sc', 'logFDq2s', 'loghillPDq2_sc',  'logm.SDIs', 
                    'logPC1s', 'logPC2s',
                    'MATave_sc', 'CMIave_sc', 
                    'logSA_sc')]
#str(AttriDF)

AttriDF_Biv <- melt(AttriDF, id.vars = 'AGMRsd')

AttriDF_Biv$variable <- factor(AttriDF_Biv$variable,
                               levels = 
                                 c( 
                                   'loghillDivq0_sc', 'logFDq0s', 'loghillPDq0_sc', 'DeadBioDmg_sc',
                                   'loghillDivq1_sc', 'logFDq1s', 'loghillPDq1_sc',  'logCVDBH_sc',                               
                                   'loghillDivq2_sc', 'logFDq2s', 'loghillPDq2_sc',  'logm.SDIs',                                                                
                                   'logPC1s', 'logPC2s',
                                   'MATave_sc', 'CMIave_sc',                                
                                   'logSA_sc') ,                                                                             
                               labels =                              
                                 c(
                                   expression(TDq0),expression(FDq0), expression(PDq0),expression(Biotic~damage),
                                   expression(TDq1), expression(FDq1), expression(PDq1),  expression(Size~inequality),                                 
                                   expression(TDq2), expression(FDq2), expression(PDq2),expression(Stand~density~index),
                                   expression('CWM[pc1]'), expression('CWM[pc2]'),
                                   expression(MAT), expression(CMI),                                                            
                                   expression(Stand~age)) )                             

tiff("2025_Bivariate of biomass mortality and variates.tiff",height=7000,width=6500,res=600,compression="lzw")
plot_custom_ggplot2(AttriDF_Biv, aes(x= value, y = AGMRsd))+
  facet_wrap(variable~.,scales="free", switch = "x", ncol=4,labeller = "label_parsed")+
  xlab('')+ylab(expression(Biomass~mortality~rate))
dev.off()




##--------------------------------Stem mortality rate-------------------------##
AttriDF <- DFA[, c( 'ANMRcsd', 
                    'loghillDivq0_sc', 'logFDq0s', 'loghillPDq0_sc', 'DeadBioDmg_sc',
                    'loghillDivq1_sc', 'logFDq1s', 'loghillPDq1_sc', 'logCVDBH_sc', 
                    'loghillDivq2_sc', 'logFDq2s', 'loghillPDq2_sc',  'logm.SDIs', 
                    'logPC1s', 'logPC2s',
                    'MATave_sc', 'CMIave_sc', 
                    'logSA_sc')]

AttriDF_Biv <- melt(AttriDF, id.vars = 'ANMRcsd')

AttriDF_Biv$variable <- factor(AttriDF_Biv$variable,
                               levels = 
                                 c( 
                                   'loghillDivq0_sc', 'logFDq0s', 'loghillPDq0_sc', 'DeadBioDmg_sc',
                                   'loghillDivq1_sc', 'logFDq1s', 'loghillPDq1_sc', 'logCVDBH_sc', 
                                   'loghillDivq2_sc', 'logFDq2s', 'loghillPDq2_sc',  'logm.SDIs', 
                                   'logPC1s', 'logPC2s',
                                   'MATave_sc', 'CMIave_sc', 
                                   'logSA_sc') ,                                                                             
                               labels =                              
                                 c(
                                   expression(TDq0),expression(FDq0), expression(PDq0),expression(Biotic~damage),
                                   expression(TDq1), expression(FDq1), expression(PDq1),  expression(Size~inequality),                                 
                                   expression(TDq2), expression(FDq2), expression(PDq2),expression(Stand~density~index),
                                   expression('CWM[pc1]'), expression('CWM[pc2]'),
                                   expression(MAT), expression(CMI),                                                             
                                   expression(Stand~age)))                              

tiff("2025_Bivariate of stem mortality and variates.tiff",height=7000,width=6500,res=600,compression="lzw")
plot_custom_ggplot2(AttriDF_Biv, aes(x= value, y = ANMRcsd))+
  facet_wrap(variable~.,scales="free", switch = "x", ncol=4,labeller = "label_parsed")+
  xlab('')+ylab(expression(Stem~mortality~rate))
dev.off()

##****************************************************************************##

