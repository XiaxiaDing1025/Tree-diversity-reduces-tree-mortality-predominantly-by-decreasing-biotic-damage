##*************Mortality--insect and pathogen damage---diversity**************##
##----------------------------------------------------------------------------##
##  Phylogenetic diversity reduces tree mortality by decreasing biotic damage ##
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

##********Table S1, Fig. S9, description of variables used before scale*******##

setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/To_Nature_202502/R")
AttriDF <- fread('BC_Data_TableS1_Fig_S1_S9_250131.csv')
library(reshape2)
number_ticks <- function(n) {function(limits) pretty(limits, n)}
AttriDF_hist <- melt(AttriDF, id.vars = 'PLOT_ID')
AttriDF_hist$variable <- factor(AttriDF_hist$variable,
                                levels = c( 'AGMR', 'ANMRc', 'DeadBioDmg','CVDBH',  'm.SDI',
                                            'hillDivq0', 'FDq0', 'hillPDq0','hillDivq1','FDq1', 'hillPDq1','hillDivq2','FDq2', 'hillPDq2',
                                            'PC1', 'PC2', 'MATave', 'CMIave', 'Age', 'Length','Census' , 'Year', 'PS') ,                 
                                labels = c(expression('Biomass mortality (%'*~yr^-1*')'), expression('Stem mortality (%'*~yr^-1*')'),
                                           expression('Biotic damage'~ '(%)'),
                                           expression(Size~inequality), expression("Stand density"~'('*m^3~ ha^-1*')'),
                                           expression(TDq0),expression(FDq0), expression(PDq0),expression(TDq1), expression(FDq1),expression(PDq1),
                                           expression(TDq2), expression(FDq2), expression(PDq2),expression(CWM[pc1]), expression(CWM[pc2]),
                                           expression("MAT"~ "(°C)"), expression("CMI"~"(cm)"), expression("Stand age"~"(years)"),
                                           expression("Census length"~"(years)"), expression(Censuses~number),expression(Average~year), 
                                           expression("Plot size" ~ '('*m^2*')')))

tiff("2025_Histogram of all BC.tiff",height=7500,width=6500,res=600,compression="lzw")
ggplot(AttriDF_hist, aes(x= value))+
  facet_wrap(variable~.,scales="free", switch = "x", ncol=4,labeller = "label_parsed")+
  geom_histogram(bins = 100, color= '#a8c4ae', alpha= 0.2)+
  scale_y_continuous(breaks=number_ticks(4))+
  scale_x_continuous(breaks=number_ticks(4))+
  xlab("")+
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
dev.off()



AttriDF <- AttriDF[, -c('PLOT_ID')]
str(AttriDF)
library(data.table)
AttributeTable <- t(sapply(AttriDF, function(x) list(Mean=mean(x),SD=sd(x),Min=min(x),Max=max(x)))) %>% as.data.frame()
AttributeTable1 <- rownames_to_column(AttributeTable)
AttributeTable1
fwrite(AttributeTable1, 'TableS1_Attributes of mortality.csv')



##************************Spatial plots***************************************##
##coordinates were not available because of data use limitation.
DFA <- fread('NEW_DF_of_Chapter2_Mortality_with_ATA_ACMIA_XXD.csv')
pts<-DFA[,c("PLOT_ID","Lat", "Long", "MATave","CMIave",   
            'm.SDI', 
            'AGMR','ANMRc' ,'hillDivq0' , 'hillPDq0', 'FDq0', 
            'DeadBioDmg', 'CVDBH', 'PC1', 'PC2', 'Age')]
pts<-unique(pts)
factor(pts$PLOT_ID)

library(bcmaps)
bec <- bec()
pts <- st_as_sf(pts, coords = c('Long', 'Lat'))
pts <- st_set_crs(pts, 4326)
pts <- transform_bc_albers(pts)

cor <- st_coordinates(pts) %>% as.data.table()

pts1 <- st_drop_geometry(pts)

pts1 <- cbind(pts1, cor)

##########Plot distribution------------
library(ggspatial)
library(paletteer)

mapCMI <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=CMIave), size= 0.6)+ 
  scale_color_gradientn(colours = c('#d73027', '#fdae61', '#fee08b', '#d9ef8b', '#91bfdb', '#4575b4'),
                        name=c("CMI (cm)")) +
  geom_jitter() +
  theme_classic() +
  theme(legend.key=element_blank()) +
  theme(legend.text = element_text(size = 9), legend.title = element_text(size= 9)) +
  theme(legend.position = "right", legend.box = "vertical", legend.background = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


mapMAT <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=MATave), size= 0.6)+ 
  scale_color_distiller("MAT (°C)", palette = "Spectral", limit=c(-2.14, 10.35))+
  geom_jitter()+theme_classic()+
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())




mapSA <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey') +
  geom_point(data = pts1, aes(x= X, y= Y, col=Age), size= 0.6) + 
  scale_color_gradientn(colours = c('#f7f7f7', '#d4e6f1', '#a8c4d9', '#6fa0bf', '#006f93', '#003f61'),
                        name = "Stand age\n(years)") +
  geom_jitter() +
  theme_classic() +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(size = 9), legend.title = element_text(size = 9)) +
  theme(legend.position = "right", legend.box = "vertical", legend.background = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


mapAGMR <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=AGMR), size= 0.6)+ 
  scale_color_gradientn(colors =c('#8ecae6', '#fcafa9', '#be1e2d','#ef4043', '#ff7e9f' ,'#ff87ab', '#e1b7ff') ,
                        name= expression(atop('Biomass \n mortality', '(%'*~yr^-1*')')))+
  geom_jitter()+theme_classic()+ 
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())




mapANMR <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=ANMRc), size= 0.6)+ 
  scale_color_gradientn(colors =c('#8ecae6', '#fcafa9', '#be1e2d','#ef4043', '#ff7e9f' ,'#ff87ab', '#e1b7ff') ,
                        name= expression(atop('Stem \n mortality', '(%'*~yr^-1*')')))+
  geom_jitter()+theme_classic()+ 
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())



mapDamage <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=DeadBioDmg), size= 0.6)+ 
  scale_color_gradientn(colors =c('#8ecae6',  '#fbdf9d', '#fbc99d', '#fbb39d', '#fba09d','#be1e2d','#ef4043') , 
                        name=expression(atop('Damage', '(%)'))) +
  geom_jitter()+theme_classic()+ 
  #theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())



# PC1: Resistance to Water Stress and Insects
mapPC1 <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=PC1), size= 0.6)+ 
  scale_color_gradientn(colors =c('#FBBEDE', '#ddffff','#d7f1ee','#b6e3e8','#8ecae6', '#219ebc',  '#126782' ) , 
                        name=expression(CWM[pc1]))+
  geom_jitter()+theme_classic()+ 
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())


# PC2: Resource Acquisition
mapPC2 <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=PC2), size= 0.6)+ 
  scale_color_gradientn(colors = c('#fbb040','#fcec52', '#99ca3c', '#208b3a' , '#066839','#0a5c36' , '#287271' ), 
                        name=expression(CWM[pc2]))+
  geom_jitter()+theme_classic()+ 
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())


mapSDI <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey') +
  geom_point(data = pts1, aes(x= X, y= Y, col=m.SDI), size= 0.6) + 
  scale_color_gradientn(colours = c('#C1E7E3', '#DCFFFB', '#FBBEDE', '#E5CCFF', '#EFBBFF', "#D896FF", '#7A5B8D'),
                        breaks = c(146, 922, 1461, 2000, 2841),
                        name = expression(atop("Stand \n density index", '('*m^3~ ha^-1*')'))) +
  geom_jitter() + theme_classic() + 
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(size = 9), legend.title = element_text(size = 9)) +
  theme(legend.position = "right", legend.box = "vertical", legend.background = element_blank()) +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())



mapCVDBH <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=CVDBH),size=0.6)+ 
  geom_point(position = "dodge")+
  scale_color_paletteer_c("ggthemes::Temperature Diverging",
                          name=c("Size \n inequality"))+
  geom_jitter()+theme_classic()+ 
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())




mapTD <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=hillDivq0), size= 0.6)+ 
  scale_color_gradientn(colors = c('#ffe1ac','#ffffa9','#d8faac',  '#208b3a' , '#066839', '#a3f2d8','#287271' ), 
                        name=c("Taxonomic \n diversity") )+
  geom_jitter()+theme_classic()+ 
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())



mapFD <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=FDq0), size= 0.6)+ 
  scale_color_gradientn(colors =c( '#f46303', '#fbb39d', '#9bf0ce', '#39e5ae', '#20b884',  '#208b3a' , '#066839') , 
                        name=c("Functional \n diversity") )+
  geom_jitter()+theme_classic()+ 
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())



mapPD <- ggplot() + 
  geom_sf(data = bec, col= 'lightgrey')+
  geom_point(data = pts1, aes(x= X, y= Y, col=hillPDq0), size= 0.6)+ 
  scale_color_gradientn(colors =c('#be1e2d','#ef4043', '#ff92b1' ,'#ff7e9f' ,'#ff87ab','#e1b7ff', '#dbaaff' ) , 
                        name=c("Phylogenetic \n diversity") )+
  geom_jitter()+theme_classic()+ 
  theme(legend.key=element_blank())+
  theme(legend.text = element_text( size = 9),legend.title = element_text(size= 9))+
  theme(legend.position = "right",legend.box = "vertical", legend.background = element_blank())+
  theme(panel.border = element_rect(colour = "black",fill = NA))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())




rbind_gtable_max <- function(...){
  gtl <- list(...)
  stopifnot(all(sapply(gtl, is.gtable)))
  bind2 <- function (x, y) 
  {
    stopifnot(ncol(x) == ncol(y))
    if (nrow(x) == 0) 
      return(y)
    if (nrow(y) == 0) 
      return(x)
    y$layout$t <- y$layout$t + nrow(x)
    y$layout$b <- y$layout$b + nrow(x)
    x$layout <- rbind(x$layout, y$layout)
    x$heights <- gtable:::insert.unit(x$heights, y$heights)
    x$rownames <- c(x$rownames, y$rownames)
    x$widths <- grid::unit.pmax(x$widths, y$widths)
    x$grobs <- append(x$grobs, y$grobs)
    x
  }
  
  Reduce(bind2, gtl)
}


tiff("2025_Spatial plot BC.tiff",height=7000,width=5500,res=600,compression="lzw")
cowplot::plot_grid(mapAGMR,mapANMR, mapDamage,
                   mapSDI, mapCVDBH, 
                   mapTD, mapFD, mapPD, 
                   mapPC1, mapPC2, 
                   mapMAT, mapCMI,  mapSA,
                   labels = c('auto'),
                   #label_x = 0.2, 
                   ncol=3, align = 'hv')
dev.off()                   
















