##*************Mortality--insect and pathogen damage---diversity**************##
##----------------------------------------------------------------------------##
#Tree diversity reduces tree mortality predominately by decreasing biotic damage# 
##                             Xiaxia Ding & Han Chen                         ##
##----------------------------------------------------------------------------##

install.packages("FD")
library(FD);library(ade4); library(ape);library(gridExtra)
library(ggfortify)


setwd("E:/XXD_LU_PhD/Dive_mortality/Mortality_writing/To_Nature_202502/R")
trt_clear <- fread('BC_Data_Traits_250131.csv')
trt_clear_m <- as.matrix(trt_clear[, -c('Genus','Family')], rownames = "Species")

Genus_traits_BC<-autoplot(prcomp(trt_clear_m,scale=TRUE),
                          data = trt,
                          loadings = TRUE, #col="Genus",
                          loadings.colour = "#a8c4ae",
                          loadings.label = T,
                          size=2,loadings.label.size = 2, scale=0,
                          shape= 'Habit')+
  scale_shape_manual(values = c("D" = 17,  # Change to numerical values for shapes
                                "E" = 19), 
                     labels = c("D" = "Broadleaves",
                                "E" = "Conifers")) + 
  geom_text_repel(aes(label= rownames(trt_clear_m)), size= 2, max.overlaps = Inf)+
  guides(size= 'none')+
  theme(panel.border = element_rect(fill = NA, size=0.5),
        panel.background = element_blank())

Genus_traits_BC
tiff("BC_1219_Plots_Genus_Species_Trait_PCA_update2_2025.tiff",width=1900,height=1600,res=300,compression="lzw")
grid.arrange(Genus_traits_BC, nrow=1)
dev.off() 



FD2 <- fread('BC_Data_CWM_Traits_250131.csv')
tiff("BC_1219_Plots_CWM_PCA_update2_2025.tiff",width=1600,height=1600,res=300,compression="lzw")
CWM_traits_BC<-autoplot(prcomp(FD2,scale=TRUE),
                        loadings = TRUE, col="#a8c4ae", loadings.colour = "black",
                        loadings.label = T,
                        size=0.5,loadings.label.size = 3, scale=0,
                        loadings.label.repel=T)+
  xlab(expression(paste('CWMpc1 (36.2%)')))+
  ylab(expression(paste('CWMpc2 (18.7%)')))+
  theme(panel.border = element_rect(fill = NA, size=0.5),
        panel.background = element_blank())
CWM_traits_BC
dev.off() 



















