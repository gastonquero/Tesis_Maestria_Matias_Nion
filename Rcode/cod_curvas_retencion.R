###########################################################################
# codigo para el analisis de Datos curvas de retencion de humedad
# ensayo de Tesis Maestria
# Ensayo MN.3
# Analisis de datos Gaston Quero - Matias Nion
# 29/07/2021
##################################################################################


getwd ()
setwd ("C:/Users/Usuario/OneDrive/Documentos/Tesis_Maestria_Matias")


# Paquetes 
library (lme4)
library (emmeans)
library ("car")
library ("nlmrt")
library ("easynls")
library (tidyverse)
library ("ggplot2")       
library ("lattice")
library ("latticeExtra")
library (multcompView)
library (multcomp)
library ("dplyr")
library (ggjoy)
library ("ggridges")
library (hrbrthemes)
library(tidyverse)
library(forcats)
library("viridis")
library("lmerTest")
library(lubridate)
library(nycflights13)
library(nlme)
library(xtable)
library(stringr)
library(data.table)
library(svMisc)
library(ggpubr)
library("ggsci")
library("FactoMineR")
library("factoextra")
library("corrplot")

## cargar datos
CC.MN.1 <- read_delim (file = "./Data/rawdata/macetas_cc.txt" , 
                                                   delim ="\t", quote = "\"", 
                                                   escape_backslash = FALSE,
                                                   escape_double = TRUE, 
                                                   col_names = TRUE, 
                                                   col_types = NULL,
                                                   locale = default_locale(), 
                                                   na = "NA")

today()
CC.MN.1$maceta <- as.factor (CC.MN.1$maceta)

CC.MN.1 <- CC.MN.1 %>%
           dplyr::mutate (date = dmy (CC.MN.1$fecha)) %>%
           dplyr::mutate (agua = peso.sust.mac - peso.sust.seco - peso.bolsa) %>%
           dplyr::mutate (hp.porc = (agua *100)/peso.sust.seco)
         


agua.x <- 17.300 - 16.20 - 0.059 

hp.porc.x <- ( 1.041 *100)/16.20

agua.WC2 <- (5*16.20)/100

SPM.2 <- 0.81 + 16.20 + 0.059 

ggscatter (CC.MN.1 , x = "dias", y = "hp.porc", facet.by = "maceta",
           #ylim = c(20,22),
           color = "maceta", shape = 16, size = 3, # Points color, shape and size
           add = "loess"  # Add regressin line
           #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
           #conf.int = TRUE, # Add confidence interval
           #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)


ggscatter (CC.MN.1 , x = "dias", y = "peso.sust.mac",
           #ylim = c(20,22),
           color = "maceta", shape = 16, size = 3, # Points color, shape and size
           add = "loess",  # Add regressin line
           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
           conf.int = TRUE, # Add confidence interval
           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)


ggscatter (CC.MN.1 , x = "dias", y = "peso.sust.mac", facet.by = "maceta",
           color = "maceta", shape = 16, size = 3, # Points color, shape and size
           add = "loess",  # Add regressin line
           add.params = list(color = "blue", fill = "lightgray")
           #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)

summary (CC.MN.1$hp.porc)
ggscatter (CC.MN.1 , x = "dias", y = "hp.porc", facet.by = "maceta",
           color = "maceta", shape = 16, size = 3, # Points color, shape and size
           add = "loess",  # Add regressin line
           add.params = list(color = "blue", fill = "lightgray")
           #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)

CC.MN.1$maceta <- as.factor (CC.MN.1$maceta)

plot ( hp.porc ~ dias, pch=16, 
       ylim =c(15,32),
       xlim =c(0,35),
       type ="n",
       data= CC.MN.1)

pot.27 <- CC.MN.1 %>%
          dplyr::filter (maceta== "27")

pot.29 <- CC.MN.1 %>%
  dplyr::filter (maceta== "29")

pot.25 <- CC.MN.1 %>%
  dplyr::filter (maceta== "25")

pot.6 <- CC.MN.1 %>%
  dplyr::filter (maceta== "6")

pot.26 <- CC.MN.1 %>%
  dplyr::filter (maceta== "26")

pot.30 <- CC.MN.1 %>%
  dplyr::filter (maceta== "30")


points (hp.porc ~ dias, pch=16, 
         col="red",type ="b", 
         data=pot.27)

points (hp.porc ~ dias, pch=17, 
        col="blue",type ="b", 
        data=pot.29)

points (hp.porc ~ dias, pch=18, 
        col="black",type ="b", 
        data=pot.25)

points (hp.porc ~ dias, pch=19, 
        col="orange",type ="b", 
        data=pot.6)

points (hp.porc ~ dias, pch=20, 
        col="green",type ="b", 
        data=pot.26)


points (hp.porc ~ dias, pch=21, 
        col="darkgreen",type ="b", 
        data=pot.30)


4.641 - 4.096


xx.0 <- CC.MN.1 %>%
        filter (dias==0)


xx.14 <- CC.MN.1 %>%
         filter (dias==14)

xx.16 <- CC.MN.1 %>%
       filter (dias==16)

xx.34 <- CC.MN.1 %>%
         filter (dias==34)

summary (xx.0)
summary (xx.14)

summary (xx.16)

summary (xx.34)


3.989 - 3.908
4.641 - 4.015
4.533 - 3.908

3.908 - 3.568

## Calculo para definir tratamiento control 

agua.x34 <- 19.82 - 16.20 - 0.059 

hp.porc.x34 <- (agua.x34 *100)/16.20

agua.x.control <- 19 - 16.20 - 0.059 

hp.porc.x.ctrl <- (agua.x.control *100)/16.20

## Calculo para definir tratamiento estres 

agua.x.deficit <- 17.5 - 16.20 - 0.059 


hp.porc.x.dfct <- (agua.x.deficit  *100)/16.20

## delta agua 

agua.x.control - agua.x.deficit

hp.porc.x.ctrl - hp.porc.x.dfct

## cargar datos
CC.MN.1.WP4C <- read_delim (file = "./Data/rawdata/curva_de_retencion.txt" , 
                       delim ="\t", quote = "\"", 
                       escape_backslash = FALSE,
                       escape_double = TRUE, 
                       col_names = TRUE, 
                       col_types = NULL,
                       locale = default_locale(), 
                       na = "NA")

CC.MN.1.WP4C$curva <- as.factor (CC.MN.1.WP4C$curva)

CC.MN.1.WP4C  <- CC.MN.1.WP4C %>%
             dplyr::mutate (agua =  placaf - pesoplaca - tierras) %>%
             dplyr::mutate (hp.porc = (agua *100)/tierras) %>%
             dplyr::mutate (pot.1 = pot*-1)

summary(CC.MN.1.WP4C)




ggscatter (CC.MN.1.WP4C , x = "hp.porc", y = "pot.1", 
           ylim = c(0,11),
           facet.by = "curva",
           color = "curva", shape = 16, size = 3, # Points color, shape and size
           add = "loess",  # Add regressin line
           add.params = list(color = "blue", fill = "lightgray")
           #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)
plot ( pot ~ hp.porc, pch=16, 
      ylim =c(0,-2),
      xlim =c(0,30),
      type ="n",
      data= CC.MN.1.WP4C)

unique(CC.MN.1.WP4C$curva)

c1 <- CC.MN.1.WP4C %>%
      dplyr::filter (curva == "c1")
c2 <- CC.MN.1.WP4C %>%
  dplyr::filter (curva == "c2")

c3 <- CC.MN.1.WP4C %>%
  dplyr::filter (curva == "c3")

c4 <- CC.MN.1.WP4C %>%
  dplyr::filter (curva == "c4")

points (pot ~ hp.porc, pch=16, col="red", type="p",
        data=c1)

points (pot ~ hp.porc, pch=16, col="blue",
        data=c2)

points (pot ~ hp.porc, pch=16, col="darkgreen",
        data=c3)

points (pot ~ hp.porc, pch=16, col="black",
        data=c4)


abline (h=-1.5, lty=2)
abline (h= 0, lty=2)
#abline (h=-1.5/2, lty=2)

#abline (v= hp.porc.x34, col="red", lwd=2.5)
abline (v= hp.porc.x.ctrl, col="blue", lwd=3)
abline (h= -0.15, col="blue", lty=3, lwd=3)

abline (v=   hp.porc.x.dfct , col="darkgreen", lwd=3)
abline (h= -0.35, col="darkgreen", lty=2, lwd=3 )

-0.15 + 0.35 

abline (v=  21 , col="red", lwd=1.5)
abline (h= 0, col="red", lwd=1.5)


abline (v=6.444444, col="blue", lwd=3 )
abline (h= -0.4, col="blue", lwd=3 )
abline (v= 2.5, col="red", lwd=3 )
abline (h=-0.9)

abline (v= 7.641975, col="blue", lwd=1.5 )


abline (v= 10, col="black", lwd=1.5, lty= 2 )
abline (h= -1, col="black", lwd=1.5, lty= 2 )
abline (v= 2.5, col="red", lwd=3 )
abline (h=-0.9)


2.5*5/ 25

