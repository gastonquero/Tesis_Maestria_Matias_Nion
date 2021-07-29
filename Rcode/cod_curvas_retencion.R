###########################################################################
# codigo para el analisis de Datos curvas de retencion de humedad
# ensayo de Tesis Maestria
# Ensayo MN.3
# Analisis de datos Gaston Quero - Matias Nion
# 29/07/2021
##################################################################################


getwd ()
setwd ("R:/Tesis_Maestria_Matias_Nion")


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

# cargar datos
# cargo los datos de la curva de retencion con el WP4C


CC.MN.3.WP4C <- read_delim (file = "./Data/rawdata/curva_de_retencion_WP4C.txt" , 
                            delim ="\t", quote = "\"", 
                            escape_backslash = FALSE,
                            escape_double = TRUE, 
                            col_names = TRUE, 
                            col_types = NULL,
                            locale = default_locale(), 
                            na = "NA")

CC.MN.3.WP4C$curva <- as.factor (CC.MN.3.WP4C$curva)

CC.MN.3.WP4C  <- CC.MN.3.WP4C %>%
                 dplyr::mutate (agua =  placaf - pesoplaca - tierras) %>%
                 dplyr::mutate (hp.porc = (agua *100)/tierras) %>%
                 dplyr::mutate (pot.1 = pot*-1)

summary(CC.MN.3.WP4C)

hist (CC.MN.3.WP4C$pot )

MPa1.5. <- CC.MN.3.WP4C %>%
           dplyr::filter (pot < -1.5 )

CC.MN.3.WP4C.1 <- CC.MN.3.WP4C %>%
                  dplyr::filter (pot < 0 ) %>%
                  dplyr::filter (pot >= -1.5 ) %>%
                  dplyr::filter (hp.porc < 29 )

write_delim (CC.MN.3.WP4C.1, file = "./Data/procdata/CC.MN.3.WP4C.1.txt",
             delim = ";", na = "NA")


ggscatter (CC.MN.3.WP4C.1 , x = "hp.porc", y = "pot.1", 
           xlim = c(0,35),
           ylim = c(0,1.6),
           color = "gray58", shape = 16, size = 3, # Points color, shape and size
           add = "loess",  # Add regressin line
           add.params = list(color = "black", fill = "lightgray")
           #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
) +
  geom_hline(yintercept = 1.5, lty=1)+
  geom_hline(yintercept = 0, lty=1) +
  geom_vline(xintercept = 17.4, lty=2) +
  geom_vline(xintercept = 7.6, lty=2)  +
  geom_hline(yintercept = 0.3, lty=2) +
  geom_hline(yintercept = 0.15, lty=2) 

#hasta aca la figura de retencion de WP4C ###

# datos curvas de suelo
CC.MN.3.suelos <- read_delim (file = "./Data/rawdata/curvas_retencion_suelos.txt" , 
                            delim ="\t", quote = "\"", 
                            escape_backslash = FALSE,
                            escape_double = TRUE, 
                            col_names = TRUE, 
                            col_types = NULL,
                            locale = default_locale(), 
                            na = "NA")

CC.MN.3.suelos$Muestra <- as.factor (CC.MN.3.suelos$Muestra)


summary(CC.MN.3.suelos)

hist (CC.MN.3.suelos$presion )



write_delim (CC.MN.3.WP4C.1, file = "./Data/procdata/CC.MN.3.WP4C.1.txt",
             delim = ";", na = "NA")


ggscatter (CC.MN.3.suelos , x = "HP", y = "presion",
           xlim = c(0,35),
           ylim = c(0,1.6),
           color = "gray58", shape = 16, size = 3, # Points color, shape and size
           add = "loess",  # Add regressin line
           add.params = list(color = "black", fill = "lightgray")
           #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           #cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
) +
  geom_hline(yintercept = 1.5, lty=1)+
  geom_hline(yintercept = 0, lty=1) 