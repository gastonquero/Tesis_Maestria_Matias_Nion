###########################################################################
# codigo para el analisis de Datos ensayo de Tesis Maestria
# Ensayo MN.1
# Analisis de datos Gaston Quero - Matias Nion
# 30/10/2019
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
consumo.1 <- read_delim (file = "./Data/rawdata/consumo_imputado.txt" , 
                                                   delim ="\t", quote = "\"", 
                                                   escape_backslash = FALSE,
                                                   escape_double = TRUE, 
                                                   col_names = TRUE, 
                                                   col_types = NULL,
                                                   locale = default_locale(), 
                                                   na = "NA")

consumo.1$maceta <- as.character (consumo.1$maceta)

consumo.1 <- consumo.1 %>%
              dplyr::mutate(pot=str_c (maceta,"." , clon, ".", tratamiento))%>%
              dplyr::mutate (Date = dmy (consumo.1$fecha))%>%
              dplyr::mutate (pesof.1 = pesof + sum.riego) %>%
              dplyr::filter (Date > "2019-11-27") %>%
              dplyr::filter (pesof != 56) %>%
              dplyr::filter (pesoi != 56)

consumo.1$pesof[consumo.1$pesof ==  7.696 ] <-  17.696
consumo.1$pesof.1[consumo.1$pesof.1 ==  7.696 ] <-  17.696

peso.obj <- 17.530 

list.pot <- unique (consumo.1$pot)

pot <- consumo.1 %>%
       dplyr::filter (pot == "1.g2.d")




run.plot.WSPM <- function ( dt = NULL, t1=NULL, peso.obj =NULL){
  
  dir.create (file.path ("Figures", "Plots.WSPM"), showWarnings = FALSE)
  
  dt <- dt
  
  t1.1 <- dmy (t1)
  
  dt.1 <- dt %>%
          dplyr::mutate (dias = Date - t1.1) 

  pesof.1 ~ dias
  wspm.curve <- ggscatter (dt.1 , x = "dias", y = "pesof.1", 
                          title = unique(dt.1$pot), 
                          xlab = "time (d)",
                          ylab = "WSPM (Kg)",
                          point=FALSE) +
    geom_line(color = "gray48", linetype =2, size = 0.5) +
    geom_point(color = "black", size = 1.5) +
    geom_hline(yintercept = mean (dt.1$pesof.1, na.rm = TRUE ), linetype =3, size = 1) #+
    #geom_hline(yintercept = peso.obj, linetype =2, col="red", size = 1) 
  
  
  
  ggexport ( wspm.curve, filename = str_c("./Figures/Plots.WSPM/",unique(dt.1$pot), ".tiff"),
            width = 700, height = 500)
  
  print (wspm.curve)          
}


X.wsp <- lapply(list.pot, function(filt.pot){
  
  print (filt.pot)
  pot <- consumo.1 %>%
         dplyr::filter (pot == filt.pot)
  
  run.plot.WSPM (dt = pot , t1= "28/11/2019")    
  
  
})


head(consumo.1)

## condicion hidrica HP%

CC.MN.1 <- read_delim (file = "./Data/rawdata/macetas_cc.txt" , 
                       delim ="\t", quote = "\"", 
                       escape_backslash = FALSE,
                       escape_double = TRUE, 
                       col_names = TRUE, 
                       col_types = NULL,
                       locale = default_locale(), 
                       na = "NA")
summary (CC.MN.1)

peso.sust.seco <- mean (CC.MN.1$peso.sust.seco, na.rm = TRUE)
peso.bolsa <- mean (CC.MN.1$peso.bolsa, na.rm = TRUE)

consumo.2 <-  consumo.1 %>%
              dplyr::mutate (agua = pesof.1 - peso.sust.seco - peso.bolsa) %>%
              dplyr::mutate (hp.porc = (agua *100)/peso.sust.seco)

consumo.2 %>%
   group_by (tratamiento) %>%
   summarise (max(hp.porc, na.rm=TRUE))


list.pot2 <- unique (consumo.2$pot)

run.plot.hp <- function ( dt = NULL, t1=NULL, hp.ref =NULL){
  
  dir.create (file.path ("Figures", "Plots.HP"), showWarnings = FALSE)
  
  dt <- dt
  
  t1.1 <- dmy (t1)
  
  dt.1 <- dt %>%
          dplyr::mutate (dias = Date - t1.1) 
  

  hp.curve <- ggscatter (dt.1 , x = "dias", y = "hp.porc", 
                           title = unique(dt.1$pot), 
                           ylim=c(0, 25),
                           xlab = "time (d)",
                           ylab = "HP (%)",
                           point=FALSE) +
    geom_line(color = "gray48", linetype =2, size = 0.5) +
    geom_point(color = "black", size = 1.5) +
    geom_hline(yintercept = mean (dt.1$hp.porc, na.rm = TRUE ), linetype =3, size = 1) #+
  #geom_hline(yintercept = peso.obj, linetype =2, col="red", size = 1) 
  
  
  
  ggexport ( hp.curve, filename = str_c("./Figures/Plots.HP/",unique(dt.1$pot), ".tiff"),
             width = 700, height = 500)
  
  print (hp.curve)          
}


X.HP <- lapply(list.pot2, function(filt.pot){
  
  print (filt.pot)
  pot <- consumo.2 %>%
    dplyr::filter (pot == filt.pot)
  
  run.plot.hp (dt = pot , t1= "28/11/2019", hp.ref =NULL)    
  
  
})


summary (consumo.2)
t1= "28/11/2019"

t1.1 <- dmy (t1)

consumo.2 <- consumo.2 %>%
             dplyr::mutate (dias = Date - t1.1)
                 
                 
                 
                 
                 
ggscatter (consumo.2 , x = "dias", y = "hp.porc", facet.by = "clon",
           color="tratamiento",
           palette = c("navyblue", "darkorange"),
           #title = unique(dt.1$pot), 
           ylim=c(0, 25),
           xlab = "time (d)",
           ylab = "HP (%)",
           point=TRUE) 



ggboxplot (consumo.2 , x ="tratamiento", y = "hp.porc", facet.by = "clon",
           color="black", fill= "tratamiento",
           palette = c("darkorange","navyblue" ),
           add = "mean_sd", 
           xlab = "trat",
           ylab = "HP (%)")

ggviolin (consumo.2 , x ="tratamiento", y = "hp.porc", facet.by = "clon",
           color="tratamiento",
           palette = c("blue", "red"),
           add = "jitter", 
           xlab = "trat",
           ylab = "HP (%)")




  
  #geom_line(color = "gray48", linetype =2, size = 0.5) +
  #geom_point(color = "black", size = 1.5) +
  geom_hline(yintercept = mean (hp.porc, na.rm = TRUE ), linetype =3, size = 1) #+

