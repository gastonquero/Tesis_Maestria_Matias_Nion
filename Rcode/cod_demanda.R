###########################################################################
# codigo para el analisis de Datos ensayo de Tesis Maestria
# Ensayo invernaculo con deficit 
# Analisis de datos Gaston Quero - Matias Nion
# 27/3/2021
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
library("data.table")
library(ggplot2)
library(ggridges)

#install.packages("modeest")
library(modeest)

#devtools::install_github("tidyverse/lubridate")
## cargar datos ####

sensor.1 <- read_delim (file = "./Data/rawdata/sensor.1.txt" , 
                        delim ="\t", quote = "\"", 
                        escape_backslash = FALSE,
                        escape_double = TRUE, 
                        col_names = TRUE, 
                        col_types = NULL,
                        locale = default_locale(), 
                        na = "NA")

unique (sensor.1$fecha)

sensor.2 <- read_delim (file = "./Data/rawdata/sensor.2.txt" , 
                        delim ="\t", quote = "\"", 
                        escape_backslash = FALSE,
                        escape_double = TRUE, 
                        col_names = TRUE, 
                        col_types = NULL,
                        locale = default_locale(), 
                        na = "NA")

unique (sensor.2$fecha)

sensor.3 <- read_delim (file = "./Data/rawdata/sensor.3.txt" , 
                        delim ="\t", quote = "\"", 
                        escape_backslash = FALSE,
                        escape_double = TRUE, 
                        col_names = TRUE, 
                        col_types = NULL,
                        locale = default_locale(), 
                        na = "NA")

unique (sensor.3$fecha)

### definicion de tiempos ###

t0= "8/11/2019"
t0.1 <- dmy (t0)


t1 = "16/12/2019"
t1.1 <- dmy (t1)

t2= "16/3/2020"
t2.1 <- dmy (t2)


t1.1 - t0.1

t2.1 - t1.1

## verificar la estructura de los datos y cambiar a unidades de tiempo

sensor.1a <- sensor.1 %>%
             dplyr::mutate (fecha.1 = dmy (fecha)) %>%
             dplyr::mutate (hora.1 = hms (hora))  %>%
             dplyr::mutate (mes = lubridate::month (fecha.1, label=TRUE, abbr=FALSE)) %>%
             dplyr::filter (mes != "Abril") %>%
             dplyr::mutate (mes =  fct_relevel (mes,  "Marzo", "Febrero", "Enero","Diciembre",
                                     "Noviembre" ))  %>%
             dplyr::mutate (hora.fact = as.factor (lubridate::hour (hora.1)))

summary (sensor.1a )

sensor.2a <- sensor.2 %>%
             dplyr::mutate (fecha.1 = dmy (fecha)) %>%
             dplyr::mutate (hora.1 = hms (hora))  %>%
             dplyr::mutate (mes = lubridate::month (fecha.1, label=TRUE, abbr=FALSE)) %>%
             dplyr::filter (mes != "Abril") %>%
             dplyr::mutate (mes =  fct_relevel (mes,  "Marzo", "Febrero", "Enero","Diciembre",
                                     "Noviembre" ))  %>%
            dplyr::mutate (hora.fact = as.factor (lubridate::hour (hora.1)))

summary (sensor.2a)


sensor.3a <- sensor.3 %>%
             dplyr::mutate (fecha.1 = dmy (fecha)) %>%
             dplyr::mutate (hora.1 = hms (hora))  %>%
             dplyr::mutate (mes = lubridate::month (fecha.1, label=TRUE, abbr=FALSE)) %>%
             dplyr::filter (mes != "Abril") %>%
             dplyr::mutate (mes =  fct_relevel (mes,  "Marzo", "Febrero", "Enero","Diciembre",
                                     "Noviembre" ))  %>%
            dplyr::mutate (hora.fact = as.factor (lubridate::hour (hora.1)))


summary (sensor.3a)

sensores <- bind_rows (sensor.1a,
                       sensor.2a,
                       sensor.3a)


unique (sensores$fecha.1 )


### selecciono por periodo
periodo_1 <- sensores  %>%
             dplyr::filter (fecha.1 >= t0.1 ) %>%
             dplyr::filter (fecha.1 <= t1.1 ) 

unique (periodo_1$fecha.1 )
head (periodo_1)

periodo_1.t <- periodo_1 %>%
               dplyr::group_by (fecha.1) %>%
               dplyr::summarise (media.temp = round (mean (temp),2) ,
                                 maxima.temp = round (max(temp), 2),
                                 minima.temp = round (min (temp), 2),
                                 media.hum = round (mean (hr),2) ,
                                 maxima.hum = round (max(hr), 2),
                                 minima.hum = round (min (hr), 2))


periodo_1.dpv <- periodo_1.t %>%
                dplyr::summarise ( t.min = round (mean (minima.temp),2),
                                   t.max = round (mean (maxima.temp),2),
                                   t.prom = ((t.max +  t.min ) / 2),
                                   hr.min = round (mean (minima.hum),2),
                                   hr.max = round (mean (maxima.hum),2),
                                   hr.prom = ((hr.max +  hr.min ) / 2)) %>%
                                   dplyr::mutate (pvs.max = 0.6108 * exp ((17.27*t.max)/(t.max + 237.3))) %>%
                                   dplyr::mutate (pvs.min = 0.6108 * exp ((17.27*t.min)/(t.min + 237.3))) %>%
                                   dplyr::mutate (pvs.med = (pvs.max + pvs.min)/2) %>%
                                   dplyr::mutate (pva = ((pvs.min * (hr.max/100)) + (pvs.max * (hr.min/100)))/2) %>%
                                   dplyr::mutate (dpv = pvs.med - pva ) %>%
                                   dplyr::mutate (periodo = "E1") %>%
                                   dplyr::select (periodo, everything())




periodo_2 <- sensores  %>%
             dplyr::filter (fecha.1 > t1.1 ) %>%
             dplyr::filter (fecha.1 <= t2.1 ) 


unique (periodo_2$fecha.1 )
head (periodo_2)

periodo_2.t <- periodo_2 %>%
               dplyr::group_by (fecha.1) %>%
               dplyr::summarise (media.temp = round (mean (temp),2) ,
                    maxima.temp = round (max(temp), 2),
                    minima.temp = round (min (temp), 2),
                    media.hum = round (mean (hr),2) ,
                    maxima.hum = round (max(hr), 2),
                    minima.hum = round (min (hr), 2))


periodo_2.dpv <- periodo_2.t %>%
  dplyr::summarise ( t.min = round (mean (minima.temp),2),
                     t.max = round (mean (maxima.temp),2),
                     t.prom = ((t.max +  t.min ) / 2),
                     hr.min = round (mean (minima.hum),2),
                     hr.max = round (mean (maxima.hum),2),
                     hr.prom = ((hr.max +  hr.min ) / 2)) %>%
                     dplyr::mutate (pvs.max = 0.6108 * exp ((17.27*t.max)/(t.max + 237.3))) %>%
                     dplyr::mutate (pvs.min = 0.6108 * exp ((17.27*t.min)/(t.min + 237.3))) %>%
                     dplyr::mutate (pvs.med = (pvs.max + pvs.min)/2) %>%
                     dplyr::mutate (pva = ((pvs.min * (hr.max/100)) + (pvs.max * (hr.min/100)))/2) %>%
                     dplyr::mutate (dpv = pvs.med - pva ) %>%
                     dplyr::mutate (periodo = "E2") %>%
                     dplyr::select (periodo, everything())

demanda_periodo <- bind_rows (periodo_1.dpv,
                              periodo_2.dpv)


demanda_periodo <- demanda_periodo %>%
                    dplyr::mutate (pvs.max= (round (pvs.max,2)))%>%
                    dplyr::mutate (pvs.min= (round (pvs.min,2))) %>%
                    dplyr::mutate (pvs.med= (round (pvs.med,2))) %>%
                    dplyr::mutate (pva = (round (pva,2))) %>%
                    dplyr::mutate (dpv = (round (dpv,2))) 
  
write_delim (demanda_periodo , file ="./Data/procdata/demanda_periodo.txt", 
             delim = ";", na = "NA")



unique (periodo_2$fecha.1 )

## primera figura de la demanda

unique (sensor.1a$hora.fact)

horas.lux <- as.factor (c( "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"))

##  nos quedamos con el periodo de 9 a 19 horas
sensor.1b  <- sensor.1a  %>%
              dplyr::filter (hora.fact %in% horas.lux ) 


list.mes <- unique (sensor.1b$mes)

dt.temp <- bind_rows (lapply (list.mes, function(filt.mes){
  
  dt.mes <- sensor.1b %>%
    dplyr::filter ( mes == filt.mes)
  
  esta <-  dt.mes %>%
    dplyr::group_by (hora.fact) %>%
    dplyr::summarise (media.1 = round (mean (temp),2) ,
                      #moda.2 = round (mfv (temp), 2),
                      maxima = round ( max(temp), 2),
                      minima = round (min (temp), 2))%>%
    dplyr::mutate (temp.promedio.3 = (maxima + minima) / 2) %>%
    dplyr::mutate (dif.1.3  = media.1 - temp.promedio.3 ) %>%
    #dplyr::mutate (dif.1.2  = media.1 - moda.2) %>%
    dplyr::ungroup() %>%
    dplyr::mutate (mes = unique (dt.mes$mes)) %>%
    dplyr::select (mes, everything())
  
  
  return (esta)
  
}))

dt.hr <- bind_rows (lapply (list.mes, function(filt.mes){
  
  dt.mes <- sensor.1b %>%
    dplyr::filter ( mes == filt.mes)
  
  esta <- hr.hora <- dt.mes %>%
    dplyr::group_by (hora.fact) %>%
    dplyr::summarise (media.1 = round (mean (hr),2) ,
                      #moda.2 = round (mfv (hr), 2),
                      maxima = round ( max(hr), 2),
                      minima = round (min (hr), 2))%>%
    dplyr::mutate (hr.promedio.3 = (maxima + minima) / 2) %>%
    dplyr::mutate (dif.1.3  = media.1 - hr.promedio.3 ) %>%
    #dplyr::mutate (dif.1.2  = media.1 - moda.2) %>%
    dplyr::ungroup() %>%
    dplyr::mutate (mes = unique (dt.mes$mes)) %>%
    dplyr::select (mes, everything())
  
  
  return (esta)
  
}))








ggplot(sensor.1b, aes(x = temp, y = hora.fact, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  theme_minimal()


ggplot(sensor.1b, aes(x =  hr, y = hora.fact, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  theme_minimal()





list.mes <- unique (sensor.1b$mes)

lapply (list.mes, function(filt.mes){
  
  dt.mes <- sensor.1b %>%
            dplyr::filter ( mes == filt.mes)
  
#table ( dt.mes$fecha)
  
  ggplot (dt.mes  , aes(x = temp, y = hora.fact, facet.by = mes, 
                         fill = 0.5 - abs(0.5 - stat(ecdf)))) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
    scale_fill_viridis_c(name = "Tail probability", direction = -1) +
    labs(title = filt.mes) + 
    theme_minimal()
  

  
})


lapply (list.mes, function(filt.mes){
  
  dt.mes <- sensor.1b %>%
    dplyr::filter ( mes == filt.mes)
  
  #table ( dt.mes$fecha)
  
  ggplot (dt.mes  , aes(x = hr, y = hora.fact, facet.by = mes, 
                        fill = 0.5 - abs(0.5 - stat(ecdf)))) +
    stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
    scale_fill_viridis_c(name = "Tail probability", direction = -1) +
    labs(title = filt.mes) + 
    theme_minimal()

})


lapply (list.mes, function(filt.mes){
  
  dt.mes <- sensor.1b %>%
    dplyr::filter ( mes == filt.mes)
  
  #table ( dt.mes$fecha)
  
  
  plot.1 <- ggboxplot (dt.mes, x = "hora.fact", y = "temp", title = filt.mes,
             notch = TRUE)
  
  ggexport( plot.1 , filename =str_c("./Figures/", "boxtemp", filt.mes, ".png"))
  
  
  # Add violin + mean_sd
  
  plot.2 <- ggviolin(dt.mes, x = "hora.fact",  y = "temp" ,title = filt.mes,
           add = "mean_sd")
  
  
  ggexport (plot.2 , filename =str_c("./Figures/", "violintemp", filt.mes, ".png"))
  
  
})




write_excel_csv (dt.temp  , file ="./Data/procdata/dt.temp.csv",
                 na = "NA",
                 append = FALSE,
                 delim = ";",
                 quote_escape = "double",
                 eol = "\n" )



write_delim (dt.temp , file ="./Data/procdata/dt.temp.txt", 
             delim = ";", na = "NA")

## hr

lapply (list.mes, function(filt.mes){
  
  dt.mes <- sensor.1b %>%
    dplyr::filter ( mes == filt.mes)
  
  #table ( dt.mes$fecha)
  
  
  plot.1 <- ggboxplot (dt.mes, x = "hora.fact", y = "hr", title = filt.mes,
                       notch = TRUE)
  
  ggexport( plot.1 , filename =str_c("./Figures/", "boxhr", filt.mes, ".png"))
  
  
  # Add violin + mean_sd
  
  plot.2 <- ggviolin(dt.mes, x = "hora.fact",  y = "hr" ,title = filt.mes,
                     add = "mean_sd")
  
  
  ggexport (plot.2 , filename =str_c("./Figures/", "violinhr", filt.mes, ".png"))
  
  
})


dt.hr <- bind_rows (lapply (list.mes, function(filt.mes){
  
  dt.mes <- sensor.1b %>%
    dplyr::filter ( mes == filt.mes)
  
  esta <- hr.hora <- dt.mes %>%
    dplyr::group_by (hora.fact) %>%
    dplyr::summarise (media.1 = round (mean (hr),2) ,
                      #moda.2 = round (mfv (hr), 2),
                      maxima = round ( max(hr), 2),
                      minima = round (min (hr), 2))%>%
    dplyr::mutate (hr.promedio.3 = (maxima + minima) / 2) %>%
    dplyr::mutate (dif.1.3  = media.1 - hr.promedio.3 ) %>%
    #dplyr::mutate (dif.1.2  = media.1 - moda.2) %>%
    dplyr::ungroup() %>%
    dplyr::mutate (mes = unique (dt.mes$mes)) %>%
    dplyr::select (mes, everything())
  
  
  return (esta)
  
}))



write_excel_csv (dt.hr  , file ="./Data/procdata/dt.hr.csv",
                 na = "NA",
                 append = FALSE,
                 delim = ";",
                 quote_escape = "double",
                 eol = "\n" )



write_delim (dt.hr , file ="./Data/procdata/dt.hr.txt", 
             delim = ";", na = "NA")



# unifico los datos de temperatura y humedad


dt.demanda <-  dt.temp %>%
               dplyr::inner_join ( dt.hr, by=c("mes","hora.fact" ))



dt.demanda.1 <- dt.demanda %>%
                dplyr::select (mes, hora.fact, media.1.x, media.1.y ) %>%
                dplyr::rename (Temp.med = media.1.x, HR = media.1.y ) %>%
                dplyr::mutate (hora.time = as.character(hora.fact)) %>%
                dplyr::mutate (hora.x = str_c (hora.time, ":00")) %>%
                dplyr::mutate (hora.x1 = hm (hora.x)) %>%
                dplyr::select (mes, hora.x1, Temp.med,   HR) %>%
                dplyr::rename (hora = hora.x1) %>%
                dplyr::mutate (pvs = 0.6108 * exp ((17.27*Temp.med)/(Temp.med + 237.3))) %>%
                dplyr::mutate (pva = (pvs * HR)/100 ) %>%
                dplyr::mutate (vpd = pvs - pva )



dt.demanda.2 <- dt.demanda.1 %>%
                dplyr::group_by ( mes) %>%
                dplyr::summarise (media.vpd  = round (mean (vpd),2) ,
                    max.vpd = round ( max(vpd), 2),
                    min.vpd = round (min (vpd), 2),
                    media.temp  = round (mean (Temp.med),2),
                    max.temp = round ( max(Temp.med), 2),
                    min.temp = round (min (Temp.med), 2),
                    media.hr  = round (mean (HR),2),
                    max.hr = round (max (HR), 2),
                    min.hr = round (min (HR), 2))%>%
               dplyr::ungroup() %>%
               dplyr::mutate (mes = unique (dt.demanda.1$mes)) %>%
               dplyr::select (mes, everything())


  
summary (dt.demanda.2)


ggplot(dt.demanda.1, aes(x = Temp.med, y = mes, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  theme_minimal()

ggplot(dt.demanda.1, aes(x = HR, y = mes, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  theme_minimal()
                           
ggplot(dt.demanda.1, aes(x = vpd, y = mes, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  theme_minimal()


###
t1= "28/11/2019"

t1.1 <- dmy (t1)


sensor.1c <- sensor.1b %>%
             dplyr::filter ( fecha.1 >= t1.1 ) %>%
             dplyr::mutate (dias = fecha.1 - t1.1) %>%
             dplyr::mutate (pvs = 0.6108 * exp ((17.27*temp)/(temp + 237.3))) %>%
             dplyr::mutate (pva = (pvs * hr)/100 ) %>%
             dplyr::mutate (vpd = pvs - pva ) %>%
             dplyr::filter (hora.fact %in% horas.lux ) 

list.fecha <- (unique (sensor.1c$fecha))

dt.fecha <- bind_rows (lapply (list.fecha, function(filt.fecha){
  
  dt.1 <- sensor.1c%>%
           dplyr::filter (fecha == filt.fecha)
  
dt.2 <-  dt.1 %>%
         dplyr::group_by (fecha) %>%
         dplyr::summarise (media.vpd  = round (mean (vpd),2) ,
                      max.vpd = round ( max(vpd), 2),
                      min.vpd = round (min (vpd), 2),
                      media.temp  = round (mean (temp),2),
                      maxima = round ( max(temp), 2),
                      minima = round (min (temp), 2),
                      media.hr  = round (mean (hr),2),
                      maxima = round (max(hr), 2),
                      minima = round (min (hr), 2))%>%
    dplyr::ungroup() %>%
    dplyr::mutate (fecha.1 = unique (dt.1$fecha.1)) %>%
    dplyr::select (fecha.1, everything())
  
  
  return (dt.2)
  
}))



dt.fecha.1 <- dt.fecha %>%
              dplyr::mutate (dias = fecha.1 - t1.1) 




ggscatter (dt.fecha.1 , x = "dias", y = "media.vpd", 
           title = "vpd", 
           ylim=c(0, 10),
           xlab = "time (d)",
           ylab = "vpd",
           point=FALSE) +
  geom_line (color = "gray48", linetype =2, size = 0.5) +
  geom_point (color = "black", size = 1.5) +
  geom_hline (yintercept = mean (dt.fecha.1$media.vpd, na.rm = TRUE ), col="red", linetype =2, size = 1.5) 

%>%
  geom_hline (yintercept = mean(dt.fecha.1$max.vpd, na.rm = TRUE), col="blue", linetype =1, size = 2) %>%
  geom_hline (yintercept = mean(dt.fecha.1$min.vpd, na.rm = TRUE), col="orange", linetype =1, size = 2) 
  




ggscatter (sensor.1b , x = "dias", y = "vpd", 
           title = "vpd", 
           #ylim=c(0, 25),
           xlab = "time (d)",
           ylab = "vpd (kPa)",
           point=FALSE) +
  geom_line(color = "gray48", linetype =2, size = 0.5) +
  geom_point(color = "black", size = 1.5) +
  geom_hline(yintercept = mean (sensor.1b$vpd, na.rm = TRUE ), col="red", linetype =1, size = 2)


ggscatter (sensor.1b , x = "dias", y = "dewpt", 
           title = "dwp", 
           #ylim=c(0, 25),
           xlab = "time (d)",
           ylab = "dwp",
           point=FALSE) +
  geom_line(color = "gray48", linetype =2, size = 0.5) +
  geom_point(color = "black", size = 1.5) +
  geom_hline(yintercept = mean (sensor.1b$dewpt, na.rm = TRUE ), col="red", linetype =1, size = 2)



ggscatter (sensor.1b , x = "dias", y = "temp", 
           title = "dwp", 
           #ylim=c(0, 25),
           xlab = "time (d)",
           ylab = "C",
           point=FALSE) +
  geom_line(color = "gray48", linetype =2, size = 0.5) +
  geom_point(color = "black", size = 1.5) +
  geom_hline(yintercept = mean (sensor.1b$temp, na.rm = TRUE ), col="red", linetype =1, size = 2)




ggscatter (consumo.2 , x = "dias", y = "hp.porc", facet.by = "clon",
           color="tratamiento",
           palette = c("navyblue", "darkorange"),
           #title = unique(dt.1$pot), 
           ylim=c(0, 25),
           xlab = "time (d)",
           ylab = "HP (%)",
           point=TRUE) 


 #+
#geom_hline(yintercept = peso.obj, linetype =2, col="red", size = 1) 









x  <- hm (9:00)

xx <- c("09:10", "09:02", "1:10")
ms(xx)
ms("7 6")
ms("6,5")
hm(c("09:10", "09:02", "1:10"))
hm("7 6")
hm("6,5")















dt.temp <- bind_rows (lapply (list.mes, function(filt.mes){
  
  dt.mes <- sensor.1b %>%
    dplyr::filter ( mes == filt.mes)
  
  esta <- temp.hora <- dt.mes %>%
    dplyr::group_by (hora.fact) %>%
    dplyr::summarise (media.1 = round (mean (temp),2) ,
                      moda.2 = round (mfv (temp), 2),
                      maxima = round ( max(temp), 2),
                      minima = round (min (temp), 2))%>%
    dplyr::mutate (temp.promedio.3 = (maxima + minima) / 2) %>%
    dplyr::mutate (dif.1.3  = media.1 - temp.promedio.3 ) %>%
    dplyr::mutate (dif.1.2  = media.1 - moda.2) %>%
    dplyr::ungroup() %>%
    dplyr::mutate (mes = unique (dt.mes$mes)) %>%
    dplyr::select (mes, everything())
  
  
  return (esta)
  
}))



write_excel_csv (dt.temp  , file ="./Data/procdata/dt.temp.csv",
                 na = "NA",
                 append = FALSE,
                 delim = ";",
                 quote_escape = "double",
                 eol = "\n" )



write_delim (dt.temp , file ="./Data/procdata/dt.temp.txt", 
             delim = ";", na = "NA")





ggplot(sensor.1a, aes(x = temp, y = mes)) +
         stat_density_ridges(quantile_lines = TRUE, quantiles = 2)


ggplot(sensor.1a, aes(x = temp, y = hora.fact)) +
     stat_density_ridges(quantile_lines = TRUE, quantiles = 2)




ggplot(sensor.1a, aes(x = temp, y = mes, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
theme_minimal()







  
unique (sensor.1a$hora.fact )





stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.975), alpha = 0.7)


  scale_fill_viridis_d(name = "Quartiles")


  iris.tr %>%
    mutate(Species = fct_reorder(Species, mSW)) %>%
    ggplot() +
    aes(Species, mSW, color = Species) +
    geom_point()

  sizes <- factor(sizes, levels = c("small", "medium", "large"))


+ geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.9, alpha = 0.7,
    vline_size = 1, vline_color = "red",
    point_size = 0.4, point_alpha = 1,
    position = position_raincloud(adjust_vlines = TRUE)
  
 


mes.1 <- lubridate::month (dmy (sensor.1$fecha), label = TRUE)


month(ymd(080101))
month(ymd(080101), label = TRUE)
month(ymd(080101), label = TRUE, abbr = FALSE)
month(ymd(080101) + months(0:11), label = TRUE)

month(ymd(080101), label = TRUE)
unique (sensor.1$fecha )



sensor.2 <- read_delim (file = "./Data/rawdata/sensor.2.txt" , 
                        delim ="\t", quote = "\"", 
                        escape_backslash = FALSE,
                        escape_double = TRUE, 
                        col_names = TRUE, 
                        col_types = NULL,
                        locale = default_locale(), 
                        na = "NA")

unique (sensor.2$fecha )



sensor.3 <- read_delim (file = "./Data/rawdata/sensor.3.txt" , 
                        delim ="\t", quote = "\"", 
                        escape_backslash = FALSE,
                        escape_double = TRUE, 
                        col_names = TRUE, 
                        col_types = NULL,
                        locale = default_locale(), 
                        na = "NA")

fechas.3 <- read_delim (file = "./Data/rawdata/fechas.3.txt" , 
                        delim ="\t", quote = "\"", 
                        escape_backslash = FALSE,
                        escape_double = TRUE, 
                        col_names = TRUE, 
                        col_types = NULL,
                        locale = default_locale(), 
                        na = "NA")

sensor.3.x1 <- sensor.3 %>%
               dplyr::filter (fecha %in% fechas.3$fecha)

unique (sensor.3.x1$fecha )



sensor.3.x2 <- sensor.3 %>%
               dplyr::filter (!fecha %in% fechas.3$fecha)

unique (sensor.3.x2$fecha )





####3 
sensor.1$fecha

sensor.1a <- sensor.1 %>%
             dplyr::mutate (date = dmy (fecha)) %>%
             dplyr::mutate (hours = hms (hora)) 



unique (sensor.1a$date)

unique (sensor.1a$hours)


sensor.2a <- sensor.2 %>%
            dplyr::mutate (date = dmy (fecha))

unique (sensor.2a$date)


sensor.3.x1a <- sensor.3.x1 %>%
                dplyr::mutate (date = mdy (fecha))


sensor.3.x2a <- sensor.3.x2 %>%
                dplyr::mutate (date = dmy (fecha))

sensor.3a <- bind_rows (sensor.3.x2a ,sensor.3.x1a )

unique (sensor.3a$date )



summary (sensor.1a$date)

sensor.1a
# Plot a subset of the data
ss1 <- subset (sensor.1a, date <= as.Date("2020-03-22"))


dt =ss1

t1= "8/11/2019"
dt <- dt

h1 = "00:00:00"

t1.1 <- dmy (t1)

h1.1 <- hms (h1)

dt.1 <- dt %>%
        dplyr::mutate (dias = date - t1.1) %>%
        dplyr::mutate (hrs = hours - h1.1)

pesof.1 ~ dias
 ggscatter (dt.1 , x = "dias", y = "temp", 
                         title = "s1", 
                         xlab = "time (d)",
                         ylab = "Temp",
                         point=FALSE) +
  geom_line(color = "gray48", linetype =2, size = 0.5) +
  #geom_point(color = "gray", size = 1.5) +
  geom_hline(yintercept = mean (dt.1$temp, na.rm = TRUE ), linetype =2, size = 1, color="red") +
   geom_hline(yintercept = median (dt.1$temp, na.rm = TRUE ), linetype =2, size = 1, color="blue") 
 
 
 ggscatter (dt.1 , x = "hora", y = "temp", 
            title = "s1", 
            xlab = "time (h)",
            ylab = "Temp",
            point=FALSE) +
   #geom_line(color = "gray48", linetype =2, size = 0.5) +
   geom_point(color = "Black", size = 1.5) +
   geom_hline(yintercept = mean (dt.1$temp, na.rm = TRUE ), linetype =2, size = 1, color="red") +
   geom_hline(yintercept = median (dt.1$temp, na.rm = TRUE ), linetype =2, size = 1, color="blue") 
 
 
 ggboxplot (dt.1 , x = "hours", y = "temp", 
            title = "s1", 
            xlab = "time (h)",
            ylab = "Temp") 
 
 ggplot(dt.1) + geom_boxplot(aes(x = hora, y = temp, group = hora))
 
 
 +
   #geom_line(color = "gray48", linetype =2, size = 0.5) +
   geom_point(color = "Black", size = 1.5) +
   geom_hline(yintercept = mean (dt.1$temp, na.rm = TRUE ), linetype =2, size = 1, color="red") +
   geom_hline(yintercept = median (dt.1$temp, na.rm = TRUE ), linetype =2, size = 1, color="blue") 
 
 
 
 
 #+
#geom_hline(yintercept = peso.obj, linetype =2, col="red", size = 1) 

dt <- ss1

t1.1 <- dmy (t1)

dt.1 <- dt %>%
        dplyr::mutate (dias = date - 2019-11-08) 

pesof.1 ~ dias
wspm.curve <- ggscatter (dt.1 , x = "dias", y = "temp")
, 
                         title = unique(dt.1$pot), 
                         xlab = "time (d)",
                         ylab = "WSPM (Kg)",
                         point=FALSE) +
  geom_line(color = "gray48", linetype =2, size = 0.5) +
  geom_point(color = "black", size = 1.5) +
  geom_hline(yintercept = mean (dt.1$pesof.1, na.rm = TRUE ), linetype =3, size = 1) #+
#geom_hline(yintercept = peso.obj, linetype =2, col="red", size = 1) 






ss1a  <- ss1 %>%
        dplyr::mutate (mes = month(as.POSIXlt(date, format="%d/%m/%Y"))) %>%
        dplyr::mutate (mes.1 = as.factor (mes ))
      

summary (ss1a )  
ss1b  <- na.omit ( ss1a )  
        
 ggboxplot (ss1b, x = "mes.1", y="temp",
                     add = "mean", rug = TRUE,
                     color = "mes", fill = "mes",
                     palette = c("#00AFBB", "#E7B800", "red"))   
        
        
            
ggplot (data = ss1, aes(x =  hora, y = temp)) + 
          #geom_point(color = "black", size = 1) + 
          geom_bar(stat="identity") + 
          scale_x_datetime(limits =c(mdy_hms("10/2/16 20:00:00"),mdy_hms("10/3/16 20:00:00")))        
        
        
        
   gghistogram (ss1a, x = "hora",
                    add = "mean", rug = TRUE,
                    color = "mes", fill = "mes",
                    palette = c("#00AFBB", "#E7B800", "red"))    
        


unique (ss1a$mes )

ss2 <- subset (sensor.2a, date <= as.Date("2020-03-22"))




ss3 <- subset (sensor.3a, date <= as.Date("2020-03-22"))






ggplot (data = ss1, aes(x =  hora, y = temp)) + 
        geom_point(aes (color = "mes", size = 1))

       geom_point(color = "black", size = 1)





ps2 <- ggplot (data = ss2, aes(x = date, y = temp)) + 
         geom_point(color = "blue", size = 1)

ps3 <-  ggplot (data = ss3, aes(x = date, y = temp)) + 
         geom_point(color = "black", size = 1)

ggarrange (ps1,ps2,ps3, ncol=1, nrow=3  )


ss <- bind_rows(ss1, ss2, ss3)

ggplot (ss, aes(x = date, y = temp)) + 
        geom_point(aes(color = sensor), size = 0.5) +
        scale_color_manual(values = c("darkred", "navyblue", "gray48")) +
  theme_minimal()


# Base plot with date axis
p <- ggplot(data = ss1, aes(x = date, y = temp)) + 
     geom_point(color = "darkred", size = 1)


p
# Set axis limits c(min, max)
min <- NA
max <- as.Date ("2020-03-22")
p + scale_x_date(date_labels = "%B/%Y")



p + stat_smooth(
  color = "#FC4E07", fill = "#FC4E07",
  method = "loess"
)


ggplot (data = ss2, aes(x = date, y = temp)) + 
  geom_point(color = "blue", size = 1)





# unificar ###
sensores <- bind_rows (sensor.1 , sensor.2 , sensor.3)


###

ggboxplot (sensores , x = "fecha", y = "temp", 
           #title = unique(dt.1$pot), 
           #xlab = "time (d)",
           #ylab = "WSPM (Kg)",
           point=TRUE)

dev.off()
devoff()
+
  geom_line(color = "gray48", linetype =2, size = 0.5) +
  geom_point(color = "black", size = 1.5) +
  geom_hline(yintercept = mean (dt.1$pesof.1, na.rm = TRUE ), linetype =3, size = 1) #+

#geom_hline(yintercept = peso.obj, linetype =2, col="red", size = 1) 




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

