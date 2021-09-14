##############################################################################
# Codigo para el analisis de los datos de asimilacion neta 
# 
# Gaston Quero - Matias Nion   
#
# 189/01/2021
################################################################################

getwd()

setwd("C:/Users/Usuario/OneDrive/Documentos/Tesis_Maestria_Matias")


# Paquetes 
library (lme4)
library (emmeans)
library ("car")
library ("nlmrt")
library ("easynls")
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
library(plantecophys)

# Se carga los datos para el analisis de las curvas A vs Ci
# Datos de diciembre

An.Ci.dic.sin_out <- read_delim (file = "./Data/procdata/Data.IRGA.dic.sin_out.txt", 
                    delim =";", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = TRUE, 
                    col_types = NULL,
                    locale = default_locale(), 
                    na = "NA")

head (An.Ci.dic.sin_out)

An.Ci.dic.sin_out <- An.Ci.dic.sin_out %>%
                     dplyr::select (c(orden, out,cond.hidrica,clon , maceta, pot, everything() ))


dt=An.Ci.dic.sin_out 
sample = "diciembre"
head (dt)

unique (An.Ci.dic.sin_out$pot )

run.fit.An.Ci <-  function (dt=NULL, sample = NULL, TPU = NULL ){
  
  dir.create (file.path ("Figures", "Plot.fit.An.Ci"), showWarnings = FALSE)
              list.pot <- unique (dt$pot)  
  
  dt.1 <- bind_rows (lapply (list.pot, function (filt.pot){
    
    #filt.pot = "gc_c_22"
    
    print (filt.pot)
    #intrinseca cond  IWUE
    #instantanea         
    
    pot.x <- dt %>%
             dplyr::filter (pot ==  filt.pot) %>%
             dplyr::select (c(Ci, Photo,Cond, Trmmol, "Ci/Ca", Tleaf, PARi)) 
             #dplyr::filter (Ci > 0) %>%
             #dplyr::filter (Photo > 0)
    
    pot.x <- na.omit (pot.x)
    
    #plot ( Photo ~ Ci, data = pot.x )
    
    X1 <- pot.x %>%
          dplyr::filter (Photo == (max (Photo))) %>%
          #dplyr::filter (Ci == (max (Ci))) %>%
          dplyr::mutate (WUEtrnsp = Photo/(1e4*Trmmol))%>%
          dplyr::mutate (WUEcond = Photo/Cond )%>%
          dplyr::mutate (EUR = Photo/PARi)
      
      
     #head (X1)
    if (TPU == TRUE){
      
      fit.x <- fitaci (pot.x, fitTPU=TRUE)
      Amax <- max (fitted (fit.x ))
      tpu <- "TPU"
      
    }
    
    if (TPU == FALSE){
      
      fit.x <- fitaci (pot.x, fitTPU=FALSE)
      Amax <- max (fitted (fit.x ))
      tpu <-"sinTPU"
      
      
    }
    
    Vcmax <-  coef (fit.x) [1]
    Jmax <-  coef (fit.x) [2]
    Rd <- coef (fit.x) [3]
    TPU <- coef (fit.x) [4]
    
    dt.coef <- tibble (pot=filt.pot, X1, Amax =Amax,  Vcmax , Jmax, Rd, TPU, )
    
    dt.coef <- dt.coef %>%
               dplyr::mutate (x= pot)%>%
               tidyr::separate (x,c("clon", "cond.hidrica", "maceta"), sep="_") %>%
               dplyr::select (clon,  cond.hidrica, maceta, everything())
     
    
    png (filename = str_c ("./Figures/Plot.fit.An.Ci/",filt.pot, "_", "fit.An.Ci", tpu ,"_", sample,".png"),
        width = 480, height = 480, units = "px", pointsize = 12)
    
    plot (fit.x, main = str_c (filt.pot,"_",tpu,"-",sample))
    dev.off()
    
    
    
    plot (fit.x, main = str_c (filt.pot,"-",tpu,sample))
    
    
    return (dt.coef)
    
  }))
  
  df.1 <- bind_rows (dt.1 )
  return (df.1)
  
} # aca termina la funcion 

fit.An.Ci.dic.sin_out <- run.fit.An.Ci (dt=An.Ci.dic.sin_out, sample = "diciembre", TPU = FALSE)

fit.An.Ci.dic.sin_out <- bind_cols(periodo = "T0-1", fit.An.Ci.dic.sin_out)

## T12
An.Ci.mar.sin_out <- read_delim (file = "./Data/procdata/Data.IRGA.mar.sin_out.txt", 
                                 delim =";", quote = "\"",
                                 escape_backslash = FALSE,
                                 escape_double = TRUE, 
                                 col_names = TRUE, 
                                 col_types = NULL,
                                 locale = default_locale(), 
                                 na = "NA")

unique (An.Ci.mar.sin_out$pot )



head (An.Ci.mar.sin_out)

An.Ci.mar.sin_out <- An.Ci.mar.sin_out %>%
                     dplyr::select (c(orden, out, cond.hidrica, clon, maceta, pot, everything() ))


gc2 <- An.Ci.mar.sin_out %>%
       dplyr::filter ( pot == "gc_c_2")


gc22 <- An.Ci.mar.sin_out %>%
        dplyr::filter ( pot == "gc_c_22")

plot ( Photo ~ Ci , data=gc2)
plot ( Photo ~ Ci , data=gc22)



unique (An.Ci.mar.sin_out$pot)
An.Ci.mar.sin_out <- An.Ci.mar.sin_out %>%
                     dplyr::filter ( pot != "gc_c_22") %>%
                     dplyr::filter ( pot != "gc_c_2")



fit.An.Ci.mar.sin_out <- run.fit.An.Ci (dt=An.Ci.mar.sin_out, sample = "marzo", TPU = FALSE)

fit.An.Ci.mar.sin_out <- bind_cols(periodo = "T1-2", fit.An.Ci.mar.sin_out)


## aca uno las matrices de los valore calculados por el modelo ####

fit.An.MM.3 <- bind_rows (fit.An.Ci.dic.sin_out, fit.An.Ci.mar.sin_out )

summary (fit.An.MM.3 )

fit.An.MM.3 <- fit.An.MM.3 %>%
               dplyr::select (-TPU)

fit.An.MM.3.a <- fit.An.MM.3 %>%
                 dplyr::group_by (periodo, pot, clon, cond.hidrica, maceta) %>%
                 summarise_all(mean) %>%
                 dplyr::ungroup()


write_delim (fit.An.MM.3.a, file ="./Data/procdata/fit.An.MM.3.a.txt", 
             delim = ";", na = "NA")



      
head ( fit.An.MM.3.a)

ggdotplot (fit.An.MM.3.a, "clon", "Photo", 
                     fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (fit.An.MM.3.a, "clon", "Cond", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (fit.An.MM.3.a, "clon", "Trmmol", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (fit.An.MM.3.a, "clon", "Ci/Ca", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))



ggdotplot (fit.An.MM.3.a, "clon", "WUEtrnsp", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (fit.An.MM.3.a, "clon", "WUEcond", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))



ggdotplot (fit.An.MM.3.a, "clon", "EUR", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (fit.An.MM.3.a, "clon", "Vcmax", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (fit.An.MM.3.a, "clon", "Jmax", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (fit.An.MM.3.a, "clon", "Rd", 
           fill = "cond.hidrica", facet.by = "periodo",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


plot.2 <- ggdotplot (fit.An.Ci.dic.sin_out.TPU, "clon", "Jmax", fill = "cond.hidrica",
           color = "cond.hidrica",
           #ylim=c(min(fit.An.Ci.dic.sin_out.TPU$Jmax)-5,max(fit.An.Ci.dic.sin_out.TPU$Jmax )+5),
           add = "mean_sd",
           palette = c("blue", "red"))

ggarrange(plot.1, plot.2)


ggscatter (fit.An.Ci.dic.sin_out.TPU, x ="Vcmax", y = "Jmax", facet.by ="clon",
          color = "cond.hidrica", fill = "cond.hidrica",
          palette = c("navyblue", "darkorange"),
          shape = 16, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)


ggscatter (fit.An.Ci.dic.sin_out.TPU, x ="Vcmax", y = "Jmax", facet.by ="cond.hidrica",
           color = "clon", fill = "clon",
           palette = c("navyblue", "darkorange", "green", "red"),
           shape = 16, size = 3, # Points color, shape and size
           add = "reg.line",  # Add regressin line
           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
           conf.int = TRUE, # Add confidence interval
           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)


fit.An.Ci.dic.sin_out.sinTPU <- run.fit.An.Ci (dt=An.Ci.dic.sin_out, sample = "diciembre", TPU=FALSE)

plot.1 <- ggdotplot (fit.An.Ci.dic.sin_out.sinTPU , "clon", "Vcmax", fill = "cond.hidrica",
                     color = "cond.hidrica",
                     ylim=c(min(fit.An.Ci.dic.sin_out.TPU$Vcmax)-5,max(fit.An.Ci.dic.sin_out.TPU$Vcmax )+5),
                     add = "mean_sd",
                     palette = c("navyblue", "darkorange"))


plot.2 <- ggdotplot (fit.An.Ci.dic.sin_out.sinTPU, "clon", "Jmax", fill = "cond.hidrica",
                     color = "cond.hidrica",
                     #ylim=c(min(fit.An.Ci.dic.sin_out.TPU$Jmax)-5,max(fit.An.Ci.dic.sin_out.TPU$Jmax )+5),
                     add = "mean_sd",
                     palette = c("blue", "red"))

ggarrange(plot.1, plot.2)


ggscatter (fit.An.Ci.dic.sin_out.sinTPU, x ="Vcmax", y = "Jmax", facet.by ="clon",
           color = "cond.hidrica", fill = "cond.hidrica",
           palette = c("navyblue", "darkorange"),
           shape = 16, size = 3, # Points color, shape and size
           add = "reg.line",  # Add regressin line
           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
           conf.int = TRUE, # Add confidence interval
           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)


ggscatter (fit.An.Ci.dic.sin_out.sinTPU, x ="Vcmax", y = "Jmax", facet.by ="cond.hidrica",
           color = "clon", fill = "clon",
           palette = c("navyblue", "darkorange", "green", "red"),
           shape = 16, size = 3, # Points color, shape and size
           add = "reg.line",  # Add regressin line
           add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
           conf.int = TRUE, # Add confidence interval
           cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
           cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)













write_delim (fit.An.Ci.dic.sin_out.TPU, file ="./Data/procdata/fit.An.Ci.dic.sin_out.TPU.txt", 
             delim = ";", na = "NA")


write_delim (fit.An.Ci.dic.sin_out.sinTPU, file ="./Data/procdata/fit.An.Ci.dic.sin_out.sinTPU.txt", 
             delim = ";", na = "NA")

head (An.Ci.dic.sin_out)
An.Ci.dic.sin_out.400 <- An.Ci.dic.sin_out %>%
                         dplyr::filter (Ci < 500)



fit.An.Ci.dic.sin_out.400.TPU <- run.fit.An.Ci (dt=An.Ci.dic.sin_out.400, sample = "diciembre", TPU = TRUE)



fit.An.Ci.dic.sin_out.sinTPU <- run.fit.An.Ci (dt=An.Ci.dic.sin_out.400, sample = "diciembre", TPU=FALSE)







unique (An.Ci.dic.sin_out$pot)





An.Ci.dic.1 <- An.Ci.dic %>%
             dplyr::mutate (pot = str_c ( clon, cond.hidrica , maceta, sep="_"))
  

length (unique (An.Ci.dic.1$pot))

table (An.Ci.dic.1$pot)


dt= An.Ci.dic.out.1 
obs = 0

 run.plot.An.Ci <-  function (dt=NULL, obs = NULL, sample = NULL){
   
dir.create (file.path ("Figures", "Plots.An.Ci"), showWarnings = FALSE)
 list.pot <- unique (dt$pot)  
   
 dt.1 <- lapply (list.pot, function (filt.pot){
   
   #filt.pot = "gc_d_5"
  
    pot.x <- dt %>%
             dplyr::filter (pot ==  filt.pot) 
    
    ord.x <- as.numeric (row.names(pot.x))
    
    pot.x1 <-  pot.x %>%
              dplyr::mutate (orden = ord.x) %>%
              dplyr::select (orden, everything()) %>%
              dplyr::filter (orden > obs)
    
    id <- unique (pot.x1$pot )
    
    plot.potx1 <- ggscatter (pot.x1, x = "Ci", y = "Photo", 
                             title = id,  
                             add = "loess", 
                             conf.int = TRUE)
    
    print (plot.potx1)
    
    ggexport(plot.potx1, filename = str_c ("./Figures/Plots.An.Ci/",id, "_", "An.Ci","_", sample, ".png"))
    
    return (pot.x1)
    
 })
 
 df.1 <- bind_rows (dt.1 )
 return (df.1)
 
 } # aca termina la funcion  
   
Data.IRGA.dic <- run.plot.An.Ci  ( dt = An.Ci.dic.1 , obs =19)


write_delim (Data.IRGA.dic, file ="./Data/procdata/Data.IRGA.dic.txt", 
             delim = ";", na = "NA")


unique ( Data.IRGA.dic$clon)


An.g2.dic  <-  Data.IRGA.dic %>%
               dplyr::filter (clon == "g2")



# reingreso los datos que pueden estar mal

An.Ci.dic.out <- read_delim (file = "./Data/rawdata/IRGA.diciembre/anco2dic_outliers.txt", 
                         delim ="\t", quote = "\"",
                         escape_backslash = FALSE,
                         escape_double = TRUE, 
                         col_names = TRUE, 
                         col_types = NULL,
                         locale = default_locale(), 
                         na = "NA")


table (An.Ci.dic.out$pot)

str (An.Ci.dic.out )

pot.out <- c(21,37,41,65,85,109,113,133,137,157,229,237,301,305,325,373,397,401,421,433,434,437,445,449,
             481,493,497,517,529,541,545,565)

length(pot.out)

nrow (An.Ci.dic.out) - length(pot.out)


#### aca se filtran los datos segun  matias

An.Ci.dic.out.1 <- An.Ci.dic.out  %>%
                   dplyr::filter (!out %in% pot.out )


table (An.Ci.dic.out.1$pot)


Data.IRGA.dic.sin_out <- run.plot.An.Ci  ( dt = An.Ci.dic.out.1 , obs =0, sample = "sin_out")


write_delim (Data.IRGA.dic.sin_out , file ="./Data/procdata/Data.IRGA.dic.sin_out.txt", 
             delim = ";", na = "NA")




######### An vs PPFD #############

An.ppfd.dic <- read_delim (file = "./Data/rawdata/IRGA.diciembre/anppfddic_GQ.txt", 
                         delim ="\t", quote = "\"",
                         escape_backslash = FALSE,
                         escape_double = TRUE, 
                         col_names = TRUE, 
                         col_types = NULL,
                         locale = default_locale(), 
                         na = "NA")

An.ppfd.dic.1 <- An.ppfd.dic %>%
                 dplyr::mutate (pot = str_c ( clon, cond.hidrica , maceta, sep="_"))


length (unique (An.ppfd.dic.1$pot))

table (An.ppfd.dic.1$pot)

dt = An.ppfd.dic.1 
obs = 19 
sample = "dic_original"


run.plot.An.ppfd <-  function (dt=NULL, obs = NULL, sample = NULL){
  
  dir.create (file.path ("Figures", "Plots.An.PPFD"), showWarnings = FALSE)
  list.pot <- unique (dt$pot)  
  
  dt.1 <- lapply (list.pot, function (filt.pot){
    
    #filt.pot = "g1_c_11"
    
    pot.x <- dt %>%
             dplyr::filter (pot ==  filt.pot) 
    
    ord.x <- as.numeric (row.names(pot.x))
    
    pot.x1 <-  pot.x %>%
      dplyr::mutate (orden = ord.x) %>%
      dplyr::select (orden, everything()) %>%
      dplyr::filter (orden > obs)
    
    id <- unique (pot.x1$pot )
    
    plot.potx1 <- ggscatter (pot.x1, x = "PARi", y = "Photo", 
                             title = id,  
                             add = "loess", 
                             conf.int = TRUE)
    
    print (plot.potx1)
    
    ggexport(plot.potx1, filename = str_c ("./Figures/Plots.An.PPFD/",id, "_", "An.ppfd","_", sample, ".png"))
    
    return (pot.x1)
    
  })
  
  df.1 <- bind_rows (dt.1 )
  return (df.1)
  
} # aca termina la funcion 


Data.IRGA.ppfd.dic <- run.plot.An.ppfd (dt = An.ppfd.dic.1, obs=15, sample = "dic_original")

write_delim (Data.IRGA.ppfd.dic,
             file ="./Data/procdata/Data.IRGA.ppfd.dic.txt", 
             delim = ";", na = "NA")

run.plot.An.time <-  function (dt=NULL, obs = NULL, sample = NULL){
  
  dir.create (file.path ("Figures", "Plots.An.time"), showWarnings = FALSE)
  list.pot <- unique (dt$pot)  
  
  dt.1 <- lapply (list.pot, function (filt.pot){
    
    #filt.pot = "g1_c_16"
    
    pot.x <- dt %>%
      dplyr::filter (pot ==  filt.pot) 
    
    
    
    ord.x <- as.numeric (row.names(pot.x))
    
    
    pot.x1 <-  pot.x %>%
      dplyr::mutate (orden = ord.x) %>%
      dplyr::select (orden, everything()) %>%
      dplyr::filter (orden < obs)
    
    id <- unique (pot.x1$pot )
    
    plot.potx1 <- ggscatter (pot.x1, x = "orden", y = "Photo", 
                             xlim= c(0,20),
                             title = id,  
                             add = "loess", 
                             conf.int = TRUE)
    
    print (plot.potx1)
    
    ggexport(plot.potx1, filename = str_c ("./Figures/Plots.An.time/",id, "_", "An.time","_", sample, ".png"))
    
    return (pot.x1)
    
  })
  
  df.1 <- bind_rows (dt.1 )
  return (df.1)
  
} # aca termina la funcion 


Data.IRGA.timeCi.dic <- run.plot.An.time (dt = An.Ci.dic.out.1, obs=20, sample = "dic_timeCi")


Data.IRGA.timeppfd.dic <- run.plot.An.time (dt = An.ppfd.dic.1, obs=20, sample = "dic_timePPFD")

write_delim (Data.IRGA.timeppfd.dic,
             file ="./Data/procdata/Data.IRGA.time.dic.ppffd.txt", 
             delim = ";", na = "NA")

write_delim (Data.IRGA.timeCi.dic,
             file ="./Data/procdata/Data.IRGA.time.dic.Ci.txt", 
             delim = ";", na = "NA")



### datos que se deberian sacar 

dat.dic.out <- read_delim (file = "./Data/rawdata/IRGA.diciembre/out_data_irga_dic.txt", 
                             delim ="\t", quote = "\"",
                             escape_backslash = FALSE,
                             escape_double = TRUE, 
                             col_names = TRUE, 
                             col_types = NULL,
                             locale = default_locale(), 
                             na = "NA")




length(unique (An.Ci.dic.1$pot))

















## analisis de Quenching

files.param.nion.1.act25 <- dir(str_c("./Data/rawdata/PAM.diciembre/parametro"), pattern = "*.asc")

pre.PAM.param.nion.1.act25.diciembre <- bind_rows (lapply(files.param.nion.1.act25, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/PAM.diciembre/parametro/", filt.raw) , 
                    delim =",", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = TRUE, 
                    col_types = list(col_character (), col_character(),col_character ()),
                    locale = default_locale(), 
                    na = "NA")
  dt <- dt %>%
    dplyr::mutate (pot=filt.raw)
  
}))



#pre.PAM.param$Parameter <- as.factor (pre.PAM.param$Parameter)
param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")


pre.PAM.param.nion.1.act25.diciembre <- pre.PAM.param.nion.1.act25.diciembre %>%
                                    dplyr::filter (Parameter %in% param)%>%
                                    dplyr::select (-X3)

write_delim(pre.PAM.param.nion.1.act25.diciembre, file ="./Data/procdata/pre.PAM.param.nion.1.act25.diciembre.txt" , 
            delim = ",", na = "NA")


### aca reingreso los datos 
# nion.1.act25.
PAM.param.nion.1.act25.diciembre <- read_delim (file ="./Data/procdata/pre.PAM.param.nion.1.act25.diciembre.txt", 
                            delim =",", quote = "\"",
                            escape_backslash = FALSE,
                            escape_double = TRUE, 
                            col_names = TRUE, 
                            col_types = NULL, 
                            locale = default_locale(), 
                            na = "NA")

unique (PAM.param.nion.1.act25.diciembre$pot)

PAM.param.nion.1.act25.diciembre.a <- PAM.param.nion.1.act25.diciembre %>%
                                  dplyr::mutate (x = pot)%>%
                                  tidyr::separate(x , c("pos","genotype","cond.hid", NA), sep="_") %>%
                                  dplyr::mutate (protocol = "nion.1.act25")



diciembre_quenching <- run_quenching_analisis ( dt = PAM.param.nion.1.act25.diciembre.a, sampling ="diciembre",protocol = "nion.1.act25", rep=NULL)

head(diciembre_quenching)
unique ()

diciembre_quenching.5.5 <- diciembre_quenching %>%
                       dplyr::filter (time.min  == 5.5)

diciembre_quenching.5.5.a <- diciembre_quenching.5.5 %>%
                       dplyr::mutate (x = pot)%>%
                       tidyr::separate(x , c("pos","genotype","cond.hid", NA), sep="_") %>%
                       dplyr::mutate (protocol = "nion.1.act25") %>%
                       dplyr::select (-c( ciclo, time.min,protocol, pot, sampling ,pos))%>%
                       dplyr::mutate (cond.hid = factor(cond.hid, levels = c("c","d")))%>%
                       dplyr::mutate (trat = str_c( genotype , cond.hid , sep="." ))

boxplot (phiPS2 ~  genotype * cond.hid, data= diciembre_quenching.5.5.a )

boxplot (phiPS2 ~  genotype + cond.hid, data= diciembre_quenching.5.5.a )

ggboxplot(diciembre_quenching.5.5.a, facet.by = "genotype",
          ylim=c(0,1), x = "cond.hid", y = "phiPS2",
          add = "jitter", fill="cond.hid")

ggboxplot(diciembre_quenching.5.5.a, facet.by = "genotype",
          ylim=c(0,1), x = "cond.hid", y = "phiNPQ",
          add = "jitter", fill="cond.hid")

ggboxplot(diciembre_quenching.5.5.a, facet.by = "genotype",
          ylim=c(0,1), x = "cond.hid", y = "phiNO",
          add = "jitter", fill="cond.hid")



## model 
lm.phiPS2.diciembre <- lm (phiPS2 ~ genotype + cond.hid + genotype * cond.hid, data= diciembre_quenching.5.5.a)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiPS2.diciembre$fitted.values, lm.phiPS2.diciembre$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiPS2.diciembre")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiPS2.diciembre$residuals, main="phiPS2.diciembre",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiPS2.diciembre$residuals),
              sd (lm.phiPS2.diciembre$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiPS2.diciembre$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiPS2.diciembre")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiPS2.diciembre$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiPS2 ~ trat, data = diciembre_quenching.5.5.a)

# ANAVA
(anava.lm.phiPS2.diciembre  <- anova (lm.phiPS2.diciembre))

#0.015486/0.012587
## phiPS2
### en g1
em.phiPS2.g1.diciembre <- emmeans (lm.phiPS2.diciembre, ~ cond.hid, 
                                   at = list (genotype = "g1"))

(cr.phiPS2.g1.diciembre <- contrast (em.phiPS2.g1.diciembre, method = "pairwise"))

cld.phiPS2.g1.diciembre <- cld (em.phiPS2.g1.diciembre, sort=FALSE)

blue.phiPS2.g1.diciembre <- cbind(genotipo= "g1", muestreo="diciembre" ,parametro ="phiPS2",cld.phiPS2.g1.diciembre )

### en g2
em.phiPS2.g2.diciembre <- emmeans (lm.phiPS2.diciembre, ~ cond.hid, 
                               at = list (genotype = "g2"))

(cr.phiPS2.g2.diciembre <- contrast (em.phiPS2.g2.diciembre, method = "pairwise"))

cld.phiPS2.g2.diciembre <- cld (em.phiPS2.g2.diciembre, sort=FALSE)

blue.phiPS2.g2.diciembre <- cbind(genotipo= "g2", muestreo="diciembre" ,parametro ="phiPS2",cld.phiPS2.g2.diciembre )

### en gc
em.phiPS2.gc.diciembre <- emmeans (lm.phiPS2.diciembre, ~ cond.hid, 
                               at = list (genotype = "gc"))

(cr.phiPS2.gc.diciembre <- contrast (em.phiPS2.gc.diciembre, method = "pairwise"))

cld.phiPS2.gc.diciembre <- cld (em.phiPS2.gc.diciembre, sort=FALSE)

blue.phiPS2.gc.diciembre <- cbind(genotipo= "gc", muestreo="diciembre" ,parametro ="phiPS2",cld.phiPS2.gc.diciembre )

### en gt
em.phiPS2.gt.diciembre <- emmeans (lm.phiPS2.diciembre, ~ cond.hid, 
                               at = list (genotype = "gt"))

(cr.phiPS2.gt.diciembre <- contrast (em.phiPS2.gt.diciembre, method = "pairwise"))

cld.phiPS2.gt.diciembre <- cld (em.phiPS2.gt.diciembre, sort=FALSE)

blue.phiPS2.gt.diciembre <- cbind(genotipo= "gt", muestreo="diciembre" ,parametro ="phiPS2",cld.phiPS2.gt.diciembre )


################ ############## ##################333 










quenching.analisis.diciembre.5.5 <- quenching.analisis.diciembre.5.5  %>%
                                 dplyr::mutate(quantum.yield = factor(quantum.yield, levels = c("phiPS2","phiNPQ","phiNO")))%>%
                                 dplyr::mutate (phi = round (phi,2))
                                 dplyr::mutate(quantum.yield = factor(quantum.yield, levels = c("phiPS2","phiNPQ","phiNO")))%>%


g2 <- quenching.analisis.diciembre.5.5 %>%
      dplyr::filter (  genotype == "g2" )




 ggbarplot (quenching.analisis.diciembre.5.5, facet.by = "genotype",
            "cond.hid", "phi",
             title ="idea.1",
             ylab="cond.hid",
             fill = "quantum.yield", color = "quantum.yield", 
             palette = c("gray18","gray38" ,"gray58"),
             label = FALSE, 
             lab.col = "white", lab.pos = "in")


 ggbarplot (quenching.analisis.diciembre.5.5, facet.by = "genotype",
            "cond.hid", "phi",
            title ="idea.1",
            ylab="cond.hid",
            fill = "quantum.yield", color = "quantum.yield", 
            palette = c("gray18","gray38" ,"gray58"),
            label = FALSE, 
            lab.col = "white", lab.pos = "in")
 
 
 
ggexport (phi.plot, filename = str_c("./Figures/Plots.PAM/",sampling,"phi.plot", unique(id$pot), ".tiff"),
          width = 700, height = 500)

## datos de Diciembre
files.diciembre <- dir(str_c("./Data/rawdata/PAM.diciembre/curva"), pattern = "*.asc")
PAM.curves.diciembre <- bind_rows (lapply(files.diciembre, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/PAM.diciembre/curva/", filt.raw) , 
                    delim =",", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = TRUE, 
                    col_types = NULL,
                    locale = default_locale(), 
                    na = "NA")
  dt <- dt %>%
    dplyr::mutate (pot=filt.raw)
  
}))


curves_diciembre <- lapply(files.diciembre, function(filt.pot){
  
  print (filt.pot)
  pot <- PAM.curves.diciembre %>%
    dplyr::filter (pot == filt.pot)
  
  run.plot.curve (dt=pot,sampling = "diciembre" )    
  
  
})


## analisis de Quenching

files.param.nion.1.act25 <- dir(str_c("./Data/rawdata/PAM.diciembre/parametro"), pattern = "*.asc")

pre.PAM.param.nion.1.act25.diciembre <- bind_rows (lapply(files.param.nion.1.act25, function (filt.raw) {
  
  dt <- read_delim (file = str_c("./Data/rawdata/PAM.diciembre/parametro/", filt.raw) , 
                    delim =",", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = TRUE, 
                    col_types = list(col_character (), col_character(),col_character ()),
                    locale = default_locale(), 
                    na = "NA")
  dt <- dt %>%
    dplyr::mutate (pot=filt.raw)
  
}))

#pre.PAM.param$Parameter <- as.factor (pre.PAM.param$Parameter)
param <-  c("Fo","Fm","Fv","Fv/Fm","Fs","Fm'","Fo'")

pre.PAM.param.nion.1.act25.diciembre <- pre.PAM.param.nion.1.act25.diciembre %>%
  dplyr::mutate (Parameter = fct_recode (Parameter, "phi.psII.h"= "èPS2"))

pre.PAM.param.nion.1.act25.diciembre <- pre.PAM.param.nion.1.act25.diciembre %>%
  dplyr::filter (Parameter %in% param)%>%
  dplyr::select (-X3)

write_delim(pre.PAM.param.nion.1.act25.diciembre, path ="./Data/procdata/pre.PAM.param.nion.1.act25.diciembre.txt" , 
            delim = ",", na = "NA")


unique (pre.PAM.param.nion.1.act25.diciembre$pot)

### aca reingreso los datos 
# nion.1.act25.
PAM.param.nion.1.act25.diciembre <- read_delim (file ="./Data/procdata/pre.PAM.param.nion.1.act25.diciembre.txt", 
                                            delim ="\t", quote = "\"",
                                            escape_backslash = FALSE,
                                            escape_double = TRUE, 
                                            col_names = TRUE, 
                                            col_types = NULL, 
                                            locale = default_locale(), 
                                            na = "NA")

unique (PAM.param.nion.1.act25.diciembre$pot)

PAM.param.nion.1.act25.diciembre.a <- PAM.param.nion.1.act25.diciembre %>%
                                      dplyr::mutate (x = pot)%>%
                                      tidyr::separate(x , c("pos","genotype", NA), sep="_") %>%
                                      dplyr::mutate (protocol = "nion.1.act25")



diciembre_quenching <- run_quenching_analisis ( dt = PAM.param.nion.1.act25.diciembre.a, sampling ="diciembre",protocol = "nion.1.act25", rep=NULL)











