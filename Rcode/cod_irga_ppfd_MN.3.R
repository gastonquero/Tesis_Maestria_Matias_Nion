##############################################################################
# Codigo para el analisis de los datos de asimilacion neta 
# 
# Gaston Quero - Matias Nion   
#
# 189/01/2021
################################################################################

getwd()

setwd ("R:/Tesis_Maestria_Matias_Nion")


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

# Se carga los datos para el analisis de las curvas A vs Ci
# Datos de diciembre

######### An vs PPFD #############

An.ppfd.dic <- read_delim (file = "./Data/rawdata/anppfddic_GQ.txt", 
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
obs = 15 
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

Data.IRGA.ppfd.dic <-Data.IRGA.ppfd.dic %>%
                     dplyr::arrange (PARi)



Data.IRGA.ppfd.dic.0 <- Data.IRGA.ppfd.dic %>%
                        dplyr::filter (PARi < 199 ) %>%
                        dplyr::mutate (PARx = 0)

Data.IRGA.ppfd.dic.200 <- Data.IRGA.ppfd.dic %>%
                          dplyr::filter (PARi > 199 ) %>%
                          dplyr::filter (PARi < 399 ) %>%
                          dplyr::mutate (PARx = 200)


Data.IRGA.ppfd.dic.400 <- Data.IRGA.ppfd.dic %>%
                          dplyr::filter (PARi > 399 ) %>%
                          dplyr::filter (PARi < 598 ) %>%
                          dplyr::mutate (PARx = 400)


Data.IRGA.ppfd.dic.600 <- Data.IRGA.ppfd.dic %>%
                          dplyr::filter (PARi > 598 ) %>%
                          dplyr::filter (PARi < 798 ) %>%
                          dplyr::mutate (PARx = 600)


Data.IRGA.ppfd.dic.800 <- Data.IRGA.ppfd.dic %>%
                          dplyr::filter (PARi > 798 ) %>%
                          dplyr::filter (PARi < 999 ) %>%
                          dplyr::mutate (PARx = 800)


Data.IRGA.ppfd.dic.1000 <- Data.IRGA.ppfd.dic %>%
                           dplyr::filter (PARi > 999 ) %>%
                           dplyr::filter (PARi < 1198 ) %>%
                           dplyr::mutate (PARx = 1000)


Data.IRGA.ppfd.dic.1200 <- Data.IRGA.ppfd.dic %>%
                           dplyr::filter (PARi > 1198 ) %>%
                           dplyr::mutate (PARx = 1200)


Data.IRGA.ppfd.dic.2 <- bind_rows (Data.IRGA.ppfd.dic.0,
                                   Data.IRGA.ppfd.dic.200,
                                   Data.IRGA.ppfd.dic.400,
                                   Data.IRGA.ppfd.dic.600,
                                   Data.IRGA.ppfd.dic.800,
                                   Data.IRGA.ppfd.dic.1000,
                                   Data.IRGA.ppfd.dic.1200)

Data.IRGA.ppfd.dic.2a <- Data.IRGA.ppfd.dic.2 %>%
                         dplyr::select (c(clon,cond.hidrica,Photo, Cond,Ci, PARx, pot)) 
  


ggline (Data.IRGA.ppfd.dic.2a, x = "PARx", y = "Photo", 
        facet.by= "clon", 
        color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))

ggscatter (Data.IRGA.ppfd.dic.2a, x = "PARx", y = "Photo",
           facet.by= "clon", 
           color = "cond.hidrica",
           palette = c("darkcyan", "firebrick1"))


ggline (Data.IRGA.ppfd.dic.2a, x = "PARx", y = "Cond",
        facet.by= "clon", 
        color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))


ggline (Data.IRGA.ppfd.dic.2a, x = "PARx", y = "Ci",
        facet.by= "clon", 
        color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))


ggscatter (Data.IRGA.ppfd.dic.2a, x = "Cond", y = "Photo",
        facet.by= "clon", 
        color = "cond.hidrica",
        palette = c("darkcyan", "firebrick1"))



gt.dic <- Data.IRGA.ppfd.dic.2 %>%
          dplyr::filter (clon == "gt") %>%
          dplyr::filter (cond.hidrica == "c") %>%
           dplyr::select (c(cond.hidrica,Photo, Cond,Ci, PARx, pot))

unique (gt.dic$pot )
ggline (g1.dic, x = "PARx", y = "Photo",color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))




ggline (g1.dic, x = "PARx", y = "Cond",color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))

ggline (g1.dic, x = "PARx", y = "Ci",color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))


g2.dic <- Data.IRGA.ppfd.dic.2 %>%
          dplyr::filter (clon == "g2") %>%
  #dplyr::filter (cond.hidrica == "c") %>%
          dplyr::select (c(cond.hidrica,Photo, Cond,Ci, PARx, pot))

ggline (g2.dic, x = "PARx", y = "Photo",color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))

ggline (g2.dic, x = "PARx", y = "Cond",color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))

ggline (g2.dic, x = "PARx", y = "Ci",color = "cond.hidrica",
        add = "mean_se", error.plot = "pointrange",  
        palette = c("darkcyan", "firebrick1"))





unique(g1.dic.c$pot )
head(g1.dic.c )
tail(g1.dic.c )

ggscatter (g1.dic.c, x = "PARx", y = "Photo",  
           add = "loess", 
           conf.int = TRUE)

ggline (g1.dic.c, "PARx", "Photo")

ggline (g1.dic.c, x = "PARx", y = "Photo",
        add = "mean_se", error.plot = "pointrange")


        palette = c("darkcyan", "firebrick1"))

ggline(df2, "dose", "len",
       linetype = "supp", shape = "supp")

head (g1.dic)

ggline (g1.dic, "PARi", "Photo")

        linetype = "cond.hidrica", shape = "cond.hidrica")


       






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











