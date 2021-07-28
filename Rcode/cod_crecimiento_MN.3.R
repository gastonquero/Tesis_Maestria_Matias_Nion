######################################################################################
# Codigo para el analisis de los Datos de crecimiento de la Tesis de Maestria     #
#  Matias Nion                                                                      #
#                                                                                    #
# Gaston Quero - Matias Nion - Sebastian Simondi                                                   #          
# 26/7/2021                                                                          #
######################################################################################

#####
# los datos correspoden a .......

getwd ()
setwd ( "R:/Tesis_Maestria_Matias_Nion")

#setwd("C:/Users/Usuario/OneDrive/Documentos/Paper_fotosintesis_dinámica")

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
library("Hmisc")
library("PerformanceAnalytics")
library (GGally)
library(ggcorrplot)


### seteo del directorio donde estan los datos

crecimiento.MN.3 <- read_delim (file ="./Data/rawdata/crecimiento.txt", 
                             delim ="\t", quote = "\"", escape_backslash = FALSE,
                             escape_double = TRUE, col_names = TRUE, col_types = NULL,
                             na = "NA")

head (crecimiento.MN.3)

crecimiento.MN.3a <- crecimiento.MN.3 %>%
                     dplyr::mutate (ensayo ="MN.3") %>%
                     dplyr::mutate (ambiente ="Inv") %>%
                     dplyr::rename (genotipo = clon)  %>%
                     dplyr::rename (cond.hidr = cond.hidrica) %>%
                     dplyr::rename (pos = maceta) %>%
                     dplyr::select (c(ensayo, ambiente, genotipo,
                                      cond.hidr, pos, everything()))

crecimiento.MN.3b <- crecimiento.MN.3a %>%
                    dplyr::mutate (date = dmy (fecha)) %>%
                    dplyr::mutate (pot = str_c (pos,genotipo,cond.hidr, sep="."))

# primera medida de altura
# corte de agua
t0 = "12/11/2019"

# segunda medida de altura
# fin de corto plazo

t1 = "20/12/2019"

## tercera medida de altura
# fin de largo plazo

t2 = "22/03/2020"

# pasar a formato date
t0 <- dmy (t0)
t1 <- dmy (t1)
t2 <- dmy (t2)
  
  
listapot <- unique(crecimiento.MN.3b$pot)
  
  df.TC <- bind_rows ((lapply (listapot, function(filt.pot){
    
    #filt.pot = "3.g1.d"
    print (filt.pot)
    
    dt.pot <- crecimiento.MN.3b %>%
              dplyr::filter (pot== filt.pot)
    
    id.pot <- unique (dt.pot$pot)
    
    W1 <- dt.pot  %>%
          dplyr::select (date , alt) %>%
          dplyr::filter (date == t0 )%>%
          dplyr::rename (alt.0 = alt)
    
    W2 <- dt.pot %>%
          dplyr::select (date, alt)  %>%
          dplyr::filter (date == t1) %>%
          dplyr::rename (date.ST = date)%>%
          dplyr::rename (alt.ST = alt)
    
    W3 <- dt.pot %>%
          dplyr::select (date, alt)%>%
          dplyr::filter (date == t2) %>%
          dplyr::rename (date.LT = date)%>%
          dplyr::rename (alt.LT = alt)
    
    
    dt.TC <- bind_cols (W1, W2, W3) 
    
    dt.TC <- dt.TC %>%
             dplyr::mutate (delta.ST = (date.ST - date)) %>%
             dplyr::mutate (delta.LT = (date.LT - date)) %>%
             dplyr::mutate (delta.LT_ST = (date.LT - date.ST)) %>%
             dplyr::mutate (delta.STnum = as.numeric (delta.ST)) %>%
             dplyr::mutate (delta.LTnum = as.numeric (delta.LT)) %>%
             dplyr::mutate (delta.LT_STnum = as.numeric (delta.LT_ST)) %>%
             dplyr::mutate (TCA.ST = (alt.ST - alt.0)/delta.STnum) %>%
             dplyr::mutate (TCA.LT = (alt.LT - alt.0)/delta.LTnum) %>%
             dplyr::mutate (TCA.ST_LT = (alt.LT - alt.ST)/delta.LT_STnum) %>%
             dplyr::mutate (TCR.ST = (log(alt.ST) - log(alt.0))/delta.STnum)%>%
             dplyr::mutate (TCR.LT = (log(alt.LT) - log(alt.0))/delta.LTnum) %>%
             dplyr::mutate (TCR.ST_LT = (log(alt.LT) - log(alt.ST))/delta.LT_STnum)
      
    dt.x <- dt.pot [1,] %>%
            dplyr::select (c(ensayo, ambiente, genotipo, cond.hidr,   pos, ))
    

    dt.TCA <- dt.TC  %>%
              dplyr::select (c( TCA.ST, TCA.LT,
                             TCA.ST_LT))
    
    dt.TCA.W <- bind_cols (dt.x, pot=id.pot,   dt.TCA )
    
    dt.TCA.L <- dt.TCA.W  %>% 
               pivot_longer(
        cols = starts_with("TC"),
        names_to = "duracion",
        names_prefix = "TCA.",
        values_to = "TCA",
        values_drop_na = TRUE)
    
    dt.TCR <- dt.TC  %>%
              dplyr::select (c( TCR.ST, TCR.LT,TCR.ST_LT))
    
    dt.TCR.W <- bind_cols (pot=id.pot,dt.TCR )
    
    dt.TCR.L <- dt.TCR.W  %>% 
      pivot_longer(
        cols = starts_with("TC"),
        names_to = "duracion",
        names_prefix = "TCR.",
        values_to = "TCR",
        values_drop_na = TRUE)
    
    dt.TC.1 <- dt.TCA.L %>%
              dplyr::inner_join( dt.TCR.L, by=c("pot","duracion") )
    
    return (dt.TC.1)
    })))
    
  head (df.TC)
  
  df.TC.ST_LT <- df.TC %>%
                 dplyr::filter (duracion != "LT") %>%
                 dplyr::mutate (cond.term = str_c (duracion,cond.hidr, sep="."))
  
  
  ggboxplot (df.TC.ST_LT, "genotipo", "TCR",
             facet.by = "duracion",
             add = "mean_sd", 
             #fill="cond.hidr", 
             color = "cond.hidr",
             palette = c("darkorange", "navyblue"))
  
  ggdotplot (df.TC.ST_LT, "genotipo", "TCA",
             facet.by ="cond.hidr",
             add = "mean_sd", 
             fill="duracion", color = "duracion",
             palette = c("gray48", "red"))
  
  ggdotplot (df.TC.ST_LT, "genotipo", "TCA",
             add = "mean_sd", 
             fill="cond.term", color = "cond.term",
             palette = c("darkorange", "navyblue", "black", "green"))
  
  
  
  
  
  
  unique (  df.TC.ST_LT$cond.term )
  
  
 
  ggdotplot (df.TC.ST_LT, "genotipo", "TCR",
             add = "mean_sd", 
             fill="cond.term", color = "cond.term",
             palette = c("darkorange", "navyblue", "black", "green"))
  
  
  
  ggboxplot (df.TC.ST_LT, "duracion", "TCA",facet.by = "genotipo"
             , color = "cond.hidr",
            palette = c("#00AFBB", "#E7B800")) 
  
  ggboxplot (df.TC.ST_LT, "duracion", "TCR",facet.by = "genotipo"
             , color = "cond.hidr",
             palette = c("#00AFBB", "#E7B800")) 
  
  ggboxplot (df.TC.ST_LT, "genotipo", "TCR",facet.by = "duracion"
             , color = "cond.hidr",
             palette = c("#00AFBB", "#E7B800")) 
  
  ggboxplot (df.TC.ST_LT, "genotipo", "TCA",facet.by = "duracion"
             , color = "cond.hidr",
             palette = c("#00AFBB", "#E7B800")) 
  

  
  ggdotplot (df.TC.ST_LT, "genotipo", "TCA",
             facet.by = "duracion",
             add = "mean_sd", 
             fill="cond.hidr", color = "cond.hidr",
             palette = c("darkorange", "navyblue"))
  
  ggdotplot (df.TC.ST_LT, "duracion", "TCA",
             facet.by =  "genotipo",
             add = "mean_sd", 
             fill="cond.hidr", color = "cond.hidr",
             palette = c("darkorange", "navyblue"))
  
  
  
  
  ggdotplot (df.TC.ST_LT, "genotipo", "TCA",
             facet.by ="cond.hidr",
             add = "mean_sd", 
             fill="duracion", color = "duracion",
             palette = c("gray48", "red"))
  
  ggdotplot (df.TC.ST_LT, "genotipo", "TCR",
             facet.by ="cond.hidr",
             add = "mean_sd", 
             fill="duracion", color = "duracion",
             palette = c("gray48", "red"))
  
  
  
  
  df.TC.LT <- df.TC %>%
              dplyr::filter (duracion != "ST_LT")
  
  
  ggboxplot (df.TC.LT, "duracion", "TCA",facet.by = "genotipo"
             , color = "cond.hidr",
             palette = c("#00AFBB", "#E7B800")) 
  
  ggboxplot (df.TC.LT, "duracion", "TCR",facet.by = "genotipo"
             , color = "cond.hidr",
             palette = c("#00AFBB", "#E7B800")) 
  
  ggboxplot (df.TC.LT, "genotipo", "TCR",facet.by = "duracion"
             , color = "cond.hidr",
             palette = c("#00AFBB", "#E7B800")) 
  
  ggboxplot (df.TC.LT, "genotipo", "TCA",facet.by = "duracion"
             , color = "cond.hidr",
             palette = c("#00AFBB", "#E7B800")) 
  
  ggdotplot (df.TC.LT, "genotipo", "TCA",
             facet.by = "duracion",
             add = "mean_sd", 
             fill="cond.hidr", color = "cond.hidr",
             palette = c("darkorange", "navyblue")) 
  
  
  ggdotplot (df.TC.LT, "genotipo", "TCA",
             facet.by ="cond.hidr",
             add = "mean_sd", 
             fill="duracion", color = "duracion",
             palette = c("gray48", "red")) 
  
  
  
      dplyr::mutate (dias = Date - date.1) %>%
      dplyr::filter (E > 0)
    
    maceta <-  unique(datos$maceta)
    clon   <-  unique(datos$clon)
    trat   <- unique(datos$tratamiento)
    
    dat.E.1 <- data.frame( maceta=maceta, clon=clon, cond.hidrica = trat, pot=id, dat.E )
    
    dat.E.1 <-  dat.E.1 %>%
      dplyr::select (maceta, clon, cond.hidrica, pot, dias, everything())
    
    
    write.table ( dat.E.1, file =paste("./Data/outpout.E/dt.E","pot_",filtro,".txt",sep=""),
                  append = FALSE, quote = TRUE, sep = ",",
                  eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                  col.names = TRUE)
    
    return (dat.E.1)
    
  }))








t1 = "12/11/2019"
t1.1 <- dmy (t1)

crecimiento.MN.3b <- crecimiento.MN.3b %>%
                      dplyr::mutate (dias = date - t1.1) 


head (crecimiento.MN.3b )

list.pot <- unique (crecimiento.MN.3b$pot)










lapply (list.pot, function (filt.pot){
  
  dt.pot <- crecimiento.MN.3b %>%
            dplyr::filter (pot =="3.g1.d")
  
  id.pot <-  unique (dt.pot$pot)
  
  list.day <- unique (dt.pot$dias )
  
  alt1.1 <- dt.pot %>%
             dplyr::filter (dias == 0 )
  
  alt1.2 <- dt.pot %>%
            dplyr::filter (dias == 21 )
  
  
  
  ggline (dt.pot, x = "dias", y = "diam",
          linetype = "genotipo", title =  id.pot,
          #shape = "genotipo",
          #add = "mean_se",
          color = "genotipo", palette = "black")
  
  
})

                     

lapply (list.pot )
 ggscatter (dt.1 , x = "date", y = "alt", 
            facet.by = "cond.hidr",
            #title = unique(dt.1$pot), 
                         xlab = "time (d)",
                         ylab = "altura (cm)",
                         point=TRUE)
 
dt= crecimiento.MN.3b 
t1 = t1 = "12/11/2019"
 
 
 
 
 head ( )
 
 
 
 ggline (dt.1, x = "dias", y = "diametro",facet.by = "cond.hidr",
         linetype = "genotipo", 
         #shape = "genotipo",
         add = "mean_se",
         color = "genotipo", palette = c("black", "pink","navyblue", "darkorange"))
 
 
 
 
 
 +
  geom_line(color = "gray48", linetype =2, size = 0.5) +
  geom_point(color = "black", size = 1.5) +
  geom_hline(yintercept = mean (dt.1$pesof.1, na.rm = TRUE ), linetype =3, size = 1


### unifico las matrices 

unique (Phi.PSII.MN.3a$ensayo)
unique (Phi.PSII.MN.3a$ambiente)
unique (Phi.PSII.MN.3a$genotipo)
unique (Phi.PSII.MN.3a$pos)
unique (Phi.PSII.MN.3a$protocolo)


# protocolos cambio a actinicas
Phi.PSII.MN.3a <- Phi.PSII.MN.3a  %>%
                  dplyr::mutate (ppfd.act = protocolo)
                       

Phi.PSII.MN.3a$ppfd.act [Phi.PSII.MN.3a$ppfd.act == "nion.1.act25"] <- 1000



Phi.PSII.MN.3a <- Phi.PSII.MN.3a %>%
                  dplyr::mutate (ppfd.amb = 450)

unique (Phi.PSII.MN.3a$date)

Phi.PSII.MN.3b  <- Phi.PSII.MN.3a  %>%
                   dplyr::select (ensayo, ambiente,date, ppfd.amb,protocolo, ppfd.act, genotipo ,cond.hidr, pos,  hoja, everything() )

write_excel_csv (Phi.PSII.MN.3b , file ="./Data/procdata/Phi.PSII.MN.3b.csv",
                  na = "NA",
                  append = FALSE,
                  delim = ";",
                  quote_escape = "double",
                  eol = "\n" )

write_delim (Phi.PSII.MN.3b , file ="./Data/procdata/Phi.PSII.MN.3b.txt", 
             delim = ";", na = "NA")


# saco el ultimo ciclo #### 
unique  (Phi.PSII.MN.3b$protocolo)

nion.1.act25.quenching <- Phi.PSII.MN.3b %>%
                          dplyr::filter (ciclo == max (ciclo) - 1)




Phi.PSII.MN.3.ciclo.q <- bind_rows (nion.1.act25.quenching)


write_excel_csv (Phi.PSII.MN.3.ciclo.q, file ="./Data/procdata/Phi.PSII.MN.3.ciclo.q.csv",
                 na = "NA",
                 append = FALSE,
                 delim = ";",
                 quote_escape = "double",
                 eol = "\n" )

write_delim (Phi.PSII.MN.3.ciclo.q, file ="./Data/procdata/Phi.PSII.MN.3.ciclo.q.txt", 
             delim = ";", na = "NA")

##################################### Plot ###
# filtro los datos fuera de rango

phiNPQslow.neg <- Phi.PSII.MN.3.ciclo.q %>%
                  dplyr::filter (phiNPQslow < 0 ) %>%
                  dplyr::arrange (phiNPQslow )

Phi.PSII.MN.3.ciclo.q.1  <- Phi.PSII.MN.3.ciclo.q %>%
                            #dplyr::filter (phiNPQslow > -0.005 ) %>%
                            dplyr::filter ( qs.fo < 1 )


Phi.PSII.MN.3.ciclo.q.1$phiNPQslow[Phi.PSII.MN.3.ciclo.q.1$phiNPQslow < 0] <- 0


Phi.PSII.MN.3.ciclo.q.1 <- Phi.PSII.MN.3.ciclo.q.1 %>%
                           dplyr::mutate (trat = str_c (ambiente, protocolo, genotipo, cond.hidr)) %>%
                           dplyr::mutate (amb.ppfd = str_c (ambiente, ppfd.amb)) 
  
##### 

# grafico de correlaciones 

phis <- c ("phiPS2", "qp.fo" , "qs.fo" , "phiNPQ", "phiNO", "phiNO.psII" ,
            "phiNO.basal" , "phiNPQfast", "phiNPQslow" ) 


phis.cor <- Phi.PSII.MN.3.ciclo.q.1 %>%
            dplyr::select (all_of(phis))

ggpairs (phis.cor)

corr <- round (cor (phis.cor), 2)

p.mat <- cor_pmat (phis.cor )

ggcorrplot(corr, hc.order = TRUE, type = "lower",
           colors = c("darkred", "white", "blue"),insig ="blank",
           lab = TRUE,p.mat = p.mat, sig.level =0.05)

Phi.PSII.MN.3.ciclo.q.1$date <- as.factor (Phi.PSII.MN.3.ciclo.q.1$date)
Phi.PSII.MN.3.ciclo.q.1$cond.hidr <- as.factor (Phi.PSII.MN.3.ciclo.q.1$cond.hidr)
Phi.PSII.MN.3.ciclo.q.1$genotipo <- as.factor (Phi.PSII.MN.3.ciclo.q.1$genotipo)


#### corremos el modelo 
## model 
lm.phiPS2 <- lm (phiPS2 ~  genotipo + cond.hidr + date + 
                           genotipo * cond.hidr * date
                     , data= Phi.PSII.MN.3.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiPS2$fitted.values, lm.phiPS2$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiPS2")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiPS2$residuals, main="phiPS2",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiPS2$residuals),
              sd (lm.phiPS2$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiPS2$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiPS2")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiPS2$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiPS2 ~ trat, data = Phi.PSII.MN.3.ciclo.q.1)

# ANAVA
(anava.lm.phiPS2  <- anova (lm.phiPS2))

############## qp ######
lm.qp.fo <- lm (qp.fo ~  genotipo + cond.hidr + date + 
                  genotipo * cond.hidr * date
                , data= Phi.PSII.MN.3.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.qp.fo$fitted.values, lm.qp.fo$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="qp.fo")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.qp.fo$residuals, main="qp.fo",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.qp.fo$residuals),
              sd (lm.qp.fo$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.qp.fo$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="qp.fo")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.qp.fo$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( qp.fo ~ trat, data = Phi.PSII.MN.3.ciclo.q.1)

# ANAVA
(anava.lm.qp.fo  <- anova (lm.qp.fo))


####
#Phi.PSII.MN.3.ciclo.q.1$protocolo  <- factor(Phi.PSII.MN.3.ciclo.q.1$protocolo, 
 #                                     levels = c("pz400.2", 
  #                                             "pz400",
   #                                           "pz800.2", "pz800"))

#Phi.PSII.MN.3.ciclo.q.1$ppfd.act  <- factor(Phi.PSII.MN.3.ciclo.q.1$ppfd.act, 
 #                                            levels = c("200", 
  #                                                      "400",
   #                                                     "850", "1700"))
### hacer los sactter

Phi.PSII.MN.3.ciclo.q.1.T0 <- Phi.PSII.MN.3.ciclo.q.1 %>%
                              dplyr::filter ( date == "T0") 


Phi.PSII.MN.3.ciclo.q.1.T1 <- Phi.PSII.MN.3.ciclo.q.1 %>%
                              dplyr::filter ( date == "T1") 


ggscatter (  Phi.PSII.MN.3.ciclo.q.1.T0, x = "qp.fo", y = "phiPS2",
             ylim=c(0,1),
             #xlim=c(0,1),
             facet.by =  "genotipo",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "cond.hidr",
             size =5,
             palette = c("darkorange", "navyblue")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)


ggscatter (  Phi.PSII.MN.3.ciclo.q.1.T1, x = "qp.fo", y = "phiPS2",
             ylim=c(0,1),
             #xlim=c(0,1),
             facet.by =  "genotipo",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "cond.hidr",
             size =5,
             palette = c("darkorange", "navyblue")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)




ggscatter ( Phi.PSII.MN.3.ciclo.q.1.T0, x = "qp.fo", y = "phiPS2",
             ylim=c(0,1),
             #xlim=c(0,1),
             facet.by =  "cond.hidr",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "genotipo",
             size =5,
             palette = c("black","red", "navyblue", "pink")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)


ggscatter ( Phi.PSII.MN.3.ciclo.q.1.T1, x = "qp.fo", y = "phiPS2",
            ylim=c(0,1),
            #xlim=c(0,1),
            facet.by =  "cond.hidr",
            #add="reg.line",
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            color= "genotipo",
            size =5,
            palette = c("black","red", "navyblue", "pink")
            
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)


############## NPQ ######

lm.phiNPQ <- lm (phiNPQ ~  genotipo + cond.hidr + date + 
                   genotipo * cond.hidr * date
                , data= Phi.PSII.MN.3.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiNPQ$fitted.values, lm.phiNPQ$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiNPQ")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiNPQ$residuals, main="phiNPQ",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiNPQ$residuals),
              sd (lm.phiNPQ$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiNPQ$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiNPQ")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiNPQ$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest ( phiNPQ ~ trat, data = Phi.PSII.MN.3.ciclo.q.1)

# ANAVA
(anava.lm.phiNPQ  <- anova (lm.phiNPQ))


############## phiNPQ fast ##################33

lm.phiNPQfast <- lm ( phiNPQfast ~  genotipo + cond.hidr + date + 
                        genotipo * cond.hidr * date
                 , data= Phi.PSII.MN.3.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiNPQfast$fitted.values, lm.phiNPQfast$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiNPQ")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiNPQfast$residuals, main="phiNPQfast",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiNPQfast$residuals),
              sd (lm.phiNPQfast$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiNPQfast$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiNPQ")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiNPQfast$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest (phiNPQfast ~ trat, data = Phi.PSII.MN.3.ciclo.q.1)

# ANAVA
(anava.lm.phiNPQfast  <- anova (lm.phiNPQfast))

  
ggscatter (  Phi.PSII.MN.3.ciclo.q.1.T0, x = "phiNPQfast", y = "phiNPQ",
             ylim=c(0,1),
             #xlim=c(0,1),
             facet.by =  "genotipo",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "cond.hidr",
             size =5,
             palette = c("darkorange", "navyblue")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)


ggscatter (  Phi.PSII.MN.3.ciclo.q.1.T1, x = "phiNPQfast", y = "phiNPQ",
             ylim=c(0,1),
             #xlim=c(0,1),
             facet.by =  "genotipo",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "cond.hidr",
             size =5,
             palette = c("darkorange", "navyblue")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)


### fast vs qp #######3

ggscatter (  Phi.PSII.MN.3.ciclo.q.1.T0, x = "qp.fo", y = "phiNPQfast",
             ylim=c(0,1),
             xlim=c(0,1),
             facet.by =  "genotipo",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "cond.hidr",
             size =5,
             palette = c("darkorange", "navyblue")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)


ggscatter (  Phi.PSII.MN.3.ciclo.q.1.T1, x = "qp.fo", y = "phiNPQfast",
             ylim=c(0,1),
             xlim=c(0,1),
             facet.by =  "genotipo",
             #add="reg.line",
             cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
             color= "cond.hidr",
             size =5,
             palette = c("darkorange", "navyblue")
             
) + geom_vline(xintercept = 0.5, lty =2) +
  geom_hline(yintercept = 0.5, lty =2)


