######################################################################################
# Codigo para el analisis de los Datos de fluorescencia de la Tesis de Maestria     #
#  Matias Nion                                                                      #
#                                                                                    #
# Gaston Quero - Matias Nion - Sebastian Simondi                                                   #          
# 19/7/2021                                                                          #
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

dir.data.base <- ("C:/Users/Usuario/OneDrive/Documentos/Paper_fotosintesis_dinamica")


#C:\Users\Usuario\OneDrive\Documentos\Paper_fotosintesis_dinamica\Data\procdata\MN_3


## MN.3
Phi.PSII.MN.3 <- read_delim (file =str_c(dir.data.base,"/Data/procdata/MN_3/Phi.PSII.MN.3.txt"), 
                             delim =";", quote = "\"", escape_backslash = FALSE,
                             escape_double = TRUE, col_names = TRUE, col_types = NULL,
                             na = "NA")

Phi.PSII.MN.3$id

Phi.PSII.MN.3a <- Phi.PSII.MN.3 %>%
                  tidyr::separate (id , c("ensayo","ambiente","genotipo","cond.hidr", "pos", "hoja","protocolo", "date"), sep="_") 




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
                  escape = "double",
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
                 escape = "double",
                 eol = "\n" )

write_delim (Phi.PSII.MN.3.ciclo.q, file ="./Data/procdata/Phi.PSII.MN.3.ciclo.q.txt", 
             delim = ";", na = "NA")

##################################### Plot ###
# filtro los datos fuera de rango
summary (Phi.PSII.MN.3.ciclo.q )

Phi.PSII.MN.3.ciclo.q["phiNPQslow"][Phi.PSII.MN.3.ciclo.q["phiNPQslow"] < 0] <- 0
dim ( Phi.PSII.MN.3.ciclo.q)

phiNPQslow.neg <- Phi.PSII.MN.3.ciclo.q %>%
                  dplyr::filter (phiNPQslow < 0 ) %>%
                  dplyr::arrange (phiNPQslow )

summary (phiNPQslow.neg )

Phi.PSII.MN.3.ciclo.q.1  <- Phi.PSII.MN.3.ciclo.q %>%
                            #dplyr::filter (phiNPQslow > -0.005 ) %>%
                            dplyr::filter ( qs.fo < 1 )
dim ( Phi.PSII.MN.3.ciclo.q.1)

summary (Phi.PSII.MN.3.ciclo.q.1 )


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


data.model = lm.phiPS2

#trait = "Phips2"
#subfase =NULL

list.geno <- unique(Phi.PSII.MN.3.ciclo.q.1$genotipo )

run_contraste_phi <- function (data.model = NULL, trait = NULL, subfase =NULL){
  
  dt.c <- bind_rows (lapply( list.geno, function (filt.geno) {
    
    #filt.geno ="g2"
    print (str_c (trait,"_",filt.geno))
    
    em.geno <- emmeans (data.model, ~cond.hidr |date,
                        at = list (genotipo = filt.geno))
    
    cr.geno <- contrast (em.geno ,
                         method = "pairwise")
    
    cr.geno.1 <- as.data.frame (cr.geno )
    
    contrastes <- tibble( genotipo = filt.geno,
                          subfase=subfase, trait =trait, cr.geno.1)
    
    write_csv2 (contrastes, file= str_c ("./Data/procdata/contrastes","_",subfase,"_", trait,"_", filt.geno, ".csv"))
    
    em.geno <- cbind (genotipo = filt.geno, trait =trait,
                      cld (em.geno, sort=FALSE))
    
    
  }))
  
  
}

contrastes_Phi.PS2 <- run_contraste_phi (data.model = lm.phiPS2, trait = "Phips2", subfase =NULL )


contrastes_Phi.PS2 <- contrastes_Phi.PS2 %>%
                      dplyr::mutate (emmean = round (emmean, 2 )) %>%
                      dplyr::mutate ( SE = round ( SE, 3 )) %>%
                      dplyr::mutate ( lower.CL  = round ( lower.CL , 2 )) %>%
                      dplyr::mutate ( upper.CL  = round ( upper.CL , 2 ))

write_csv2 (contrastes_Phi.PS2 , file= "./Data/procdata/contrastes_Phi.PS2.csv")







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


contrastes_phiNPQ <- run_contraste_phi (data.model = lm.phiNPQ, trait = "phiNPQ", subfase =NULL )

contrastes_phiNPQ  <- contrastes_phiNPQ  %>%
  dplyr::mutate (emmean = round (emmean, 2 )) %>%
  dplyr::mutate ( SE = round ( SE, 3 )) %>%
  dplyr::mutate ( lower.CL  = round ( lower.CL , 2 )) %>%
  dplyr::mutate ( upper.CL  = round ( upper.CL , 2 ))

write_csv2 (contrastes_phiNPQ , file= "./Data/procdata/contrastes_phiNPQ.csv")




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

  
head (Phi.PSII.MN.3.ciclo.q.1 )
############## phiNO ##################33

lm.phiNO <- lm (phiNO ~  genotipo + cond.hidr + date + 
                        genotipo * cond.hidr * date
                      , data= Phi.PSII.MN.3.ciclo.q.1)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.phiNO$fitted.values, lm.phiNO$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiNO")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.phiNO$residuals, main="phiNO",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.phiNO$residuals),
              sd (lm.phiNO$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.phiNO$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="phiNO")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.phiNO$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest (phiNO ~ trat, data = Phi.PSII.MN.3.ciclo.q.1)

# ANAVA
(anava.lm.phiNO  <- anova (lm.phiNO))


contrastes_phiNO <- run_contraste_phi (data.model = lm.phiNO, trait = "phiNO", subfase =NULL )

contrastes_phiNO  <- contrastes_phiNO  %>%
  dplyr::mutate (emmean = round (emmean, 2 )) %>%
  dplyr::mutate ( SE = round ( SE, 3 )) %>%
  dplyr::mutate ( lower.CL  = round ( lower.CL , 2 )) %>%
  dplyr::mutate ( upper.CL  = round ( upper.CL , 2 ))

write_csv2 (contrastes_phiNO , file= "./Data/procdata/contrastes_phiNO.csv")









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


