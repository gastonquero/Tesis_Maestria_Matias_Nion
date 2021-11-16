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
library(plantecophys)

# Se carga los datos de los cromatogramas

########### HPLC ########## 

pig_HPLC  <- read_delim (file = "./Data/rawdata/eucahplc.txt", 
                                  delim ="\t", quote = "\"",
                                  escape_backslash = FALSE,
                                  escape_double = TRUE, 
                                  col_names = TRUE, 
                                  col_types = NULL,
                                  locale = default_locale(), 
                                  na = "NA")
str(pig_HPLC)
pig_HPLC  <- pig_HPLC %>%
             dplyr::mutate (pot = str_c (maceta,clon,cond.hidrica,momento, sep="_"))%>%
             dplyr::mutate (trat = str_c (clon,cond.hidrica,momento, sep="_"))


  
  
unique (pig_HPLC$Name )

clor.tot  <- c("clorb", "clora","clorar" ) 
#clor.tot  <- c("clorb", "clora") 

carot.tot  <-  c("neo" , "lut", "zea1","zea2", "zea3",
                "carot","viola","antera" ) 

pig_HPLC.clor.tot  <- pig_HPLC %>%
                      dplyr::filter (Name %in% clor.tot) 


pig_HPLC.carot.tot  <- pig_HPLC %>%
                       dplyr::filter (Name %in% carot.tot)


gt <- pig_HPLC.clor.tot %>%
      dplyr::filter (clon =="gt") %>%
      dplyr::filter ( cond.hidrica == "c")



#### sumatoria

pig_HPLC.clor.tot.1 <-  pig_HPLC.clor.tot %>%
                        group_by (pot) %>%
                        summarise (clor_ab.porc = sum (area), clor_ab.pixels = sum (Area)) %>%
                        dplyr::ungroup()



pig_HPLC.carot.tot.1 <-  pig_HPLC.carot.tot %>%
                        group_by (pot) %>%
                        summarise (carotenos.totales.porc = sum (area), carotenos.totales.pixels = sum (Area)) %>%
                        dplyr::ungroup()


pig_HPLC.2 <-  pig_HPLC.clor.tot.1 %>%
               dplyr::inner_join(pig_HPLC.carot.tot.1, by= "pot") %>%
               dplyr::mutate (x= pot)%>%
               tidyr::separate (x, c("maceta","clon","cond.hidrica","momento")) %>%
               dplyr::select (c( maceta,clon,cond.hidrica,momento, pot, everything())) %>%
  dplyr::mutate (trat = str_c (clon,cond.hidrica,momento, sep="_"))




write_delim (pig_HPLC.2, file= "./Data/procdata/pig_HPLC.2.txt",
             delim = ",",na = "NA")

pig_HPLC.2$clon <- factor(pig_HPLC.2$clon , levels=c("g1", "g2", "gc", "gt"))

summary (pig_HPLC.2 )

pig_HPLC.2 <- pig_HPLC.2 %>%
              dplyr::filter ( clor_ab.porc < 55.197 )


ggboxplot (pig_HPLC.2, "clon", "clor_ab.porc",
           facet.by = "momento",
           add = "mean_sd", 
           #fill="cond.hidr", 
           color = "cond.hidrica",
           palette = c("darkcyan", "firebrick1"))

####

lm.clorab <- lm (clor_ab.porc ~  clon + cond.hidrica + momento + 
                   clon * cond.hidrica * momento, data= pig_HPLC.2)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.clorab$fitted.values, lm.clorab$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiNO")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.clorab$residuals, main="clorab",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.clorab$residuals),
              sd (lm.clorab$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.clorab$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="clorab")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.clorab$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest (clor_ab.porc ~ trat, data = pig_HPLC.2)

# ANAVA
(anava.lm.clorab  <- anova (lm.clorab))

Anova_clorab <- anova (lm.clorab)

Anova_clorab <- Anova_clorab %>%
                dplyr::mutate (variable = "clorab") %>%
                dplyr::mutate (FDV = row.names(Anova_clorab)) %>%
               dplyr::select (variable, FDV, everything())




list.geno <- unique (pig_HPLC.2$clon )

run_contraste_pig <- function (data.model = NULL, trait = NULL, subfase =NULL){
  
  dt.c <- bind_rows (lapply( list.geno, function (filt.geno) {
    
    #filt.geno ="g2"
    print (str_c (trait,"_",filt.geno))
    
    em.geno <- emmeans (data.model, ~ cond.hidrica| momento,
                        at = list (clon = filt.geno))
    
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

contrastes_clorab <- run_contraste_pig (data.model = lm.clorab, trait = "clorab", subfase =NULL )

contrastes_clorab  <- contrastes_clorab  %>%
                      dplyr::mutate (emmean = round (emmean, 2 )) %>%
                      dplyr::mutate ( SE = round ( SE, 3 )) %>%
                      dplyr::mutate ( lower.CL  = round ( lower.CL , 2 )) %>%
                      dplyr::mutate ( upper.CL  = round ( upper.CL , 2 ))

write_csv2 (contrastes_clorab , file= "./Data/procdata/contrastes_clorab.csv")


ggboxplot (pig_HPLC.2, "clon", "carotenos.totales.porc",
           facet.by = "momento",
           add = "mean_sd", 
           #fill="cond.hidr", 
           color = "cond.hidrica",
           palette = c("darkcyan", "firebrick1"))

lm.carot.tot <- lm (carotenos.totales.porc ~  clon + cond.hidrica + momento + 
                   clon * cond.hidrica * momento, data= pig_HPLC.2)

# GRAFICO DE RESIDUOS vs. PREDICHOS 
plot (lm.carot.tot$fitted.values, lm.carot.tot$residuals, pch=19, 
      col="darkblue", xlab="Predichos", ylab="Residuos",
      main="phiNO")
abline(h=0,col="red", lwd=2, lty=3)
grid(ny=20,nx=20)

# HISTOGRAMA DE LOS RESIDUALES
hist (lm.carot.tot$residuals, main="carot.tot",freq=F, 
      xlab="Residuos", ylab="Frecuencia",col="gray")
x=seq(-5e-15,9e-15,5e-15)
curve (dnorm (x,mean (lm.carot.tot$residuals),
              sd (lm.carot.tot$residuals)),add=T,lwd=2, col="red", lty=1)

# Q-Q PLOT (requiere libreria)
qqPlot (lm.carot.tot$residuals, pch=19, col="darkblue",
        xlab="Cuantiles teoricos", ylab="Cuantiles observados",
        main="carot.tot")

# PRUEBA DE NORMAILDAD
shapiro.test (lm.carot.tot$residuals)

# HOMOSCEDASTICIDAD: TEST DE LEVENE Y BARTLETT 
leveneTest (clor_ab.porc ~ trat, data = pig_HPLC.2)

# ANAVA
(anava.lm.carot.tot  <- anova (lm.carot.tot))

Anova_carot.tot <- anova (lm.carot.tot)

Anova_carot.tot <- Anova_carot.tot %>%
                   dplyr::mutate (variable = "carot.tot") %>%
                   dplyr::mutate (FDV = row.names(Anova_carot.tot)) %>%
                   dplyr::select (variable, FDV, everything())

contrastes_carot.tot <- run_contraste_pig (data.model = lm.carot.tot, trait = "carot.tot", subfase =NULL )

contrastes_carot.tot  <- contrastes_carot.tot  %>%
                         dplyr::mutate (emmean = round (emmean, 2 )) %>%
                         dplyr::mutate ( SE = round ( SE, 3 )) %>%
                         dplyr::mutate ( lower.CL  = round ( lower.CL , 2 )) %>%
                         dplyr::mutate ( upper.CL  = round ( upper.CL , 2 ))

write_csv2 (contrastes_carot.tot , file= "./Data/procdata/contrastes_carot.tot.csv")

## tabla de ANOVAS

ANOVAS_pig <- bind_rows ( Anova_clorab,
                          Anova_carot.tot)

colnames (ANOVAS_pig )
colnames (ANOVAS_pig) <-c("variable", "FDV",  "DF", "SS", "MS", "F_val" , "Pr(>Fval)")


ANOVAS_pig <- ANOVAS_pig %>%
  dplyr::mutate ( SS = round (SS, 2) ) %>%
  dplyr::mutate ( MS = round (MS, 2) )  %>%
  dplyr::mutate ( F_val = round (F_val, 1))

write_csv2 (ANOVAS_pig, file= "./Data/procdata/ANOVAS_pig.csv")









pig_espectro_HPLC  <- read_delim (file = "./Data/rawdata/Pigmentos.hplc.txt", 
                                  delim ="\t", quote = "\"",
                                  escape_backslash = FALSE,
                                  escape_double = TRUE, 
                                  col_names = TRUE, 
                                  col_types = NULL,
                                  locale = default_locale(), 
                                  na = "NA")

unique (pig_espectro_HPLC$pigment )

pigm <-  c("neo" , "lut", "zea1","zea2", "zea3","clorb",
           "carot", "clora","viola","antera",  "clorar" )    




times.pig <- pig_espectro_HPLC.1 %>%
             group_by(pigment) %>%
             dplyr::summarise( media.t = mean (rtime), sd.time = sd (rtime )) %>%
             dplyr::ungroup() %>%
             dplyr::arrange(media.t)

###

files.cromatogramas <- dir(str_c("./Data/rawdata/cromatogramas"), pattern = "*.txt")

HPLC.crom <- bind_rows (lapply(files.cromatogramas, function (filt.raw) {
  
  dt.crom <-  read.table (file = str_c("./Data/rawdata/cromatogramas/", filt.raw) ,
                         header = FALSE, sep = "\t",dec = ".", skip=27)
                
  dt.crom <- dt.crom %>%
             dplyr::mutate (pot=filt.raw)
  
}))

HPLC.crom.1 <- HPLC.crom %>%
             dplyr::rename (time = V1 ) %>%
             dplyr::rename (Y = V2) %>%
             dplyr::mutate (YmAU = Y/1000) %>%
             dplyr::mutate (x = pot)%>%
             tidyr::separate (x , c("maceta","clon", "cond.hidrica", "momento", NA))%>%
             dplyr::filter (time < 30)

### t2 ##########
t2 <- HPLC.crom.1 %>%
     dplyr::filter (momento == "t2")

list.pot.t2 <- unique ( t2$pot)


plot.crom.2 <- lapply (list.pot.t2, function (filt.pot){
  
  print (filt.pot)
  
  potx <- t2 %>%
          dplyr::filter (pot ==filt.pot)
  
  px <- ggplot (potx, aes(x = time, y = YmAU)) + 
    geom_line (size = 1, linetype =1,color= "gray48") +
    #scale_color_manual(values = c("blue", "darkgreen", "red")) +
    theme_minimal()+
    geom_vline(xintercept = times.pig [[1, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[2, 2]],linetype =3)+ 
    geom_vline(xintercept = times.pig [[3, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[4, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[5, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[6, 2]], linetype =3)+ 
    geom_vline(xintercept = times.pig [[7, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[8, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[9, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[10, 2]], linetype =4)+
    geom_vline(xintercept = times.pig [[11, 2]], linetype =3)+
    ggtitle(filt.pot)
  
  px
  
})

names (plot.crom.2) <- list.pot.t2

## t2
# g1 control 
p16.t2 <- plot.crom.2[["16.g1.c.t2.txt"]]
p11.t2 <- plot.crom.2[["11.g1.c.t2.txt"]]
p6.t2 <- plot.crom.2[["6.g1.c.t2.txt"]]

ggarrange (p16.t2, p11.t2, p6.t2, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

# g1 deficit
p13.t2 <- plot.crom.2[["13.g1.d.t2.txt"]]
p23.t2 <- plot.crom.2[["23.g1.d.t2.txt"]]
p3.t2 <- plot.crom.2[["3.g1.d.t2.txt"]]

ggarrange (p13.t2, p23.t2, p3.t2, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

# g2 control 
p7.t2 <- plot.crom.2[["7.g2.c.t2.txt"]]
p9.t2 <- plot.crom.2[["9.g2.c.t2.txt"]]
p19.t2 <- plot.crom.2[["19.g2.c.t2.txt"]]

ggarrange (p7.t2, p9.t2, p19.t2, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

#g2 deficit
p12.t2 <- plot.crom.2[["12.g2.d.t2.txt"]]
p1.t2 <- plot.crom.2[["1.g2.d.t2.txt"]]
p21.t2 <- plot.crom.2[["21.g2.d.t2.txt"]]

ggarrange (p12.t2, p1.t2, p21.t2, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

### gc control 
p2.t2 <- plot.crom.2[["2.gc.c.t2.txt"]]
p17.t2 <- plot.crom.2[["17.gc.c.t2.txt"]]
p22.t2 <- plot.crom.2[["22.gc.c.t2.txt"]]

ggarrange (p2.t2, p17.t2, p22.t2, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

### gc  deficit 
p5.t2 <- plot.crom.2[["5.gc.d.t2.txt"]]
p8.t2 <- plot.crom.2[["8.gc.d.t2.txt"]]
p15.t2 <- plot.crom.2[["15.gc.d.t2.txt"]]

ggarrange (p15.t2, p8.t2, p5.t2, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

### gt control 
p14.t2 <- plot.crom.2[["14.gt.c.t2.txt"]]
p24.t2 <- plot.crom.2[["24.gt.c.t2.txt"]]
p4.t2 <- plot.crom.2[["4.gt.c.t2.txt"]]

ggarrange (p14.t2, p24.t2, p4.t2, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

### gt  deficit 
p18.t2 <- plot.crom.2[["18.gt.d.t2.txt"]]
p10.t2 <- plot.crom.2[["10.gt.d.t2.txt"]]
p20.t2 <- plot.crom.2[["20.gt.d.t2.txt"]]

ggarrange (p18.t2, p10.t2, p20.t2, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)


####
### t1 ##########
t1 <- HPLC.crom.1 %>%
      dplyr::filter (momento == "t1")

list.pot.t1 <- unique ( t1$pot)


plot.crom.1 <- lapply (list.pot.t1, function (filt.pot){
  
  print (filt.pot)
  
  potx <- t1 %>%
          dplyr::filter (pot ==filt.pot)
  
  px <- ggplot (potx, aes(x = time, y = YmAU)) + 
        geom_line (size = 1, linetype =1,color= "gray48") +
    #scale_color_manual(values = c("blue", "darkgreen", "red")) +
    theme_minimal()+
    geom_vline(xintercept = times.pig [[1, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[2, 2]],linetype =3)+ 
    geom_vline(xintercept = times.pig [[3, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[4, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[5, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[6, 2]], linetype =3)+ 
    geom_vline(xintercept = times.pig [[7, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[8, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[9, 2]], linetype =3)+
    geom_vline(xintercept = times.pig [[10, 2]], linetype =4)+
    geom_vline(xintercept = times.pig [[11, 2]], linetype =3)+
    ggtitle(filt.pot)
  
  px
  
})

names (plot.crom.1) <- list.pot.t1

## t1
# g1 control 
p16.t1 <- plot.crom.1[["16.g1.c.t1.txt"]]
p11.t1 <- plot.crom.1[["11.g1.c.t1.txt"]]
p6.t1 <- plot.crom.1[["6.g1.c.t1.txt"]]

ggarrange (p16.t1, p11.t1, p6.t1, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

# g1 deficit
p13.t1 <- plot.crom.1[["13.g1.d.t1.txt"]]
p23.t1 <- plot.crom.1[["23.g1.d.t1.txt"]]
p3.t1 <- plot.crom.1[["3.g1.d.t1.txt"]]

ggarrange (p13.t1, p23.t1, p3.t1, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

# g2 control 
p7.t1 <- plot.crom.1[["7.g2.c.t1.txt"]]
p9.t1 <- plot.crom.1[["9.g2.c.t1.txt"]]
p19.t1 <- plot.crom.1[["19.g2.c.t1.txt"]]

ggarrange (p7.t1, p9.t1, p19.t1, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

#g2 deficit
p12.t1 <- plot.crom.1[["12.g2.d.t1.txt"]]
p1.t1 <- plot.crom.1[["1.g2.d.t1.txt"]]
p21.t1 <- plot.crom.1[["21.g2.d.t1.txt"]]

ggarrange (p12.t1, p1.t1, p21.t1, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

### gc control 
p2.t1 <- plot.crom.1[["2.gc.c.t1.txt"]]
p17.t1 <- plot.crom.1[["17.gc.c.t1.txt"]]
p22.t1 <- plot.crom.1[["22.gc.c.t1.txt"]]

ggarrange (p2.t1, p17.t1, p22.t1, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

### gc  deficit 
p5.t1 <- plot.crom.1[["5.gc.d.t1.txt"]]
p8.t1 <- plot.crom.1[["8.gc.d.t1.txt"]]
p15.t1 <- plot.crom.1[["15.gc.d.t1.txt"]]

ggarrange (p15.t1, p8.t1, p5.t1, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

### gt control 
p14.t1 <- plot.crom.1[["14.gt.c.t1.txt"]]
p24.t1 <- plot.crom.1[["24.gt.c.t1.txt"]]
p4.t1 <- plot.crom.1[["4.gt.c.t1.txt"]]

ggarrange (p14.t1, p24.t1, p4.t1, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)

### gt  deficit 
p18.t1 <- plot.crom.1[["18.gt.d.t1.txt"]]
p10.t1 <- plot.crom.1[["10.gt.d.t1.txt"]]
p20.t1 <- plot.crom.1[["20.gt.d.t1.txt"]]

ggarrange (p18.t1, p10.t1, p20.t1, ncol = 1, nrow = 3, legend  = "top", common.legend = TRUE)










g2.d.t1 <- HPLC.crom.1 %>%
      dplyr::filter (clon =="g2")%>%
      dplyr::filter (cond.hidrica == "d") %>%
      dplyr::filter (momento == "t1")

g2.c.t1 <- HPLC.crom.1 %>%
  dplyr::filter (clon =="g2")%>%
  dplyr::filter (cond.hidrica == "c") %>%
  dplyr::filter (momento == "t1") 

g2.c.t1.19 <- g2.c.t1 %>%
              dplyr::filter (maceta == 19)

g2.c.t1.7 <- g2.c.t1 %>%
              dplyr::filter (maceta == 7)

g2.c.t1.9 <- g2.c.t1 %>%
             dplyr::filter (maceta == 9)

px <- ggplot (g2.c.t1.9, aes(x = time, y = YmAU)) + 
        geom_line (size = 1, linetype =1,color= "gray48") +
        #scale_color_manual(values = c("blue", "darkgreen", "red")) +
        theme_minimal()+
        geom_vline(xintercept = times.pig [[1, 2]], linetype =3)+
        geom_vline(xintercept = times.pig [[2, 2]],linetype =3)+ 
        geom_vline(xintercept = times.pig [[3, 2]], linetype =3)+
        geom_vline(xintercept = times.pig [[4, 2]], linetype =3)+
        geom_vline(xintercept = times.pig [[5, 2]], linetype =3)+
        geom_vline(xintercept = times.pig [[6, 2]], linetype =3)+ 
        geom_vline(xintercept = times.pig [[7, 2]], linetype =3)+
        geom_vline(xintercept = times.pig [[8, 2]], linetype =3)+
        geom_vline(xintercept = times.pig [[9, 2]], linetype =3)+
        geom_vline(xintercept = times.pig [[10, 2]], linetype =3)






p7 <- ggplot (g2.c.t1.7, aes(x = time, y = YmAU)) + 
  geom_line(aes(color = maceta), size = 1, color= "black") +
  #scale_color_manual(values = c("blue", "darkgreen", "red")) +
  theme_minimal()+
  geom_vline(xintercept = 25)+
  geom_vline(xintercept = 14)


p19 <- ggplot (g2.c.t1.19, aes(x = time, y = YmAU)) + 
  geom_line(aes(color = maceta), size = 1, color= "black") +
  #scale_color_manual(values = c("blue", "darkgreen", "red")) +
  theme_minimal()+
  geom_vline(xintercept = 25)+
  geom_vline(xintercept = 14)

ggarrange (p7, p9,p19, ncol = 1, nrow = 3, common.legend = FALSE)


g2.d.t1.1 <- g2.d.t1 %>%
             dplyr::filter (maceta == 1)

g2.d.t1.12 <- g2.d.t1 %>%
              dplyr::filter (maceta == 12)

g2.d.t1.21 <- g2.d.t1 %>%
              dplyr::filter (maceta == 21)

p1 <- ggplot (g2.d.t1.1, aes(x = time, y = YmAU)) + 
      geom_line( size = 1, color= "black") +
      #scale_color_manual(values = c("blue", "darkgreen", "red")) +
      theme_minimal()+
      geom_vline(xintercept = 25)+
      geom_vline(xintercept = 14) +
      ggtitle("p1")

p12 <- ggplot (g2.d.t1.12, aes(x = time, y = YmAU)) + 
       geom_line( size = 1, color= "black") +
        #scale_color_manual(values = c("blue", "darkgreen", "red")) +
       theme_minimal()+
       geom_vline(xintercept = 25)+
       geom_vline(xintercept = 14)+
       ggtitle("p12")


p21 <- ggplot (g2.d.t1.21, aes(x = time, y = YmAU)) + 
       geom_line(size = 1, color= "black") +
       #scale_color_manual(values = c("blue", "darkgreen", "red")) +
       theme_minimal()+
       geom_vline(xintercept = 25)+
       geom_vline(xintercept = 14)+
       ggtitle("p21")







unique (HPLC.crom.1$pot)


x1 <- ggplot(g2.d.t1, aes(x = time, y = YmAU)) + 
      geom_line(aes(color = maceta), size = 1) +
      scale_color_manual(values = c("blue", "darkgreen", "red")) +
      theme_minimal() +
      geom_vline(xintercept = 25)+
      geom_vline(xintercept = 14)








x2 <- ggplot(g2.c.t1, aes(x = time, y = YmAU)) + 
  geom_line(aes(color = maceta), size = 1) +
  scale_color_manual(values = c("blue", "darkgreen", "red")) +
  theme_minimal()+
  geom_vline(xintercept = 25)+
  geom_vline(xintercept = 14)


ggarrange (x2, x1, ncol = 2, nrow = 1)


pot <- HPLC.crom.1 %>%
       dplyr::filter (pot == "1.g2.d.t1.txt")





svg (filename= str_c ("./Figures/power_dist.espectral/dist.espectral_",id.s,"_", px,".svg"), 
     width=7, 
     height=5, 
     pointsize=12)

limy <- max (pot$YmAU)+ 5
plot (YmAU ~ time,
      ylim = c(0,limy),
      #xlim=  c(dt.lmbd.int$l1,dt.lmbd.int$l2),
      col="black",
      #pch=16,cex=1,
      axes=TRUE,
      #xlab="", 
      #ylab="", 
      type="l",
      #main = str_c (id.s, "power", id.punto ,sep="_",collapse = TRUE),
      data= pot)

lapply (list.band.1, function (filt.band) { 
  dt.lmbd.x <- dt.lmbd %>%
    dplyr::filter (bandwidth == filt.band)
  
  abline ( v=dt.lmbd.x$l1, lty=2, col="gray48" )
  abline ( v=dt.lmbd.x$l2 , lty=2, col="gray48" )
  
  text(x = (dt.lmbd.x$l1 + dt.lmbd.x$l2)/2,
       y= max (dt.1.int.phot.total$uW_nm.cm2 - 10) ,
       label = str_c(dt.lmbd.x$l1 ,"_",dt.lmbd.x$l2),
       cex = 0.9,  col="gray48")
})

dev.off ()

unique (HPLC.crom$pot )




pigmentos_espectro  <- read_delim (file = "./Data/rawdata/pigmentos_espectro.txt", 
                    delim ="\t", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = TRUE, 
                    col_types = NULL,
                    locale = default_locale(), 
                    na = "NA")


pigmentos_espectro <- pigmentos_espectro %>%
                      dplyr::mutate (Cb.prim = replace(Cb, which(Cb < 0), 0)) %>%
                      dplyr::mutate (Cab = Ca/Cb.prim) %>%
                      dplyr::mutate (Cab.plus = Ca+Cb.prim)


head (pigmentos_espectro )



ggdotplot (pigmentos_espectro, "clon", "Ca", 
           fill = "cond.hidric", facet.by = "momento",
           color = "cond.hidric",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pigmentos_espectro, "clon", "Cb.prim", 
           fill = "cond.hidric", facet.by = "momento",
           color = "cond.hidric",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pigmentos_espectro, "clon", "Carot", 
           fill = "cond.hidric", facet.by = "momento",
           color = "cond.hidric",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pigmentos_espectro, "clon", "Cab", 
           fill = "cond.hidric", facet.by = "momento",
           color = "cond.hidric",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (pigmentos_espectro, "clon", "Cab.plus", 
           fill = "cond.hidric", facet.by = "momento",
           color = "cond.hidric",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


pigmentos_espectro.T0 <- pigmentos_espectro %>%
                         dplyr::filter (momento == "t1") %>%
                         dplyr::filter (Cb.prim < 8)  


ggdotplot (pigmentos_espectro.T0, "clon", "Cb.prim", 
           fill = "cond.hidric", facet.by = "momento",
           color = "cond.hidric",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


pigmentos_espectro.T1 <- pigmentos_espectro %>%
                         dplyr::filter (momento == "t2") 


ggdotplot (pigmentos_espectro.T1, "clon", "Cb.prim", 
           fill = "cond.hidric", facet.by = "momento",
           color = "cond.hidric",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

########### HPLC ########## 

pig_espectro_HPLC  <- read_delim (file = "./Data/rawdata/pigmentos_HPLC.txt", 
                                   delim ="\t", quote = "\"",
                                   escape_backslash = FALSE,
                                   escape_double = TRUE, 
                                   col_names = TRUE, 
                                   col_types = NULL,
                                   locale = default_locale(), 
                                   na = "NA")
carotenos <- c ( "lut" , "neo" , "zea1" ,
                  "zea2" , "zea3" , "antera" , "carot" , "viola")

pig_espectro_HPLC.carot  <- pig_espectro_HPLC %>%
                            dplyr::filter (pigmento %in% carotenos)


x  <-  pig_espectro_HPLC.carot %>%
                             group_by (maceta, clon,  con.hidrica, momento )%>%
                             summarise(sum.carot = sum (area , na.rm = TRUE))

ggdotplot (x , "clon", "sum.carot", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))




pig_espectro_HPLC.W  <- pig_espectro_HPLC %>%
                        pivot_wider(names_from = pigmento, values_from = area)


pig_espectro_HPLC.W.1 <- pig_espectro_HPLC.W %>%
                         dplyr::mutate (car.totales = lut + neo + zea1 +
                                          zea2 + zea3 + antera + carot + viola) %>%
                         dplyr::mutate (clox= ((cloro.a +cloro.b)*100)/ car.totales )

ggdotplot (pig_espectro_HPLC.W.1, "clon", "car.totales", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (pig_espectro_HPLC.W, "clon", "viola", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggscatter (pig_espectro_HPLC.W, "zea3", "viola", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (pig_espectro_HPLC.W, "clon", "cloro.b", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))


ggdotplot (pig_espectro_HPLC.W, "clon", "clor.a.rar", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pig_espectro_HPLC.W, "clon", "lut", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))



ggdotplot (pig_espectro_HPLC.W, "clon", "neo", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pig_espectro_HPLC.W, "clon", "zea1", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pig_espectro_HPLC.W, "clon", "zea2", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pig_espectro_HPLC.W, "clon", "zea3", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pig_espectro_HPLC.W, "clon", "antera", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))

ggdotplot (pig_espectro_HPLC.W, "clon", "carot", 
           fill = "con.hidrica", facet.by = "momento",
           color = "con.hidrica",
           #ylim=c(min(fit.An.MM.3$Photo)-5,max(fit.An.MM.3$Photo )+5),
           add = "mean_sd",
           palette = c("navyblue", "darkorange"))



antera.rar <dbl>, cloro.b <dbl>,
#   lut <dbl>, neo <dbl>, zea1 <dbl>,
#   zea2 <dbl>, zea3 <dbl>, antera <dbl>,
#   carot <dbl>, clor.a.rar <dbl>,
#   cloro.a <dbl>, viola <dbl>

head (pig_espectro_HPLC.1)


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











