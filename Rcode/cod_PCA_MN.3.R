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
library("Hmisc")
library("PerformanceAnalytics")
library (GGally)
library(ggcorrplot)
library (readr)
library (igraph)
library (ggraph)
library(tidygraph)

#install.packages("PerformanceAnalytics")
### seteo del directorio donde estan los datos

dir.data.base <- ("C:/Users/Usuario/OneDrive/Documentos/Paper_fotosintesis_dinamica")

## MN.3
Phi.PSII.MN.3 <- read_delim (file =str_c(dir.data.base,"/Data/procdata/MN_3/Phi.PSII.MN.3.txt"), 
                             delim =";", quote = "\"", escape_backslash = FALSE,
                             escape_double = TRUE, col_names = TRUE, col_types = NULL,
                             na = "NA")


Phi.PSII.MN.3a <- Phi.PSII.MN.3 %>%
                  tidyr::separate (id , c("ensayo","ambiente","genotipo","cond.hidr", "pos", "hoja","protocolo", "date"), sep="_") 




### unifico las matrices 

unique (Phi.PSII.MN.3a$ensayo)
unique (Phi.PSII.MN.3a$ambiente)
unique (Phi.PSII.MN.3a$genotipo)
unique (Phi.PSII.MN.3a$pos)
unique (Phi.PSII.MN.3a$protocolo)

unique (Phi.PSII.MN.3a$date)

# protocolos cambio a actinicas
Phi.PSII.MN.3a <- Phi.PSII.MN.3a  %>%
                  dplyr::mutate (ppfd.act = protocolo)%>%
                  dplyr::mutate (ppfd.amb = 450)

Phi.PSII.MN.3a$ppfd.act [Phi.PSII.MN.3a$ppfd.act == "nion.1.act25"] <- 1000


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


# reingreso los datos y saco el ultimo ciclo #### 

Phi.PSII.MN.3b <- read_delim (file ="./Data/procdata/Phi.PSII.MN.3b.txt", 
                  delim =";", quote = "\"", escape_backslash = FALSE,
                  escape_double = TRUE, col_names = TRUE, col_types = NULL,
                  na = "NA")

unique  (Phi.PSII.MN.3b$protocolo)

nion.1.act25.quenching <- Phi.PSII.MN.3b %>%
                          dplyr::filter (ciclo == max (ciclo) - 1)

Phi.PSII.MN.3.ciclo.q <- nion.1.act25.quenching


write_excel_csv (Phi.PSII.MN.3.ciclo.q, file ="./Data/procdata/Phi.PSII.MN.3.ciclo.q.csv",
                 na = "NA", append = FALSE, delim = ";", escape = "double",eol = "\n" )

write_delim (Phi.PSII.MN.3.ciclo.q, file ="./Data/procdata/Phi.PSII.MN.3.ciclo.q.txt", 
             delim = ";", na = "NA")


##################################### Plot ###

Phi.PSII.MN.3.ciclo.q <- read_delim (file ="./Data/procdata/Phi.PSII.MN.3.ciclo.q.txt", 
                              delim =";", quote = "\"", escape_backslash = FALSE,
                              escape_double = TRUE, col_names = TRUE, col_types = NULL,
                              na = "NA")



# filtro los datos fuera de rango

phiNPQslow.neg <- Phi.PSII.MN.3.ciclo.q %>%
                  dplyr::filter (phiNPQslow < 0 ) %>%
                  dplyr::arrange (phiNPQslow )

Phi.PSII.MN.3.ciclo.q.1  <- Phi.PSII.MN.3.ciclo.q %>%
                            dplyr::filter (phiNPQslow > -0.009 ) %>%
                            dplyr::filter ( qs.fo < 1 )


Phi.PSII.MN.3.ciclo.q.1$phiNPQslow[Phi.PSII.MN.3.ciclo.q.1$phiNPQslow < 0] <- 0


Phi.PSII.MN.3.ciclo.q.1 <- Phi.PSII.MN.3.ciclo.q.1 %>%
                           dplyr::mutate (trat = str_c (ambiente, protocolo, genotipo)) %>%
                           dplyr::mutate (amb.ppfd = str_c (ambiente, ppfd.amb)) 
  
unique ( Phi.PSII.MN.3.ciclo.q.1$amb.ppfd  )

Phi.PSII.MN.3.ciclo.q.1 <- Phi.PSII.MN.3.ciclo.q.1  %>%
                           dplyr::mutate (vr = str_c (genotipo,cond.hidr   ))

unique (Phi.PSII.MN.3.ciclo.q.1$vr)

# primer grafo
str(Phi.PSII.MN.3.ciclo.q.1 )

# PCA con la matrices del 
dt.PCA.MN.3.q1 <- Phi.PSII.MN.3.ciclo.q.1 %>%
              dplyr::select(-c (ensayo, ambiente, ppfd.amb,
                                protocolo, ppfd.act, pos, hoja,ciclo,
                                filt.pot,sphi,time.act.On,  time.act.Off,
                                sphiNPQ,trat,amb.ppfd,vr))

MN.3.pca.1 <- PCA ( dt.PCA.MN.3.q1, scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                    quanti.sup = NULL, quali.sup =c(1,2,3), row.w = NULL, 
                    col.w = NULL, graph = TRUE)

plot (MN.3.pca.1, choix = c("ind"),invisible ="ind")
plot (MN.3.pca.1, choix = c("ind"))

control <- dt.PCA.MN.3.q1 %>%
           dplyr::filter (cond.hidr == "c") %>%
           dplyr::select (-cond.hidr)


MN.3.pca.1c <- PCA ( control , scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                    quanti.sup = NULL, quali.sup =c(1,2), row.w = NULL, 
                    col.w = NULL, graph = TRUE)

plot (MN.3.pca.1c, choix = c("ind"),invisible ="ind")


deficit <- dt.PCA.MN.3.q1 %>%
  dplyr::filter (cond.hidr == "d") %>%
  dplyr::select (-cond.hidr)


MN.3.pca.1d<- PCA ( deficit , scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                     quanti.sup = NULL, quali.sup =c(1,2), row.w = NULL, 
                     col.w = NULL, graph = TRUE)

plot (MN.3.pca.1d, choix = c("ind"),invisible ="ind")


T0 <- dt.PCA.MN.3.q1 %>%
      dplyr::filter (date == "T0") %>%
       dplyr::select (-date)


MN.3.pca.1T0 <- PCA ( T0 , scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                    quanti.sup = NULL, quali.sup =c(1,2), row.w = NULL, 
                    col.w = NULL, graph = TRUE)

plot (MN.3.pca.1T0, choix = c("ind"),invisible ="ind")


T1 <- dt.PCA.MN.3.q1 %>%
  dplyr::filter (date == "T1") %>%
  dplyr::select (-date)


MN.3.pca.1T1 <- PCA ( T1 , scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                      quanti.sup = NULL, quali.sup =c(1,2), row.w = NULL, 
                      col.w = NULL, graph = TRUE)

plot (MN.3.pca.1T1, choix = c("ind"),invisible ="ind")



MN.3.pca.1 <- PCA ( dt.PCA.MN.3.q1, scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                        quanti.sup = NULL, quali.sup =c(1,2,3), row.w = NULL, 
                        col.w = NULL, graph = FALSE)


MN.3.pca.1a <- PCA ( dt.PCA.MN.3.q1, scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                    quanti.sup = NULL, quali.sup =c(1,2,3), row.w = NULL, 
                    col.w = NULL, graph = FALSE)






plot (MN.3.pca.1, choix = c("var"), title ="PCA_MN3")

fviz_pca_var(MN.3.pca.1 , col.var = "black" )

fviz_contrib(MN.3.pca.1, choice="var", axes=1, top= 10 )
fviz_contrib(MN.3.pca.1, choice="var", axes=2, top= 10 )

fviz_contrib(MN.3.pca.1, choice="var", axes=1:2, top= 10 )

fviz_pca_var (MN.3.pca.1, col.var ="contrib", 
              gradient.cols =c("navyblue", "orange", "red") )


# Create a random continuous variable of length 10 
set.seed(123) 
my.cont.var <- rnorm(9)


# Color variables by the continuous variable 
fviz_pca_var (MN.3.pca.1, 
              col.var = my.cont.var, gradient.cols = c("blue", "yellow", "red"), 
              legend.title = "Cont.Var")

set.seed(123) 
MN.3.pca.1$var


res.km <- kmeans (MN.3.pca.1$var$coord, centers = 3, nstart = 25) 
grp <- as.factor(res.km$cluster)
# Color variables by groups f
fviz_pca_var(MN.3.pca.1, 
             col.var = grp, 
             palette = c("blue", "red", "darkgreen"), legend.title = "Cluster")


# The variable Species (index = 5) is removed  before PCA analysis 

MN.3.pca.1 <- PCA ( dt.PCA.MN.3.q1, scale.unit = TRUE, ncp = 2, ind.sup = NULL, 
                    quanti.sup = NULL, quali.sup =c(1,2,3), row.w = NULL, 
                    col.w = NULL, graph = FALSE)


geno.pca <- PCA(dt.PCA.MN.3.q1 [,-2], 
                quali.sup =c(1,2),
                graph = FALSE)




#In the R code below: the argument habillage or col.ind can be used to specify the factor variable for coloring the individuals by groups.
#To add a concentration ellipse around each group, specify the argument addEllipses = TRUE. The argument palette can be used to change group colors.

fviz_pca_ind (geno.pca, geom.ind = "point", # show points only (nbut not "text") 
             col.ind = dt.PCA.MN.3.q1$genotipo, # color by groups 
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "black"), addEllipses = TRUE, # Concentration ellipses 
             legend.title = "Groups" )
             Individuals



  fviz_pca_biplot (MN.3.pca.1, label = "var",
                  habillage=1,
                  addEllipses=TRUE, ellipse.level=0.95,
                  ggtheme = theme_minimal())      






summary ( Phi.PSII.MN.3.ciclo.q.1 )

unique (Phi.PSII.MN.3.ciclo.q.1$genotipo )
unique (Phi.PSII.MN.3.ciclo.q.1$date)
unique (Phi.PSII.MN.3.ciclo.q.1$cond.hidr )
## 
#G1.T0.c <- Phi.PSII.MN.3.ciclo.q.1 %>%
 #         dplyr::filter (genotipo == "g1")%>%
  #        dplyr::filter (date == "T0") %>%
   #       dplyr::filter (cond.hidr ==  "c")

#G1.T0.d <- Phi.PSII.MN.3.ciclo.q.1 %>%
 #          dplyr::filter (genotipo == "g1")%>%
  #         dplyr::filter (date == "T0") %>%
   #        dplyr::filter (cond.hidr ==  "d")




run.plot.graphos <- function ( dt= NULL ){
  
list.x <- unique (dt$vr)

lapply (list.x, function (filt.x) {
  dtG <- dt %>%
         dplyr::filter ( vr == filt.x)
  
  phis.2 <- c ("qp.fo" , "qs.fo" , "phiNO.psII" ,
               "phiNO.basal" , "phiNPQfast", "phiNPQslow" ) 
  
  phis.cor <-  dtG %>%
    dplyr::select (all_of(phis.2))
  
  
  # hago las correlaciones 
  res.phi.2 <- rcorr (as.matrix (phis.cor))
  #res.phi.2 <- cor (as.matrix (phis.cor))
  #res.phi.2 
  # esta es la funcion que transforma la matriz de correlaciones 
  # en un data.frame   
  
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  
  # genero el data frame con las correlaciones
  res.phi.3 <- flattenCorrMatrix (res.phi.2$r, res.phi.2$P) ### esto lo puedo usar para el grafo
  
  #Create nodes list
  res.phi.3 <- as_tibble (res.phi.3)
  
  #  Get distinct source names
  origen <- res.phi.3 %>%
    dplyr::distinct (row) %>%
    dplyr::rename(label = row)
  
  # Get distinct destination names
  destino <- res.phi.3 %>%
    dplyr::distinct (column) %>%
    dplyr::rename (label = column)
  
  vertices <- full_join (origen, destino, by = "label") 
  
  vertices <- vertices %>%
    dplyr::mutate(id = 1:nrow(vertices)) %>%
    select(id, everything())
  
  # Rename the n.call column to weight
  res.phi.3 <- res.phi.3 %>%
               rename(weight = cor)
  
  # (a) Join nodes id for source column
  aristas <- res.phi.3 %>% 
    left_join (vertices, by = c("row" = "label")) %>% 
    dplyr::rename (from = id) 
  
  # (b) Join nodes id for destination column
  #fujoE <- tibble (energia = c("F.C", "F.Cb", "C.Cb" ) )
  
  aristas <- aristas %>% 
    left_join (vertices, by = c("column" = "label")) %>% 
    dplyr::rename(to = id) %>%
    dplyr::mutate (weight = round (weight, 2))%>%
    dplyr::select (from, to, weight)
  
  set.seed(123)

  # tidygraph and ggraph
  
  #Create a network object using tidygraph:
  #  Key function: tbl_graph().
  # key arguments: nodes, edges and directed.

  
  net.tidy.principal.phi <- tbl_graph (
    nodes = vertices, edges = aristas, directed = FALSE
  )
  
  # Not specifying the layout - defaults to "auto"
  
  layout.x <- create_layout (net.tidy.principal.phi, layout = 'auto')
  
  layout.x$x[1] <- 2
  layout.x$x[2] <- 4
  layout.x$x[3] <- 1
  layout.x$x[4] <- 5
  layout.x$x[5] <- 2
  layout.x$x[6] <- 4
  
  layout.x$y[1] <- 3
  layout.x$y[2] <- 3
  layout.x$y[3] <- 2
  layout.x$y[4] <- 2
  layout.x$y[5] <- 1
  layout.x$y[6] <- 1
  
  id <- unique (dtG$vr)
  
  ggraph (layout.x) + 
    geom_edge_link ( alpha = 0.25, 
                     aes(width = abs( weight), label = weight, repel = TRUE)) +
    geom_node_point(aes ( size =3, alpha =0.5)) +
    geom_node_text(aes(label = label), repel = TRUE) +
    theme_graph() +
    ggtitle (id)

})
}

run.plot.graphos ( dt= Phi.PSII.MN.3.ciclo.q.1 )

run.plot.graphos ( dt= G2C5 )





# selecciono las variables a correlacionar
phis.cor. <- Phi.PSII.MN.3.ciclo.q.1 %>%
            dplyr::select (all_of(phis.2))

# hago las correlaciones 
res.phi.2 <- rcorr (as.matrix (phis.cor))
res.phi.2 

# esta es la funcion que transforma la matriz de correlaciones en un data.frame
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# genero el data frame con las correlaciones
res.phi.3 <- flattenCorrMatrix (res.phi.2$r, res.phi.2$P) ### esto lo puedo usar para el grafo


# To visualize the network graph, we need to create two data frames from the data sets:

# nodes list: containing nodes labels and other nodes attributes
# edges list: containing the relationship between the nodes. It consists of the edge list and any additional edge attributes.


#Create nodes list

#Create nodes list

res.phi.3 <- as_tibble (res.phi.3)
#  Get distinct source names
origen <- res.phi.3 %>%
  dplyr::distinct (row) %>%
  dplyr::rename(label = row)

# Get distinct destination names
destino <- res.phi.3 %>%
  dplyr::distinct (column) %>%
  dplyr::rename (label = column)


# Join the two data to create node
# Add unique ID for each country
vertices <- full_join (origen, destino, by = "label") 

vertices <- vertices %>%
  dplyr::mutate(id = 1:nrow(vertices)) %>%
  select(id, everything())

head (vertices , 3)

#Take the phone.call data, which are already in edges list format,
# showing the connection between nodes. 
# Rename the column “n.call” to “weight”.
# Join the node IDs to the edges list data
# Do this for the “source” column and rename the id column that are brought over from nodes. New name: “from”.
# Do this for the “destination” column and rename the id column. New name: “to”
# Select only the columns “from” and “to” in the edge data. 
# We don’t need to keep the column “source” and “destination” containing the names of countries. 
# These information are already present in the node data.

# Rename the n.call column to weight
res.phi.3 <- res.phi.3 %>%
             rename(weight = cor)

# (a) Join nodes id for source column
aristas <- res.phi.3 %>% 
           left_join (vertices, by = c("row" = "label")) %>% 
           dplyr::rename (from = id) 

# (b) Join nodes id for destination column
#fujoE <- tibble (energia = c("F.C", "F.Cb", "C.Cb" ) )

aristas <- aristas %>% 
           left_join (vertices, by = c("column" = "label")) %>% 
           dplyr::rename(to = id) %>%
           dplyr::mutate (weight = round (weight, 2))%>%
           dplyr::select (from, to, weight)




head(aristas, 3)

set.seed(123)

#### igraph #######
## Create an igraph network object:
# Key R function: graph_from_data_frame().

# Key arguments:
# d: edge list
# vertices: node list
# directed: can be either TRUE or FALSE depending on whether the data is directed or undirected.

net.principal.phi <- graph_from_data_frame(
  d = aristas, vertices = vertices, 
  directed = FALSE
)

# Create a network graph with igraph

set.seed(321)
plot (net.principal.phi, 
      edge.arrow.size = 0.2,
      layout = layout_with_graphopt)



# tidygraph and ggraph

#Create a network object using tidygraph:
#  Key function: tbl_graph().
# key arguments: nodes, edges and directed.
#library(tidygraph)

net.tidy.principal.phi <- tbl_graph (
  nodes = vertices, edges = aristas, directed = FALSE
)

# Not specifying the layout - defaults to "auto"



layout.x <- create_layout (net.tidy.principal.phi, layout = 'auto')

layout.x$x[1] <- 2
layout.x$x[2] <- 4
layout.x$x[3] <- 1
layout.x$x[4] <- 5
layout.x$x[5] <- 2
layout.x$x[6] <- 4

layout.x$y[1] <- 3
layout.x$y[2] <- 3
layout.x$y[3] <- 2
layout.x$y[4] <- 2
layout.x$y[5] <- 1
layout.x$y[6] <- 1

ggraph (layout.x) + 
  geom_edge_link ( alpha = 0.25, 
                   aes(width = abs( weight), label = weight, repel = TRUE)) +
  geom_node_point(aes ( size =3, alpha =0.5)) +
  geom_node_text(aes(label = label), repel = TRUE) +
  theme_graph() 

############## aca termina el cdogp para los grafos ################
