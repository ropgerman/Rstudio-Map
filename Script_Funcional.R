  library(maptools)
library(ggplot2)
library(scatterpie)
args=commandArgs(TRUE)

#listsnps=args[1]
#listsnps=read.table(listsnps,header=FALSE)
#listsnps_F<-as.character(listsnps$V1)

table_pop<- read.csv("C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\table_pop_AllData__NoMex.txt", sep="\t", header=T)
View(table_pop)

#nuevo_vector<- as.data.frame(table_pop$pop_size)
#View(nuevo_vector)
#View(nuevo_vector)


#CORRER UNA CONDICIONAL A LA VEZ
table_pop$pop_size[table_pop$pop_size < 2] <- 1.5
table_pop$pop_size[(table_pop$pop_size > 2)  & (table_pop$pop_size < 5)] <- 2.5
table_pop$pop_size[(table_pop$pop_size >5) & (table_pop$pop_size < 10)] <- 3.5
table_pop$pop_size[(table_pop$pop_size >10) & (table_pop$pop_size < 15)] <- 4.5
table_pop$pop_size[(table_pop$pop_size >15) & (table_pop$pop_size < 20)] <- 5.5
table_pop$pop_size[(table_pop$pop_size >20) & (table_pop$pop_size < 30)] <- 6.5
table_pop$pop_size[table_pop$pop_size > 30] <- 7.5




setwd("C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\modern_frq")
getwd()


print("##### 1 #######")

rs_List<- read.table("C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\ancient_frq\\List_SNPS_EDAR.txt", header=FALSE)
View(rs_List)

listsnps_F<-as.character(rs_List$V1)
View(listsnps_F)

dat<- read.table('modern_Clusters.dat')
View(dat)

A_dat<- read.table('modern_Clusters__ONLY_America.dat')
View(A_dat)

phenotype<- c(rep("SNPs", 1))  #Lista de fenotipos asociados a casa RS
View(phenotype)

setwd("C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\ancient_frq")
getwd()

rs_ancientList<- listsnps_F
View(rs_ancientList)



print("##### 2 #######")

ancient_dat<- read.table ("C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\ancient_frq\\ancient_Clusters_noMex.dat")
View(ancient_dat)

full_dat<- rbind(dat, ancient_dat)
View(full_dat)

View(rs_List)

print("##### 3 #######")


#rsID1 <-read.table('C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\modern_frq\\rs365060.frq.strat')
#rsID2 <-read.table('C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\modern_frq\\rs3827760.frq.strat')
#rsID <-rbind(rsID1, rsID2)
#View(rsID)

################################################################################


ref<- read.csv("C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\ref.csv", sep=";", header=T)
View(ref)

#par(ps = 20, cex = 1, cex.main = 1)


labelG <- ref$label(sep=";")
View(labelG)



ref1.5 <- ref[1, 2]
View(ref1.5)
ref2.5 <- ref[2, 2]
ref3.5 <- ref[3, 2]
ref4.5 <- ref[4, 2]
ref5.5 <- ref[5, 2]
ref6.5 <- ref[6, 2]
ref7.5 <- ref[7, 2]









# theme_update(text = element_text(size=10))


custom_geom_scatterpie_legend <- 
  function (radius, x, y, n = 7, labeller,textsize=1){
    if (length(radius) > n) {
      radius <- unique(sapply(seq(min(radius), max(radius), 
                                  length.out = n), scatterpie:::round_digit))
    }
    label <- FALSE
    if (!missing(labeller)) {
      if (!inherits(labeller, "function")) {
        stop("labeller should be a function for converting radius")
      }
      label <- TRUE
    }
    dd <- data.frame(r = radius, start = 0, end = 2 * pi, x = x, 
                     y = y + radius - max(radius), maxr = max(radius))
    if (label) {
      dd$label <- labeller(dd$r)
    }
    else {
      dd$label <- dd$r
    }
    list(ggforce:::geom_arc_bar(aes_(x0 = ~x, y0 = ~y, r0 = ~r, r = ~r, 
                                     start = ~start, end = ~end), data = dd, inherit.aes = FALSE), 
         geom_segment(aes_(x = ~x, xend = ~x + maxr * 1.5, y = ~y + 
                             r, yend = ~y + r), data = dd, inherit.aes = FALSE), 
         geom_text(aes_(x = ~x + maxr * 1.6, y = ~y + r, label = ~label), 
                   data = dd, hjust = "left", inherit.aes = FALSE,size=textsize))
  }



par()$lwd

par (lwd = .5)





################################################################################

freq.plot<- function(rsID, phen) {
  rs_ancient<- read.table(paste('rs3827760', ".frq.strat", sep = ""), header = T, colClasses = c("numeric", rep("character", 4), rep("numeric", 3))) # Leer los archivos de Freqs antiguas de plink
  rs<-read.table(paste('C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\modern_frq\\', 'rs3827760', ".frq.strat", sep= ""),  header = T, colClasses = c("numeric", rep("character", 4), rep("numeric", 3))) # Leer los archivos de Freqs
  
  rs <- rbind(rs, rs_ancient) #Unir ambos archivos de freqs
  #radius<- rep(3, nrow(table_pop)) # Radio de los pies
  radius <- rep(table_pop$pop_size)
  rs<- cbind(rs, 1-rs[6], radius)
  names(rs)[9]<- "MajAF"
  MajAF<-rs$MajAF
  MAF<-rs$MAF
  minor_Allele<-rs[1,4]
  major_Allele<-rs[1,5]
  
  # Esta es una función que 'grafica'
  
  
  
  ind<-c()
  
  for (i in rs$CLST) {    # for (i in 1:length(files))
    n_Ind<-length(which(full_dat[,3]==i))
    ind[i]<-(paste(i, 'N=', n_Ind, sep=' '))}
  
  
  table_pop<-cbind(table_pop, radius, MAF, MajAF, ind)
  table_pop<-as.data.frame(table_pop)
  
  pdf('rs3827760','.pdf')
  #pdf(file = "C:\\Users\\ropge\\Documents\\INMEGEN-ADN\\Data_for_Rplots\\modern_frq\\rs365060.pdf")
  
  
  
  
  # Graficar
  worldmap<- map_data ("world")
  mapplot1<- ggplot(worldmap) +
    
    geom_map(data = worldmap, map = worldmap, aes(x = long, y = lat, map_id = region), fill = "gray70") +
    
    geom_scatterpie(aes(x = x_coord, y = y_coord, r=radius), data = table_pop, cols = c("MAF", "MajAF")) +
    
    coord_fixed() + 
    
    custom_geom_scatterpie_legend(1.5, -180, -70, n = 7, labeller= function(a) a = (ref1.5), textsize=1)+
    custom_geom_scatterpie_legend(2.5, -180, -63, n = 7, labeller= function(a) a = (ref2.5), textsize=1)+
  custom_geom_scatterpie_legend(3.5, -180, -55, n = 7, labeller= function(a) a = (ref3.5), textsize=1)+
  custom_geom_scatterpie_legend(4.5, -180, -47, n = 7, labeller= function(a) a = (ref4.5), textsize=1)+
  custom_geom_scatterpie_legend(5.5, -180, -37, n = 7, labeller= function(a) a = (ref5.5), textsize=1)+
  custom_geom_scatterpie_legend(6.5, -180, -25.5, n = 7, labeller= function(a) a = (ref6.5), textsize=1)+
  custom_geom_scatterpie_legend(7.5, -180, -11.5, n = 7, labeller= function(a) a = (ref7.5), textsize=1)+
  
  
    
    
    
    
   
  
  
  scale_fill_manual(values = c("#153E7E", "#6698FF"), labels=c(minor_Allele, major_Allele)) +
    
    labs(title = paste('rs3827760', "(",minor_Allele, (")"), phen, sep=" "), x = "Longitude", y = "Latitude", fill="", subtitle= "Allelic Frequency", caption = 'El tamano de los graficos corresponde al tamano de la muestra para cada poblacion') +
    
    # theme(legend.position = "top", plot.subtitle=element_text(color="#153E7E", face="bold"), plot.title=element_text(face="bold", size=11))+
    
    theme(legend.text=element_text(size=12))
  # theme(legend.text = element_text(colour="black", size=18))
  
  rs_Map<- mapplot1 + geom_text(data=table_pop, aes(label=ifelse( date != 0 , paste(ind, "\n", date, "years"), paste(ind, "")),  x=plot_x, y=plot_y), size=0.87)
  
  

  
  
  print(rs_Map)
  dev.off()
  
  
}

# labs(title = "Diamonds are forever...", 
#     subtitle = "Carat weight by Price", 
#     caption = "H. Wickham. ggplot2: Elegant Graphics for Data Analysis Springer-Verlag New York, 2009."







print('#### 4 ######')

mapply(rs_List[1:length(rs_List)], phenotype,FUN = freq.plot) # Pasar a la función la lista de RSs y el fenotipo de cada uno