rm(list=ls())

# setwd("D:/Dropbox/Artigos em andamento/FoodWeb")

#I. Load Packages and data preparation-------------------------------------------------------------
###1. Packages and functions----
# List the packages
used.packages <- c("tidyverse", # a collection of packages, including ggplot, dplyr and more
                   # "plyr",
                   "reshape2", # data transformation
                   "ggrepel", # avoid overlap on figures
                   "vegan", # spatial dissimilarity
                   "fluxweb", #energy flux calculation
                   "ggplot2",
                   # "dplyr",
                   "maps",
                   "ggpubr",
                   "mice") #data imputation

# Load all packages in the working space
lapply(used.packages, library, character.only = TRUE)

# Load the functions made for FaunaServices analyses
source('wilian/FaunaServices_functions.R')

###2. Load data ----
#download the fauna dataset
Raw_Dataset <- read.table("wilian/faunaServices_transects.txt", quote = "", sep = "\t", h = T, na = c("na","NA","nd"))

#download additional landuse
landuse<-read.csv("wilian/faunaServices_transects_landuse.csv",sep="\t")
landuse<-subset(landuse,select=c("transect_id","lu_cover1","lu_cover2"))
Raw_Dataset<-merge(Raw_Dataset,landuse,by="transect_id")
Raw_Dataset$lu_cover1[which(Raw_Dataset$lu_cover1=="na")]<-NA
Raw_Dataset$lu_cover2[which(Raw_Dataset$lu_cover2=="")]<-NA

###3. Prepare data for analyses--------
####Select Tr with core groups------
taxa_selection_included <- c("hymenoptera_total","arachnida_total","chilopoda","coleoptera_total","diplopoda",
                             "diptera_total","earthworms","hemiptera_total","isopoda","termites",
                             "gastropoda_total","cockroaches","dermaptera","lepidoptera_total",
                             "orthoptera")

#Run the code for data preparate to clean data and select transects with all taxa from selection
Data_Transect_selection<-dataprepare.coregroups(Raw_Dataset,taxa_selection_included)

###Select max and min

Data_Transect_selection <- Data_Transect_selection %>%
  filter(biomass=="yes" &
           latitude.y >= -23.7 & latitude.y <= -22.67 &
           longitude.y >= -52.54 & longitude.y <= -50.83)

dim(Data_Transect_selection)


#########Ignore it########
#Check this worked ok
Check_Transect <- Data_Transect_selection %>%
  group_by(transect_id,biome,lu_cover1) %>%
  dplyr::summarise(Key_taxa = all(taxa_selection_included %in% Group))
table(Check_Transect$lu_cover1)
table(Check_Transect$biome)
table(Check_Transect$Key_taxa)

####Select Tr with biomass of core groups----
#Selection of transects with biomass data for all groups
Data_Transect_BodyMass <- subset(Data_Transect_selection, !is.na(Biomass_fresh_g_sqm))

####Ignore it####
Check_Transect <- Data_Transect_BodyMass %>%
  group_by(transect_id,biome,lu_cover1) %>%
  dplyr::summarise(Key_taxa = all(taxa_selection_included %in% Group))

Check_Transect_selection_biomass <- subset(Check_Transect,Key_taxa==T)
table(Check_Transect_selection_biomass$biome)
table(Check_Transect$lu_cover1)

#Data_Transect_BodyMass <- subset(Data_Transect_selection, transect_id %in% Check_Transect_selection_biomass$transect_id)

####Select all data with body mass of core groups------
Data_BodyMass_selectedGroup<-dataprepare.bodymasscore(Raw_Dataset,taxa_selection_included)

#II. Basic data visualization-------------------------------------------------------------

###1. Geographic location---------------
world_map <- map_data("world")
sa.countries <- c('Argentina', 'Bolivia', 'Brazil',
                  'Colombia', 'Chile','Ecuador','French Guiana','Guyana',
                  'Paraguay', 'Peru', 'Suriname',
                  'Trinidad and Tobago', 'Uruguay', 'Venezuela')
countries.maps <- map_data("world", region = sa.countries)

# Plot the sample locations on the world map using ggplot2
ggplot() + geom_polygon(data = countries.maps, aes(x = long, y = lat, group = group),
                        fill = "gray80", color = "gray20", size = 0.1) +
  geom_point(data = Data_Transect_selection, aes(x = longitude.y, y = latitude.y),
             size = 2,color="white") + geom_point(data = Data_Transect_BodyMass, color='red', aes(x = longitude.y, y = latitude.y),
                                                  size = 2)

# Plot the sample locations on the world map using ggplot2
ggplot() + geom_polygon(data = countries.maps, aes(x = long, y = lat, group = group),
                        fill = "gray80", color = "gray20", size = 0.1) +
  geom_point(data = Data_Transect_selection, aes(x = longitude.y, y = latitude.y),
             size = 2,color="black")



###2. Distribution among LU---------------
Check_Transect <- Data_Transect_selection %>%
  group_by(transect_id,biome,lu_cover1) %>%
  dplyr::summarise(Key_taxa = all(taxa_selection_included %in% Group))

p <- ggplot(data = Check_Transect, aes(x = lu_cover1,fill=lu_cover1)) +  geom_bar()+ facet_wrap(~biome)+theme_bw()
p+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Check_Transect <- Data_Transect_BodyMass %>%
  group_by(transect_id,biome,lu_cover1) %>%
  dplyr::summarise(Key_taxa = all(taxa_selection_included %in% Group))

p <- ggplot(data = Check_Transect, aes(x = lu_cover1,fill=lu_cover1)) +  geom_bar()+ facet_wrap(~biome)+theme_bw()
p+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

###3. Body mass, biomass, density---------------

p <- ggplot(data = Data_BodyMass_selectedGroup, aes(x = log10(BodyMass_g))) + geom_histogram()
p + facet_wrap(~Group)

p <- ggplot(data = Data_BodyMass_selectedGroup, aes(x = log10(BodyMass_g),y=log10(Density_sqm))) + geom_point()
p + facet_wrap(~Group)

p <- ggplot(data = Data_BodyMass_selectedGroup, aes(x = log10(Biomass_fresh_g_sqm),y=log10(Density_sqm))) + geom_point()
p + facet_wrap(~Group)

p <- ggplot(data = Data_BodyMass_selectedGroup, aes(x = log10(BodyMass_g),fill=biome)) +  geom_density(alpha=0.4)
p + facet_wrap(~Group)

toto<-subset(Data_BodyMass_selectedGroup,Group=="arachnida_total" & log10(BodyMass_g)< -4)
view(toto)

#remove very weird biomasses and body mass data
Data_BodyMass_selectedGroup$Biomass_fresh_g_sqm[which(Data_BodyMass_selectedGroup$datasetid=="MF_016")]<-NA
Data_BodyMass_selectedGroup$BodyMass_g[which(Data_BodyMass_selectedGroup$datasetid=="MF_016")]<-NA
Data_BodyMass_selectedGroup$Biomass_fresh_g_sqm[which(Data_BodyMass_selectedGroup$Biomass_fresh_g_sqm<0.0000001 & Data_BodyMass_selectedGroup$Density_sqm>0)]<-NA
Data_BodyMass_selectedGroup$BodyMass_g[which(Data_BodyMass_selectedGroup$BodyMass_g<0.0000001 & Data_BodyMass_selectedGroup$Density_sqm>0)]<-NA
Data_BodyMass_selectedGroup<-subset(Data_BodyMass_selectedGroup,!is.na(BodyMass_g))

p <- ggplot(data = Data_BodyMass_selectedGroup, aes(x = log10(BodyMass_g))) + geom_histogram()
p + facet_wrap(~Group)

p <- ggplot(data = Data_BodyMass_selectedGroup, aes(x = log10(BodyMass_g),y=log10(Density_sqm))) + geom_point()
p + facet_wrap(~Group)

p <- ggplot(data = Data_BodyMass_selectedGroup, aes(x = log10(Biomass_fresh_g_sqm),y=log10(Density_sqm))) + geom_point()
p + facet_wrap(~Group)

p <- ggplot(data = Data_BodyMass_selectedGroup, aes(x = log10(BodyMass_g),fill=biome)) +  geom_density(alpha=0.4)
p + facet_wrap(~Group)


#III. Mean body mass or data imputation of missing body mass -----

###1. File for biomass only----
Data_Transect_Biomass_selection <- Data_Transect_BodyMass

###2. Mean body mass-----
Data_Transect_Meanbodymass_selection <-Data_Transect_selection
Data_Transect_Meanbodymass_selection$Mean_BodyMass_g<-NA


for (i in taxa_selection_included){
  Data_Transect_Meanbodymass_selection$Mean_BodyMass_g[which(Data_Transect_Meanbodymass_selection$Group==i)]<-
    mean(Data_BodyMass_selectedGroup$BodyMass_g[which(Data_BodyMass_selectedGroup$Group==i & Data_BodyMass_selectedGroup$BodyMass_g>0)],na.rm=TRUE)
}

Data_Transect_Meanbodymass_selection$Biomass_fresh_g_sqm<-Data_Transect_Meanbodymass_selection$Density_sqm*Data_Transect_Meanbodymass_selection$Mean_BodyMass_g
Data_Transect_Meanbodymass_selection$BodyMass_g<-Data_Transect_Meanbodymass_selection$Mean_BodyMass_g

###3. Data imputation of missing body mass ----------
Data_Transect_selection_noZeros<-subset(Data_Transect_selection,Density_sqm>0)
Data_Transect_noZeros_BodyMass<-subset(Data_BodyMass_selectedGroup,Density_sqm>0 & BodyMass_g>0)

Data_Transect_selection_noZeros$lu_cover1[which(is.na(Data_Transect_selection_noZeros$lu_cover1))]<-"na"
Data_Transect_noZeros_BodyMass$lu_cover1[which(is.na(Data_Transect_noZeros_BodyMass$lu_cover1))]<-"na"

Select_enviro_variables_noNA<-Data_Transect_selection_noZeros %>%
  select(where(~!any(is.na(.))))
Select_enviro_variables_noNA2<-Data_Transect_noZeros_BodyMass %>%
  select(where(~!any(is.na(.))))
names(Select_enviro_variables_noNA)
names(Select_enviro_variables_noNA2)

selection_for_imputation<-c("biome","transect_id","latitude.y","longitude.y",
                            "lu_cover1","sample_year","temp","prec","elev","rad",
                            "evap","climKG","Group","Density_sqm","Biomass_fresh_g_sqm","BodyMass_g")

Data_Transect_selection_noZeros<-subset(Data_Transect_selection_noZeros,select=selection_for_imputation)
Data_Transect_noZeros_BodyMass<-subset(Data_Transect_noZeros_BodyMass,select=selection_for_imputation)

#imputation
Data_Final_imputed<-bodymass.imputation(Data_Transect_selection_noZeros,Data_Transect_noZeros_BodyMass,
                                        taxa_selection_included)

#here imputation BodyMass and Biomass independently, so two cases:
#either use imputation of BodyMass and calculate Biomass
#or use imputation of Biomass and calculate BodyMass
#we would need to compare both options later, for now use BodyMass imputation
Data_Final_imputed$Biomass_fresh_g_sqm<-Data_Final_imputed$Density_sqm*Data_Final_imputed$BodyMass_g

#IV. Build energy fluxes food webs ----------

###1. Select the dataset ----------
#Data with biomass available = Data_Transect_Biomass_selection
#Data with body mass estimated as mean body mass of the group = Data_Transect_Meanbodymass_selection
#Data with imputed body mass = Data_Final_imputed

Data<-Data_Transect_selection
Data<-subset(Data,Density_sqm>0)

#Load trait table for setting food web preferences
Raw_GuildsSoil <- read.table(file = "wilian/Raw_GuildsSoil_corrected.txt", header = T,
                           dec = ".",  # the decimal point
                           sep = "\t")  # the separator symbol

Raw_GuildsSoil
###2. Build preferences from Anton's approach ----------
Foodweb_matrices <-foodweb.preferences(Data,Raw_GuildsSoil)
Data_traits<-Foodweb_matrices$DataTraits
Data_Resources<-Foodweb_matrices$DataResources
Foodweb_matrices<-Foodweb_matrices$PreferenceMat

D###3. Energy fluxes estimation-----------
#Two cases possible:
#1. "Brose team" approach = only use metabolic as loss -> function foodweb.energyflux
#2. Add "de Ruiter" approach = also use natural death rates as loss (in addition to metabolic loss)
#   -> function foodweb.energyflux.plusdeath

#Go to Rob code



#Case 1:
Foodweb_energy<-foodweb.energyflux(Data_traits,Data_Resources,Foodweb_matrices)

#Case 2:
deathrates<-read.table("deathrates.txt",header=TRUE,sep="\t")
Foodweb_energy<-foodweb.energyflux.plusdeath(Data_traits,Data_Resources,Foodweb_matrices,deathrates)

Energyweb_matrices<-Foodweb_energy$EnergyMat
Firstdata<-Foodweb_energy$FirstOutput

Check_Transect <- Data %>%
  group_by(transect_id,biome,lu_cover1) %>%
  dplyr::summarise(Key_taxa = all(taxa_selection_included %in% Group))

Firstdata<-merge(Firstdata,Check_Transect,by="transect_id")
summary(Firstdata)
Firstdata$sumflux<-as.numeric(Firstdata$sumflux)
Firstdata$sumOM<-as.numeric(Firstdata$sumOM)
Firstdata$sumlitter<-as.numeric(Firstdata$sumlitter)
Firstdata$Ngroups<-as.numeric(Firstdata$Ngroups)
Firstdata$maxTL<-as.numeric(Firstdata$maxTL)
Firstdata$lu_cover1<-as.factor(Firstdata$lu_cover1)
p1 <- ggplot(data = Firstdata, aes(x = lu_cover1,y=log10(sumflux),fill=lu_cover1)) +
  geom_boxplot() + facet_wrap(~biome)+theme_bw()+ theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 <- ggplot(data = Firstdata, aes(x = lu_cover1,y=log10(sumOM),fill=lu_cover1)) +
  geom_boxplot() + facet_wrap(~biome)+theme_bw()+ theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 <- ggplot(data = Firstdata, aes(x = lu_cover1,y=log10(sumlitter),fill=lu_cover1)) +
  geom_boxplot() + facet_wrap(~biome)+theme_bw()+ theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p4 <- ggplot(data = Firstdata, aes(x = lu_cover1,y=Ngroups,fill=lu_cover1)) +
  geom_boxplot() + facet_wrap(~biome)+theme_bw()+ theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p5 <- ggplot(data = Firstdata, aes(x = lu_cover1,y=maxTL,fill=lu_cover1)) +
  geom_boxplot() + facet_wrap(~biome)+theme_bw()+ theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- ggplot(data = Firstdata, aes(x = lu_cover1,fill=lu_cover1)) +
  geom_bar()+ facet_wrap(~biome)+theme_bw()+ theme(legend.position = "none")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggarrange(p,p4,p5,p1,p2,p3,nrow=2,ncol=3)

