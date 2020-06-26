##############################################################################################################################################
####### MAKROPHYTE DATA ANALYSE SCRIPT FOR FGG PUBLICATION ######
##############################################################################################################################################

## Set WD
setwd("C:/Users/anl85ck/Desktop/PhD/5_Macrophytes-Bavaria/4_FGG-Project/BavairanMacrophytesFGG")

#### Load Packages ####
#Analysis
library(tidyverse)

library(gamm4) #gamm analysis
library(BiodiversityR) #rankabund plot
#library(vegan) #for NMDS

#Plots
# library(ggpmisc)
library(ggpubr)  #customizing 'ggplot2'
library(ggrepel) #adds text directly to the plot
library(directlabels) #add direct labels to a plot, and hide the color legend
library(corrplot) #correlation plot

#Maps
library(raster) #maps
library(ggspatial) #maps
library(sf) #maps



##############################################################################################################################################
#### Load data ####
##############################################################################################################################################
# Macrophytes
Makroph_raw <- read.csv("./rawdata/Makrophyten_WRRL_05-17_nurMakrophytes.csv", header=TRUE, sep=";")
SPEC <- read.csv("./rawdata/Species_bis2017.csv", header=TRUE, sep=";")

# Morphology
morph <- read.csv("./rawdata/lake-data.csv",skip=0 , dec=",",header=TRUE, sep=";")

# Shapes
lakes_shape <- st_read("C:/Users/anl85ck/Desktop/PhD/3_Daten/Daten_Iffeldorf/Kartenmaterial/BAYERN/LfW/by_seen.shp")
lakes_bavaria <- subset(lakes_shape, (lakes_shape$ST == "9") | (lakes_shape$ST =="6"))
bavaria_shape <- st_read("C:/Users/anl85ck/Desktop/PhD/3_Daten/Daten_Iffeldorf/Kartenmaterial/BAYERN/LfW/by_grenz.shp")
cities <- st_read("C:/Users/anl85ck/Desktop/PhD/3_Daten/Daten_Iffeldorf/Kartenmaterial/BAYERN/Vektor/stadt_neu.shp") %>%
  st_set_crs(st_crs(lakes_shape))
rivers_shape <- st_read("C:/Users/anl85ck/Desktop/PhD/3_Daten/Daten_Iffeldorf/Kartenmaterial/BAYERN/LfW/fl_gew.shp")
rivers_important <- subset(rivers_shape, subset = rivers_shape$GEW_NAME %in% c("Donau", "Main","Roter Main","Weiäer Main","Lech",
                                                                               "Wertach","Wärnitz","Paar","Loisach","Amper","Wärm",
                                                                               "Main-Donau-Kanal","Mangfall","Alz","Tiroler Achen",
                                                                               "Salzach","Ilz","Naab","Isar","Mainkanal","Ammer",
                                                                               "Iller","Inn","Lech","Regen","Pegnitz","Altmähl",
                                                                               "Regnitz"))
rivers_important<-st_intersection(rivers_important,bavaria_shape)


# CHEM Data import
#Chem_table <- read.csv("./01_Input/Chem.Mean.Apr_AugDF.csv") %>% dplyr::rename(YEAR = Var2, Gewässer = Var1) #Import Chemical Dataset: Mean for the productive Summer Month between Apr - Aug
Chem_table <- read.csv("./rawdata/Chem.Mean.YearDF_ALL.csv") %>% dplyr::rename(YEAR = Var2, Gewässer = Var1) #Import Chemical Dataset: Mean for the productive Summer Month between Apr - Aug










##############################################################################################################################################
#### Data preparation Macrophytes ####
##############################################################################################################################################
# Calculate quantities
Makroph_raw$Messwert <- (Makroph_raw$Messwert)^3

#Selection of Lakes I am looking at
Makroph_raw <- Makroph_raw %>%
  filter(!(Gewässer=="Chiemsee" & (YEAR==2011))) %>% filter(!(Gewässer=="Chiemsee" & YEAR==2012)) %>% # 1 plot per year -> wrong
  filter(!(Gewässer=="Chiemsee" & (YEAR==2014))) %>% filter(!(Gewässer=="Chiemsee" & (YEAR==2015))) %>%
  filter(!(Gewässer=="Chiemsee" & (YEAR==2017))) %>% filter(!(Gewässer=="Staffelsee - Sued" & (YEAR==2012))) %>%
  filter(!(Gewässer=="Gr. Alpsee" & (YEAR==2012))) %>% filter(!(Gewässer=="Gr. Alpsee" & (YEAR==2013))) %>%
  filter(!(Gewässer=="Pilsensee" & (YEAR==2015))) %>% filter(!(Gewässer=="Langbuergner See" & (YEAR==2014))) %>%
  filter(!(Gewässer=="Pelhamer See" & (YEAR==2017))) %>%  filter(!(Gewässer=="Weitsee" & (YEAR==2017))) %>%
  filter(!(Gewässer=="Igelsbachsee" & (YEAR==2016))) %>%  filter(!(Gewässer=="Igelsbachsee" & (YEAR==2013))) %>%
  filter(!(Gewässer=="Grosser Brombachsee" & (YEAR==2015))) %>%  filter(!(Gewässer=="Rothsee" & (YEAR==2016))) %>%
  filter(!(Gewässer=="Altmuehlsee" & (YEAR==2013))) %>%
  distinct()

Makroph_raw$Probestelle <- plyr::revalue(Makroph_raw$Probestelle, c("0-1 m"="0-1", "1-2 m"="1-2", "2-4 m"="2-4",">4 m"="4-x" ))


## To get a dataset with all possible PLOTS
Makroph_dataset <- Makroph_raw %>% group_by(Gewässer, MST_NR, YEAR) %>%
  dplyr::select(Gewässer, MST_NR, YEAR) %>% #3590 * Gewässer, MST_NR, YEAR, Probestelle IIII 1013 *4 => 4052 mässtens eigentlich mind sein
  distinct()
Probestelle <- tibble(Probestelle=c("4-x","0-1","1-2", "2-4")) # tibble(Probestelle)
Makroph_dataset <- merge(Makroph_dataset, Probestelle, by=NULL) #996 * 4 = 3984

## Subselection of datasets
MakrophE <- Makroph_raw %>%
  filter(Form=="Em" | Form=="F-SB") %>%
  filter(str_detect(Taxon, " ")) %>%
  filter(Taxon != "Ranunculus, aquatisch")

MakrophSF <- Makroph_raw %>%
  filter(Form=="S") %>%
  filter(str_detect(Taxon, " ")) %>%
  filter(Taxon != "Ranunculus, aquatisch")

Makroph <- Makroph_raw %>%
  filter(Form!="NA") %>%
  filter(str_detect(Taxon, " ")) %>%
  filter(Taxon != "Ranunculus, aquatisch")

### PRODUCE COMMUNITY TABLE ALL
Makroph_comm_Comb2 <- Makroph %>% group_by(Gewässer, MST_NR, DATUM, Probestelle, Taxon) %>%
  filter(Messwert == max(Messwert)) %>% #get max value if we have multiple growth forms
  ungroup()%>% group_by(Gewässer, MST_NR, Probestelle, YEAR, Taxon) %>%
  summarise(Messwert = mean(Messwert)) %>% #get rid of double values for DATUM
  dplyr::select(Gewässer, MST_NR, YEAR, Probestelle, Taxon, Messwert) %>%  #duplicated %>% which %>% #check for duplications
  spread(Taxon, Messwert)%>%
  select_if(~sum(!is.na(.)) > 0)

Makroph_comm_Comb <-  dplyr::right_join(Makroph_comm_Comb2, Makroph_dataset, by=c("Gewässer", "MST_NR", "YEAR", "Probestelle"))
Makroph_comm_Comb[is.na(Makroph_comm_Comb)]<-0.0
Makroph_comm_Comb <- Makroph_comm_Comb[c(1:4,4 + which(colSums(Makroph_comm_Comb[-(c(1:4))]) !=0))]
Makroph_Lake_Trans <- Makroph_comm_Comb %>% group_by(Gewässer, MST_NR,YEAR) %>%
  summarise_at(vars(-group_cols(),-"Probestelle"), sum, na.rm=TRUE)
L<-length(Makroph_Lake_Trans)
Makroph_Lake_Trans$ALPHA <- specnumber(Makroph_Lake_Trans[4:L])

Makroph_Lake_ALPHA_SD <- Makroph_Lake_Trans %>% group_by(Gewässer, YEAR) %>% summarise(Alpha_SD = sd(ALPHA))
Makroph_Lake <- Makroph_Lake_Trans %>% group_by(Gewässer, YEAR) %>% mutate(N_Transects=n()) %>%
  summarise_at(vars(-group_cols(),-"MST_NR"), mean, na.rm=TRUE)
Makroph_Lake <- merge(Makroph_Lake, Makroph_Lake_ALPHA_SD)
Makroph_Lake$GAMMA <- specnumber(Makroph_Lake[3:(L-1)])


### PRODUCE COMMUNITY TABLE SUBMERGED
Makroph_comm_CombSF2 <- MakrophSF %>% group_by(Gewässer, MST_NR, DATUM, Probestelle, Taxon) %>%
  filter(Messwert == max(Messwert)) %>% #get max value if we have multiple growth forms
  ungroup()%>% group_by(Gewässer, MST_NR, Probestelle, YEAR, Taxon) %>%
  summarise(Messwert = mean(Messwert)) %>% #get rid of double values for DATUM
  dplyr::select(Gewässer, MST_NR, YEAR, Probestelle, Taxon, Messwert) %>%  #duplicated %>% which %>% #check for duplications
  spread(Taxon, Messwert)%>%
  select_if(~sum(!is.na(.)) > 0)

Makroph_comm_CombSF <-  dplyr::right_join(Makroph_comm_CombSF2, Makroph_dataset, by=c("Gewässer", "MST_NR", "YEAR", "Probestelle"))
Makroph_comm_CombSF[is.na(Makroph_comm_CombSF)]<-0.0
Makroph_comm_CombSF <- Makroph_comm_CombSF[c(1:4,4 + which(colSums(Makroph_comm_CombSF[-(c(1:4))]) !=0))]
Makroph_Lake_TransSF <- Makroph_comm_CombSF %>% group_by(Gewässer, MST_NR,YEAR) %>%
  summarise_at(vars(-group_cols(),-"Probestelle"), sum, na.rm=TRUE)
LSF<-length(Makroph_Lake_TransSF)
Makroph_Lake_TransSF$ALPHA_SF <- specnumber(Makroph_Lake_TransSF[4:LSF])

Makroph_Lake_ALPHA_SDSF <- Makroph_Lake_TransSF %>% group_by(Gewässer, YEAR) %>% summarise(Alpha_SD_SF = sd(ALPHA_SF))
Makroph_LakeSF <- Makroph_Lake_TransSF %>% group_by(Gewässer, YEAR) %>% mutate(N_Transects=n()) %>%
  summarise_at(vars(-group_cols(),-"MST_NR"), mean, na.rm=TRUE)
Makroph_LakeSF <- merge(Makroph_LakeSF, Makroph_Lake_ALPHA_SDSF)
Makroph_LakeSF$GAMMA_SF <- specnumber(Makroph_LakeSF[3:(LSF-1)])



### Community table emergent species
Makroph_comm_CombE2 <- MakrophE %>% group_by(Gewässer, MST_NR, DATUM, Probestelle, Taxon) %>%
  filter(Messwert == max(Messwert)) %>% #get max value if we have multiple growth forms
  ungroup()%>% group_by(Gewässer, MST_NR, Probestelle, YEAR, Taxon) %>%
  summarise(Messwert = mean(Messwert)) %>% #get rid of double values for DATUM
  dplyr::select(Gewässer, MST_NR, YEAR, Probestelle, Taxon, Messwert) %>%  #duplicated %>% which %>% #check for duplications
  spread(Taxon, Messwert)%>%
  select_if(~sum(!is.na(.)) > 0)

Makroph_comm_CombE <-  dplyr::right_join(Makroph_comm_CombE2, Makroph_dataset, by=c("Gewässer", "MST_NR", "YEAR", "Probestelle"))
Makroph_comm_CombE[is.na(Makroph_comm_CombE)]<-0.0
Makroph_comm_CombE <- Makroph_comm_CombE[c(1:4,4 + which(colSums(Makroph_comm_CombE[-(c(1:4))]) !=0))]
Makroph_Lake_TransE <- Makroph_comm_CombE %>% group_by(Gewässer, MST_NR,YEAR) %>%
  summarise_at(vars(-group_cols(),-"Probestelle"), sum, na.rm=TRUE)
LE<-length(Makroph_Lake_TransE)
Makroph_Lake_TransE$ALPHA_E <- specnumber(Makroph_Lake_TransE[4:LE])

Makroph_Lake_ALPHA_SDE <- Makroph_Lake_TransE %>% group_by(Gewässer, YEAR) %>% summarise(Alpha_SD_E = sd(ALPHA_E))
Makroph_LakeE <- Makroph_Lake_TransE %>% group_by(Gewässer, YEAR) %>% mutate(N_Transects=n()) %>%
  summarise_at(vars(-group_cols(),-"MST_NR"), mean, na.rm=TRUE)
Makroph_LakeE <- merge(Makroph_LakeE, Makroph_Lake_ALPHA_SDE)
Makroph_LakeE$GAMMA_E <- specnumber(Makroph_LakeE[3:(LE-1)])


### Combine Results
resultALL <- Makroph_Lake %>% dplyr::select(Gewässer, YEAR, GAMMA)
resultSF <- Makroph_LakeSF %>% dplyr::select(Gewässer, YEAR, GAMMA_SF) #N_Transects, ALPHA_SF, Alpha_SD_SF,
resultE <- Makroph_LakeE %>% dplyr::select(Gewässer,YEAR,  GAMMA_E) #ALPHA_E, Alpha_SD_E,
result <- merge(resultSF, resultE, by=c("Gewässer","YEAR"))
result <- merge(result, resultALL, by=c("Gewässer","YEAR"))
LAKES <- result$Gewässer

result <- merge(result, morph, by.x="Gewässer", by.y="Name_Makro_short")


##############################################################################################################################################
#### Data preparation Chemical data ####
##############################################################################################################################################

Chem <- Chem_table[c(2:5)] %>%  group_by(Gewässer, YEAR)%>% spread(key = Var3, value = value) #%>% select(-Chloroph..äg.l...0.0.m.Tiefe.)
Chem <- Chem[rowSums(is.na(Chem[,3:13]))<10,]
Chem_surf <- Chem %>% rowwise() %>%
  transmute(Gewässer, YEAR,
            Chloride = Chlorid..mg.l...0.0.m.Tiefe.,
            Conductivity = LF..20..U.00B0.C..vor.Ort...U.00B5.S.cm...0.0.m.Tiefe.,
            N_tot = N.ges...mg.l...0.0.m.Tiefe.,
            NH4N = NH4.N..mg.l...0.0.m.Tiefe.,
            NO3N = NO3.N..mg.l...0.0.m.Tiefe.,
            O2_diss = O2.gel.f6.st..mg.l...0.0.m.Tiefe.,
            P_tot = P.ges...mg.l...0.0.m.Tiefe.,
            pH = pH.Wert..vor.Ort.......0.0.m.Tiefe.,
            SiO2 = SiO2..mg.l...0.0.m.Tiefe.,
            Temp = Wassertemp..vor.Ort....U.00B0.C...0.0.m.Tiefe.,
            Transparency = Sichttiefe..cm...0.0.m.Tiefe.) %>%
  filter_at(vars(contains("_")), any_vars(!is.na(.)))

result_Chem <- merge(result, Chem_surf, by.x=c("Name_chem","YEAR"), by.y=c("Gewässer","YEAR"))

result_Chem <- result_Chem[complete.cases(result_Chem[c(8,10:11,13,32:42)]),]
result_Chem <- result_Chem %>% group_by(Gewässer) %>%  filter(YEAR==max(YEAR))

result_Chem$Fläche_log <- log10(result_Chem$Area_ha+1)
result_Chem$Tiefe_log <- log10(result_Chem$maxDepth_m+1)
result_Chem$Geländehähe_log <- log10(result_Chem$Altitude_masl+1)
result_Chem$LF_log <- log10(result_Chem$Conductivity+1)
result_Chem$Chlorid_log <- log10(result_Chem$Chloride+1)
result_Chem$N_log <- log10(result_Chem$N_tot+1)
result_Chem$NO3_log <- log10(result_Chem$NO3N+1)
result_Chem$NH4_log <- log10(result_Chem$NH4N+1)
result_Chem$O2_log <- log10(result_Chem$O2_diss+1)
result_Chem$P_log <- log10(result_Chem$P_tot+1)
result_Chem$pH_log <- log10(result_Chem$pH+1)
result_Chem$SiO2_log <- log10(result_Chem$SiO2+1)
result_Chem$Temp_log <- log10(result_Chem$Temp+1)
result_Chem$Sichttiefe_log <- log10(result_Chem$Transparency+1)

Fin_dataset<- result_Chem %>% dplyr::select(Gewässer, YEAR)

Lakes_Type <- merge(lakes_shape, morph, by.x = "SEE_NAME", by.y="Name_shape")
Lakes_Type <- merge(Lakes_Type, result_Chem, by.x = "Name_Makro_short", by.y="Gewässer")



##############################################################################################################################################
#### Data preparation rare species dataset ####
##############################################################################################################################################

#### RARE SPECIES DATASET | ALL
Makroph_Lake <- merge(Makroph_Lake, Fin_dataset, by=c("Gewässer", "YEAR"))
Makroph_Lake_h <-Makroph_Lake %>% dplyr::select("Acorus calamus":"Zannichellia palustris") %>%
  select_if(colSums(., na.rm=T) != 0)
Makroph_Lake <- cbind(Makroph_Lake[1:2], Makroph_Lake_h,Makroph_Lake[125:128])
PRESABS <- cbind(Makroph_Lake[,1:2], apply(Makroph_Lake[,3:96], 2, function(x) ifelse(x > 0.0, 1.0, x)))

SPECIESNUMBER_ALL <- PRESABS %>% dplyr::select(`Acorus calamus`:"Zannichellia palustris") %>% colSums()
sort(SPECIESNUMBER_ALL)
sum(SPECIESNUMBER_ALL==1)

boxplot(SPECIESNUMBER_ALL)

RARE <- PRESABS %>% dplyr::select(`Alisma plantago-aquatica`:"Zannichellia palustris") %>% select_if(colSums(.)<=3 )

RARE <- cbind(Makroph_Lake[,1:2], RARE)
RARE$RARESPECIES <- specnumber(RARE[3:52])

RARE_ALL <- merge(RARE, Lakes_Type, by.x=c("Gewässer"), by.y=c("Name_Makro_short"))

colnames(RARE)

##### RARE SPECIES DATASET Submerged
Makroph_LakeSF <- merge(Makroph_LakeSF, Fin_dataset, by=c("Gewässer", "YEAR"))
Makroph_LakeSF_h <-Makroph_LakeSF %>% dplyr::select("Alisma plantago-aquatica":"Zannichellia palustris") %>%
  select_if(colSums(., na.rm=T) != 0)
Makroph_LakeSF <- cbind(Makroph_LakeSF[1:2], Makroph_LakeSF_h,Makroph_LakeSF[94:97])
PRESABS_S <- cbind(Makroph_LakeSF[,1:2], apply(Makroph_LakeSF[,3:73], 2, function(x) ifelse(x > 0.0, 1.0, x)))

SPECIESNUMBER_S <- PRESABS_S %>% dplyr::select(`Alisma plantago-aquatica`:"Zannichellia palustris") %>% colSums()
sort(SPECIESNUMBER_S)
sum(SPECIESNUMBER_S==1)

RARE_S <- PRESABS_S %>% dplyr::select(`Alisma plantago-aquatica`:"Zannichellia palustris") %>% select_if(colSums(.)<=3 )

RARE_S <- cbind(Makroph_LakeSF[,1:2], RARE_S)
RARE_S$RARESPECIES_S <- specnumber(RARE_S[3:39])

RARE_ALL <- merge(RARE_ALL, RARE_S, by.x=c("Gewässer"), by.y=c("Gewässer"))


###### RARE SPECIES DATASET Emergent
Makroph_LakeE <- merge(Makroph_LakeE, Fin_dataset, by=c("Gewässer", "YEAR"))
Makroph_LakeE_h <-Makroph_LakeE %>% dplyr::select("Acorus calamus":"Typha latifolia") %>%
  select_if(colSums(., na.rm=T) != 0)
Makroph_LakeE <- cbind(Makroph_LakeE[1:2], Makroph_LakeE_h,Makroph_LakeE[55:58])
PRESABS_E <- cbind(Makroph_LakeE[,1:2], apply(Makroph_LakeE[,3:36], 2, function(x) ifelse(x > 0.0, 1.0, x)))
SPECIESNUMBER_E <- PRESABS_E %>% dplyr::select("Acorus calamus":"Typha latifolia") %>% colSums()
sort(SPECIESNUMBER_E)
sum(SPECIESNUMBER_E==1)

Lakes_Type_centers <- st_centroid(Lakes_Type)
Lakes_Type_centers <- cbind(Lakes_Type, st_coordinates(st_centroid(Lakes_Type$geometry)))

usethis::use_data(Makroph_Lake,overwrite = TRUE)
usethis::use_data(Makroph_LakeSF,overwrite = TRUE)
usethis::use_data(Makroph_LakeE,overwrite = TRUE)
usethis::use_data(PRESABS,overwrite = TRUE)
usethis::use_data(PRESABS_S,overwrite = TRUE)
usethis::use_data(PRESABS_E,overwrite = TRUE)

usethis::use_data(result_Chem,overwrite = TRUE)

usethis::use_data(Lakes_Type_centers,overwrite = TRUE)
usethis::use_data(bavaria_shape,overwrite = TRUE)
usethis::use_data(lakes_bavaria,overwrite = TRUE)
usethis::use_data(rivers_important,overwrite = TRUE)
usethis::use_data(cities,overwrite = TRUE)




##############################################################################################################################################
##############################################################################################################################################
#### RESULTS | PLOTS ####
##############################################################################################################################################
##############################################################################################################################################
###################################################################################
###############################################################################################
###############################################################################################
# MAPS
###############################################################################################
###############################################################################################
# JUST BAVARIA










###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
# # ELODEAS MAP
#
#
# df <- cbind(Makroph_LakeSF[,1:2], apply(Makroph_LakeSF[,3:73], 2, function(x) ifelse(x > 0.0, 1.0, x)))
#
# #Häufigste Arten
# df %>% select(`Alisma plantago-aquatica`:"Zannichellia palustris") %>% colSums()
#
#
# SUBM<-merge(df,result_Chem,by=c("Gewässer", "YEAR"))
#SUB_SHAPE<-merge(Lakes_Type,df,by.y=c("Gewässer"),by.x=c("Lake.y"))
#
#



# ##########################################################################
# ## Species distribution maps - SUBMERGED
#
# df <- cbind(Makroph_LakeSF[,1:2], apply(Makroph_LakeSF[,3:73], 2, function(x) ifelse(x > 0.0, 1.0, x)))
#
# #Häufigste Arten
# df %>% select(`Alisma plantago-aquatica`:"Zannichellia palustris") %>% colSums()
#
#
# SUBM<-merge(df,result_Chem,by=c("Gewässer", "YEAR"))
# SUB_SHAPE<-merge(Lakes_Type,df,by.y=c("Gewässer"),by.x=c("Lake.y"))
#
#
# ### SUBMERGED
# species <- colnames(Makroph_LakeSF[3:73])
#
# for (i in 1:length(species)){
#
#   S<-ggplot(SUB_SHAPE)+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white")+
#     geom_sf(data = rivers_important, size = 1, color = "lightblue") +
#     geom_sf(aes(fill = as.factor(SUB_SHAPE[[89+i]])),color=NA) +
#     scale_fill_manual(values = c("grey", "red"),name= "Präsenz")+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", alpha=0)+
#     annotation_scale(location = "br", width_hint = 0.25)+
#     geom_sf(data=cities)+
#     #geom_sf_text(aes(label = YEAR.y), size=5,nudge_x = 3000, nudge_y = 000,check_overlap=T) +
#     geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.0,nudge_y = -3000.5,check_overlap=T)+
#     coord_sf(xlim = c(4355000,4585000),ylim= c(5250000,5340000)) + xlab("")+ylab("")+
#     theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
#                                         axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#                                         panel.border = element_rect(colour = "black", fill=NA, size=2))
#
#   M <- ggplot(SUB_SHAPE)+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white")+
#     geom_sf(data = rivers_important, size = 1, color = "lightblue") +
#     geom_sf(aes(fill = as.factor(SUB_SHAPE[[89+i]])),color=NA) +
#     scale_fill_manual(values = c("grey", "red"),name= "Präsenz")+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", alpha=0)+
#     annotation_scale(location = "br", width_hint = 0.25)+#ggtitle(species[i])+
#     geom_sf(data=cities)+
#     geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.5,nudge_y = -3000.5,check_overlap=T)+
#     coord_sf(xlim = c(4400000,4563000),ylim= c(5411000,5485000)) +xlab("")+ylab("") +theme(legend.position="none")+
#     theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
#           axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#           panel.border = element_rect(colour = "black", fill=NA, size=2))
#
#   N <- ggplot(SUB_SHAPE)+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white")+
#     geom_sf(data = rivers_important, size = 1, color = "lightblue") +
#     geom_sf(aes(fill = as.factor(SUB_SHAPE[[89+i]])),color=NA) +
#     scale_fill_manual(values = c("grey", "red"),name= "Präsenz")+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", alpha=0)+
#     annotation_scale(location = "br", width_hint = 0.25)+#ggtitle("Elodea canadensis")+
#     geom_sf(data=cities)+
#     geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.5,nudge_y = -1000.5,check_overlap=T)+
#     coord_sf(xlim = c(4491000,4500000),ylim= c(5568000,5579000)) +xlab("")+ylab("")+
#     theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
#                                         axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#                                         panel.border = element_rect(colour = "black", fill=NA, size=2))
#
#   all<-ggarrange(ggarrange(N,M, ncol=2,widths = c(0.8, 2.1),labels =c("A","B")),S,nrow=2,common.legend=T,
#                  labels =c("","C"),legend="bottom",heights =  c(1.3, 2))
#
#   plot<-annotate_figure(all,top = text_grob(species[i], color = "black", face = "bold", size = 14))
#
#   jpeg(file = paste("./03_Results/species_maps/", species[i], '.jpeg', sep = ''), width = 22, height = 18, res=300, units="cm")
#   print(plot)
#   dev.off()
# }
#
#
# ##########################################################################
# ## Species distribution maps - EMERGENT & FLOATING & FL-LEAF
#
# species_E <- colnames(Makroph_LakeE[3:36])
#
# dfE <- cbind(Makroph_LakeE[,1:2], apply(Makroph_LakeE[,3:36], 2, function(x) ifelse(x > 0.0, 1.0, x)))
#
# #Häufigste Arten
# dfE %>% select(`Acorus calamus`:"Typha latifolia") %>% colSums() %>% sort
#
#
# SUBME<-merge(dfE,result_Chem,by=c("Gewässer", "YEAR"))
# SUB_SHAPEE<-merge(Lakes_Type,dfE,by.y=c("Gewässer"),by.x=c("Lake.y"))
#
# for (i in 1:length(species_E)){
#
#   S<-ggplot(SUB_SHAPEE)+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white")+
#     geom_sf(data = rivers_important, size = 1, color = "lightblue") +
#     geom_sf(aes(fill = as.factor(SUB_SHAPEE[[86+i]])),color=NA) +
#     scale_fill_manual(values = c("grey", "red"),name= "Präsenz")+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", alpha=0)+
#     annotation_scale(location = "br", width_hint = 0.25)+
#     geom_sf(data=cities)+
#     #geom_sf_text(aes(label = YEAR.y), size=5,nudge_x = 3000, nudge_y = 000,check_overlap=T) +
#     geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.0,nudge_y = -3000.5,check_overlap=T)+
#     coord_sf(xlim = c(4355000,4585000),ylim= c(5250000,5340000)) + xlab("")+ylab("")+
#     theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
#                                         axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#                                         panel.border = element_rect(colour = "black", fill=NA, size=2))
#
#   M <- ggplot(SUB_SHAPEE)+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white")+
#     geom_sf(data = rivers_important, size = 1, color = "lightblue") +
#     geom_sf(aes(fill = as.factor(SUB_SHAPEE[[86+i]])),color=NA) +
#     scale_fill_manual(values = c("grey", "red"),name= "Präsenz")+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", alpha=0)+
#     annotation_scale(location = "br", width_hint = 0.25)+#ggtitle(species[i])+
#     geom_sf(data=cities)+
#     geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.5,nudge_y = -3000.5,check_overlap=T)+
#     coord_sf(xlim = c(4400000,4563000),ylim= c(5411000,5485000)) +xlab("")+ylab("") +theme(legend.position="none")+
#     theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
#           axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#           panel.border = element_rect(colour = "black", fill=NA, size=2))
#
#   N <- ggplot(SUB_SHAPEE)+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white")+
#     geom_sf(data = rivers_important, size = 1, color = "lightblue") +
#     geom_sf(aes(fill = as.factor(SUB_SHAPEE[[86+i]])),color=NA) +
#     scale_fill_manual(values = c("grey", "red"),name= "Präsenz")+
#     geom_sf(data = bavaria_shape, size = 1, color = "black", alpha=0)+
#     annotation_scale(location = "br", width_hint = 0.25)+#ggtitle("Elodea canadensis")+
#     geom_sf(data=cities)+
#     geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.5,nudge_y = -1000.5,check_overlap=T)+
#     coord_sf(xlim = c(4491000,4500000),ylim= c(5568000,5579000)) +xlab("")+ylab("")+
#     theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
#                                         axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
#                                         panel.border = element_rect(colour = "black", fill=NA, size=2))
#
#   all<-ggarrange(ggarrange(N,M, ncol=2,widths = c(0.8, 2.1),labels =c("A","B")),S,nrow=2,common.legend=T,
#                  labels =c("","C"),legend="bottom",heights =  c(1.3, 2))
#
#   plot<-annotate_figure(all,top = text_grob(species_E[i], color = "black", face = "bold", size = 14))
#
#   jpeg(file = paste("./03_Results/species_maps_E/", species_E[i], '.jpeg', sep = ''), width = 22, height = 18, res=300, units="cm")
#   print(plot)
#   dev.off()
# }

### INVASIVE
# SUB_SHAPE$`Elodea canadensis`
# ggplot(SUB_SHAPE, aes(group=`Elodea canadensis`, y=log(Area_ha.x)))+geom_boxplot()
# ggplot(SUB_SHAPE, aes(group=`Elodea canadensis`, y=log(Temp)))+geom_boxplot()
#
# C<-cor(SUBM[c(22,23,118:131)])
# corrplot(C, type = "upper", method = "number", insig = "blank")
#
# corrplot(cor(SUBM[c(22,23,118:131)])[1:2,1:16, drop=FALSE], method = "number",cl.pos='n')


# #########################################################################################
# ####### RARE SPECIES MAP
# #Lakes_Type_centers <- st_centroid(RARE_ALL)
# RARE_ALL_centers <- cbind(RARE_ALL, st_coordinates(st_centroid(RARE_ALL$geometry)))
#
# x1<-ggplot(data = subset(RARE_ALL_centers, Region.x=="Alpine") , aes(geometry=geometry)) +
#   geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white") +
#   geom_sf(data = rivers_important, size = 0.05, color = "lightblue")  +
#   geom_sf(data=cities)+
#   geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.5,nudge_y = -3000.5,check_overlap=F) +
#   geom_sf(data = subset(RARE_ALL_centers, Region.x=="Alpine"),aes(fill = RARESPECIES), color=NA) +
#
#   scale_fill_gradient(low = "grey", high = "darkgreen", name="N species",limits = c(0,10))+
#   annotation_scale(location = "br", width_hint = 0.25) +
#   geom_text_repel(mapping=aes(x=X, y=Y,label=paste(RARESPECIES_S ,"|",RARESPECIES_E),color="black"),
#                   size=4, box.padding = unit(0.5, "lines"),color="black")+
#   #scale_color_continuous(low = "white", high = "purple", name="N species",limits = c(0,36))+
#   coord_sf(xlim = c(4355000,4585000),ylim= c(5250000,5340000), expand = FALSE)+
#   xlab("") + ylab("")+
#   theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=1))
#
#
# x2<-ggplot(data = subset(RARE_ALL_centers, Region.x=="Low mountain") , aes(geometry=geometry)) +
#   geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white") +
#   geom_sf(data = rivers_important, size = 0.05, color = "lightblue")  +
#   geom_sf(data=cities)+
#   geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.5,nudge_y = -3000.5,check_overlap=F) +
#   #geom_sf(aes(fill = RARE_ALL_centers), color=NA) +
#   geom_sf(data = subset(RARE_ALL_centers, Region.x=="Low mountain"),aes(fill = RARESPECIES), color=NA) +
#   scale_fill_gradient(low = "grey", high = "darkgreen", name="N species",limits = c(0,10))+
#   annotation_scale(location = "br", width_hint = 0.25) +
#   geom_text_repel(mapping=aes(x=X, y=Y,label=paste(RARESPECIES_S ,"|",RARESPECIES_E),color="black"),
#                   size=4, box.padding = unit(0.5, "lines"),color="black")+
#   #scale_color_continuous(low = "white", high = "purple", name="N species",limits = c(0,36))+
#   coord_sf(xlim = c(4400000,4563000),ylim= c(5411000,5485000)) +
#   xlab("") + ylab("")+
#   theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=1))
#
# x3<-ggplot(data = subset(RARE_ALL_centers, Region.x=="Low mountain") , aes(geometry=geometry)) +
#   geom_sf(data = bavaria_shape, size = 1, color = "black", fill = "white") +
#   geom_sf(data = rivers_important, size = 0.05, color = "lightblue")  +
#   geom_sf(data=cities)+
#   geom_sf_text(data=cities, aes(label = cities$name), nudge_x = 0.5,nudge_y = -500.5,check_overlap=F) +
#   geom_sf(aes(fill = RARESPECIES), color=NA) +
#
#   scale_fill_gradient(low = "grey", high = "darkgreen", name="N species",limits = c(0,10))+
#   annotation_scale(location = "br", width_hint = 0.25) +
#   geom_text_repel(mapping=aes(x=X, y=Y,label=paste(RARESPECIES_S ,"|",RARESPECIES_E),color="black"),
#                   size=4, box.padding = unit(0.5, "lines"),color="black")+
#   #scale_color_continuous(low = "white", high = "purple", name="N species",limits = c(0,36))+
#   coord_sf(xlim = c(4489000,4500000),ylim= c(5569000,5578000)) +
#   xlab("") + ylab("")+
#   theme(legend.position="none",axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=1))
#
#
#
#
#
#
# XALL <- ggarrange(ggarrange(x3,x2,ncol=2,widths = c(1, 2), labels =c("A","B"), common.legend=T,legend="right"),x1,nrow=2,labels =c("","C"),
#                   heights =  c(1.5, 2) )
# XALL <- annotate_figure(XALL,top ="Seltene Arten: submers (N) | emers (N)")
# XALL
#
#
# #RARE_ALL$RARESPECIES <- RARE_ALL$RARESPECIES_S + RARE_ALL$RARESPECIES_E
# dev.off()
# corrplot(cor(RARE_ALL[c(97,95,96,53,187,212,134:147)])[1:6,1:20, drop=FALSE], method = "number",cl.pos='n')
#
# gam_r_allE <- gam(RARESPECIES ~ s(Area_ha_log, bs = 'cr', k = 3)+s(Transparency_log, bs = 'cr', k = 3), #+s(Transparency_log, bs = 'cr', k = 3)+s(N_tot_log, bs = 'cr', k = 3)+
#                                       data=RARE_ALL,method="REML")
# summary(gam_r_allE)
#
#
# gam_r_allE <- gam(RARESPECIES_S ~ s(Area_ha_log, bs = 'cr', k = 3)+
#
#                                         s(O2_diss_log, bs = 'cr', k = 3), #+s(Transparency_log, bs = 'cr', k = 3)+s(N_tot_log, bs = 'cr', k = 3)+
#                                       data=RARE_ALL,method="REML")
# summary(gam_r_allE)
#
# # gam_r_allE <- gam(RARESPECIES ~ s(Area_ha_log, bs = 'cr', k = 3)+s(maxDepth_log, bs = 'cr', k = 3)+s(Altitude_masl_log, bs = 'cr', k = 3)+
# #
# #                     s(Chloride_log, bs = 'cr', k = 3)+s(Conductivity_log, bs = 'cr', k = 3)+
# #                     s(SiO2_log, bs = 'cr', k = 3)+s(NO3N_log, bs = 'cr', k = 3)+s(NH4N_log, bs = 'cr', k = 3)+
# #                     s(O2_diss_log, bs = 'cr', k = 3)+s(P_tot_log, bs = 'cr', k = 3)+s(pH_log, bs = 'cr', k = 3)+
# #                     s(Temp_log, bs = 'cr', k = 3), #+s(Transparency_log, bs = 'cr', k = 3)+s(N_tot_log, bs = 'cr', k = 3)+
# #                   data=RARE_ALL,method="REML")
#
# # gam_r_allE <- gam(RARESPECIES ~ s(Area_ha_log, bs = 'cr', k = 3)+s(maxDepth_log, bs = 'cr', k = 3)+
# #                     s(O2_diss_log, bs = 'cr', k = 3)+
# #                     s(Transparency_log, bs = 'cr', k = 3), #s(Area_ha_log)+
# #                   data=RARE_ALL,method="REML")
# # summary(gam_r_allE)
# par(mfrow=c(1,2))
# plot(gam_r_allE, rug=T)
# vis.gam(gam_r_allE,view=c("Area_ha_log","O2_diss_log"),plot.type="contour",color="bw")
#
#
#
#
