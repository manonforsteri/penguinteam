library(R2jags)
library(Hmisc)
Sys.setenv(MAKEFLAGS = "-j8")  

setwd("C:/Users/Pierre-Yves Henry/Desktop/nouveau JAGS")
spagarder<-read.table("spagarder.txt")
data12<-read.table("data12.txt", sep = " ",header=T, na="", dec=".")
DATA<-data12[which(substr(data12$ETAT,1,2)=="AD" 
                   & data12$ANNEE<2015 
                   & data12$meansitesp>=5
                   & data12$SP %in% spagarder$SP),]

##SUPPR SP ZONES HUMIDES AVANT 2001
hum<-c("PANBIA","CETCET","EMBSCH","LUSSVE", "ACRSCI", "ACRRIS", "ACRSCH") 
DATA<-DATA[which(DATA$SP %nin% hum
                 | (DATA$SP %in% hum & DATA$ANNEE<2001)     ),]

#ne garder un site que s'il a au moins 2 espèces : 
tab<-table(DATA$SITE, DATA$SP)
tab2=(tab!=0)
siteretenu=rownames(tab2)[rowSums(tab2)>=2]
DATA<-DATA[(DATA$SITE %in% siteretenu),]
DATA<-droplevels(DATA)

#matrice 0 et 1
ind<-unique(DATA$INDSITE)
CH= table(DATA$INDSITE, DATA$ANNEE)
CH<-as.data.frame.matrix(CH)
CH<-CH[which(rownames(CH) %in% ind),]
CH[CH!="0"]<-1
CH<-as.matrix(CH)

#nombre d'individus
I<-nrow(CH)
#nombre max d'occasion
Kmax<-(ncol(CH))


#liste des individus et leur site/espèce/sexe
individuals<-as.data.frame(rownames(CH))
individuals$site<-DATA$SITE[match(individuals[,1],DATA$INDSITE)]
individuals$sex<-DATA$SEXE2[match(individuals[,1],DATA$INDSITE)]
individuals$covhetero<-DATA$COVHETERO[match(individuals[,1],DATA$INDSITE)]
individuals$species<-DATA$SP[match(individuals[,1],DATA$INDSITE)]

#Sex effect
#sexe
sex<-as.data.frame(individuals$sex)
sex[,1]<-as.character(sex[,1])
sex[,1][sex[,1]=="M"]<-"aM"
list_sex<-as.data.frame(sort(unique(sex[,1])))
list_sex$row<-row.names(list_sex)
sex$numsex<-list_sex$row[match(sex[,1],list_sex[,1])]
sex<-sex[,2]
sex<-as.numeric(sex)
#nsexe
nsex<-as.numeric(nrow(list_sex))

#site effect
#site
site<-as.data.frame(individuals$site)
count_site<-stack(table(site))
colnames(count_site) <- c("nb","site")
site[,1]<-as.character(site[,1])
site_ref<-as.character(count_site$site[count_site$nb==max(count_site$nb)][1])
site[,1][site[,1]==site_ref]<-paste("0",site_ref,sep="")
list_site<-as.data.frame(sort(unique(site[,1])))
list_site$row<-row.names(list_site)
site$numsite<-list_site$row[match(site[,1],list_site[,1])]
site<-site[,2]
site<-as.numeric(as.character(site))
#nsite
nsite<-as.numeric(nrow(list_site))

#espece
species<-as.data.frame(individuals$species)

count_species<-stack(table(species))
colnames(count_species) <- c("nb","species")
species[,1]<-as.character(species[,1])
species_ref<-as.character(count_species$species[count_species$nb==max(count_species$nb)][1])
species[,1][species[,1]==species_ref]<-paste("aa",species,sep="")

list_species<-as.data.frame(sort(unique(species[,1])))
list_species$row<-row.names(list_species)
species$numspecies<-list_species$row[match(species[,1],list_species[,1])]
species<-species[,2]
species<-as.numeric(species)
#nespece
nspecies<-as.numeric(nrow(list_species))

#hétérogénéité individuelle
covhetero<- individuals$covhetero

#Create vector with occasion of marking

get.first<-function(x)min(which(x!=0))
f<-apply(CH,1,get.first)

#matrices "âges" : inc = inconnu = âge 1, res=resident = âge 2
x<-matrix(NA, ncol=dim(CH)[2] , nrow=dim(CH)[1])

for (i in 1:nrow(x)) {
  for(t in f[i]: (ncol(x)) ){  
    x[i,t]<-2
    x[i,f[i]]<-1
  }
}

#Construire la matrice pour les dummy variables sur p
DATA$SITEANNEE<-paste(DATA$SITE, "       ", DATA$ANNEE)
tab<-table(DATA$SITE, DATA$ANNEE)
tab<-as.data.frame(tab)
tab$SITEANNEE<-paste(tab$Var1, "       ", tab$Var2)
tab2<-tab[which(tab$Freq!=0),]


X<-table(DATA$INDSITE, DATA$SITEANNEE)
X<-as.data.frame(X)
X$ANNEE<-tab$Var2[match(X$Var2,tab$SITEANNEE)]
X$CAPT<-NA
X$CAPT[X$Var2 %in% tab2$SITEANNEE]<-1
X<-X[which(X$CAPT==1),]
X<-X[which(substr(X$Var1, 1, 5) == substr(X$Var2, 1, 5)),]
A<-table(X$Var1, X$ANNEE)
