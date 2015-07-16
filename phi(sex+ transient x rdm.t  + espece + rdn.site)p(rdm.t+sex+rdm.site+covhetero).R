library(R2jags)
library(Hmisc)
Sys.setenv(MAKEFLAGS = "-j8")  


setwd("C:/Users/MGGCBRUNOY/Desktop/nouveau JAGS")
data12<-read.table("data12.txt", sep = " ",header=T, na="", dec=".")
DATA<-data12[which(substr(data12$ETAT,1,2)=="AD" 
                   & data12$ANNEE<2014 
                   & data12$meansitesp>=80),]


##SUPPR SP ZONES HUMIDES AVANT 2001
hum<-c("PANBIA","CETCET","EMBSCH","LUSSVE", "ACRSCI", "ACRRIS", "ACRSCH") 
DATA<-DATA[which(DATA$SP %nin% hum
              | (DATA$SP %in% hum & DATA$ANNEE<2001)     ),]

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

A<-matrix(data = 1,nrow = I,ncol=Kmax)       #1 ou 0 si année avec baguage ou non

sink("cjs.jags")
cat("
    
    model{
    
    #Model
    for(i in 1:I){
      for(t in (f[i]):(Kmax-1)) {
        logit(phi[i,t])<- rt.phi.trans[x[i,t],t] + phi.sex[sex[i]] + phi.species[species[i]] + random.phi.site[site[i]]
        lp[i,t]<-  random.p.temps[t] + p.sex[sex[i]] + random.p.site[site[i]] + b1*covhetero[i] + p.species[species[i]]
        p[i,t]<-   A[i,t]/(1+exp(-lp[i,t]))                         #inverse logit
      
      } #t
    }#i
    
    #####Priors#####
    

    #for survival parameters

    phi.sex[1] <- 0
    phi.sex[2] ~ dnorm(0,0.01)I(-10,10)
    phi.sex[3] ~ dnorm(0,0.01)I(-10,10)
    
    for(a in 1:2){
      for(t in 1:(Kmax-1)){
        rt.phi.trans[a,t]<-mu.phi.rttrans[a] + epsilon.phi.rttrans[a,t]
        epsilon.phi.rttrans[a,t] ~ dnorm(0, tau.phi.rttrans[a])
      }
      mu.phi.rttrans[a] ~ dnorm(0,0.01)I(-10,10)         
      sigma.phi.rttrans[a] ~ dunif(0,10) 
      tau.phi.rttrans[a] <- pow(sigma.phi.rttrans[a], -2)
    }

    phi.species[1]<- 0
    for(s in 2:nspecies){
    phi.species[s] ~ dnorm(0,0.01)I(-10,10)
    }

    for(u in 1:nsite){                                             
     random.phi.site[u] ~ dnorm(0, tau.phi.site)
    }    
    sigma.phi.site   ~ dunif(0,10)
    tau.phi.site    <-  pow(sigma.phi.site, -2) 


    #for recapture parameters
   
    b1 ~ dnorm(0,0.01)I(-10,10)
    
    phi.sex[1] <- 0
    for (s in 2:3){
        phi.sex[s] ~ dnorm(0,0.01)I(-10,10)
      }
    
    for(u in 1:nsite){                                          
      random.p.site[u] ~ dnorm(0, tau.p.site)                      
    }
    sigma.p.site     ~ dunif(0,10)
    tau.p.site       <-  pow(sigma.p.site, -2)    
    
    p.species[1]<- 0
    for(s in 2:nspecies){
      p.species[s] ~ dnorm(0,0.01)I(-10,10)
    }
    
    for(t in 1:(Kmax-1)){
      random.p.temps[t] ~ dnorm(0,tau.p.temps)
    }
    sigma.p.temps     ~ dunif(0,10)
    tau.p.temps       <-  pow(sigma.p.temps, -2)

    ####Likelihood####
    for(i in 1:I){
    
      #Latent state at first capture
      z[i,f[i]]<-1
      
      for(t in (f[i]+1):Kmax){
        
        #State process
        z[i,t] ~ dbern(mu1[i,t])
        mu1[i,t]<-phi[i,t-1]*z[i,t-1]
        
        #Obs process
        y[i,t] ~ dbern(mu2[i,t])
        mu2[i,t]<-p[i,t-1]*z[i,t]
      }#t
    }#i   
    
    }
    ",fill=TRUE)
sink()

#Functiun to create a matrix with information about known latent state z
known.state.cjs<-function(ch){
  state<-ch
  for(i in 1:dim(ch)[1]){
    n1<-min(which(ch[i,]==1))
    n2<-max(which(ch[i,]==1))
    state[i,n1:n2]<-1
    state[i,n1]<-NA
  }
  state[state==0]<-NA
  return(state)
}

#Bundle.data
data <- list(I=nrow(CH), Kmax=ncol(CH), f=f, y=CH, sex=sex, site=site, nsite=nsite
                   , z=known.state.cjs(CH),covhetero=covhetero, x=x , species=species, nspecies = nspecies,
                     A=A)

#Matrix of initial values for latent state z
cjs.init.z<-function(ch,f){
  
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2<- max (which(ch[i,]==1))
    ch[i,f[i]:n2]<-NA  
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]]<-NA
  }
  return(ch)
}


#Initial values
inits<-function(){
  list(  
    z=cjs.init.z(CH,f),
    phi.sex=c(NA,rnorm(nsex-1)),
    p.sex=c(NA,rnorm(nsex-1)),
    b1=rnorm(1),
    p.species=c(NA,rnorm(nspecies-1)),
    phi.species=c(NA,rnorm(nspecies-1))
  )}


parameters<-c("phi.sex" , "p.sex"  ,"b1"  ,"p.species","phi.species",                              #effets fixes
              "random.phi.site", "random.p.site"    , "sigma.phi.site", "sigma.p.site",
              "random.p.temps" , "sigma.p.temps"    ,
              "rt.phi.trans"   , "sigma.phi.rttrans")

#modele:
modA<- jags(data=data, inits=inits, parameters=parameters, "cjs.jags", n.thin=10, n.chains=3, n.burnin=500, n.iter=1000)

save(modmultispeA, file="C:/Users/MGGCBRUNOY/Desktop/nouveau JAGS/modmultispeB")
