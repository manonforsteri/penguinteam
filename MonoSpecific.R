#----------------------------------------------------------
#--- modele de cormack-jolly-seber multi-site multi-espece
#----------------------------------------------------------
rm=list(ls())

setwd('C:/Users/Manon Ghislain L/Documents/Scripts/JAGS')

#---- charge les packages dont on aura besoin
library(R2jags) # pour appeler jags depuis R
library(FSA) # pour calculer les m-array (statistiques reduites)
library(abind) # pour combiner des arrays (generalise r/cbind)

#---- lit les donnees et les explore
dat <- read.table('datadefinitivesTURMER.txt',header=T,stringsAsFactors=T)

head(dat)
tail(dat)
names(dat)
dim(dat)

#dat[,"X.COV.espece"]<-NULL

attach(dat)

sum(S.) # nb total d individus
levels(X.COV.site) # modalites var site
levels(X.COV.sexe) # modalites var sexe
levels(X.COV.espece) # modalites var espece
nsites = length(levels(X.COV.site)) # nb sites 
nsex =  length(levels(X.COV.sexe)) # nb sexes
sexes<-levels(X.COV.sexe) # modalites var sexe

#---- met les donnees dans un format utilisable pour l analyse bayesienne
# pre-allocation memoire 
lsex = vector()
lsite = vector()
N = NULL # liste des marray par espece x site

uni<-1#unique marray index
M = list() # pre-alloue memoire pour les marray du sex courant
ind = 1
for (sex in 1:length(sexes))
{
  focsex <- sexes[sex]# recupere nom sex courant
  tempsex<-subset(dat,X.COV.sexe == focsex)
  for (s in 1:nsites)# boucle sur les sites
  { 
    i = levels(X.COV.site)[s] # recupere nom site courant
    temp <- subset(tempsex, X.COV.site == i) # selectionne le jeu de donnees pour site i espece k
    
    # les histoires sont groupees, degroupe pour revenir au format une ligne est un individu
    temp.converted <- capHistConvert(temp,freq="S.",in.type="frequency",out.type="individual")
    # convertit les histoires de capture en un m-array
    m.array.temp <- capHistSum(temp.converted,cols2use=-c(26:28))$m.array
    m.array.temp[is.na(m.array.temp)] = 0 # on remplace les NA par des 0
    m.array.temp <- m.array.temp[,-1] # on enleve la col nb de relaches
    #rowSums(m.array.temp) # pour info, calcule le nb de relaches par annee
    M[[ind]]= m.array.temp # stocke le marray pour le site courant dans la liste M
    lsex[uni]<-focsex
    lsite[uni]<-i
    ind <- ind+1
    uni <- uni+1
  }# end for (s in 1:nsites)
  lgr <- length(M) # stocke le nb de sites*sexes de l espece courante
  N <- c(N,M) # concatene les sites de l espece courante avec ceux des especes precedentes
}#end for (sex in 1:length(sexes))
#lgr[j] <- length(M) # stocke le nb de sites*sexes de l espece courante
#N <- c(N,M) # concatene les sites de l espece courante avec ceux des especes precedentes

# transforme la liste en array qui sera lu par JAGS
Marray=abind(M, along = 3)
dim(Marray) # 25 25  318 // 318 sites avec pour chacun un marray de dim 25x25

nsp<-1
vecsp = matrix(NA,nrow=max(lgr),ncol=1)
somcum = cumsum(lgr)
temp = 1:somcum[1]
vecsp[1:length(temp),1] = temp


#Si on suppose que dans I, m?me proportion de sexes que dans datas identifi?es
M<-dat[which(dat[,"X.COV.sexe"]=="M"),]
MALE<-sum(M[,"S."])
ALL<-sum(S.)
PROPMALES<-MALE/ALL

slopeLsex<-lsex
slopeLsex[which(slopeLsex=="I")]<-PROPMALES-0.5
slopeLsex[which(slopeLsex=="M")]<-0.5
slopeLsex[which(slopeLsex=="F")]<- -0.5

lsite<-rep(1: nsites, nsex)
#Si on suppose que dans I, 50% M et 50%F
#slopeLsex<-(as.numeric(as.factor(lsex))-2)/2 # -0.5 pour femelles, 0.5 pour males, et 0 pour indetermines (suppose un melange proche de 50%)

save.image(file = "TURMERarraysloaded.RData")
# on charge le package qui permet d appeler jags depuis R
library(R2jags)

### BUGS analysis

sink("cjs-mnl.jags")
cat("
    model {
    #-- effets sur les parametres
    for (t in 1:(n.occasions-1)){
    for (m in 1:nbmarrays){
      logit(phi[m,t]) <-  mu.phi+slopeLsex[vecsex[m]]*beta_phi_sex
      logit(p[m,t]) <- mu.p +slopeLsex[vecsex[m]]*beta_p_sex + eta[vecsite[m]] #c'est tordu vecsite[m] ici = s plus bas
    } 
    }
    # Priors and constraints
    
    mean.phi ~ dunif(0, 1)             
    mu.phi <- log(mean.phi / (1-mean.phi)) 
    mean.p ~ dunif(0, 1)             
    mu.p <- log(mean.p / (1-mean.p)) 
        
    beta_phi_sex~dnorm(0,0.01)
    beta_p_sex~dnorm(0,0.01)
    
    for (s in 1:nsites){
    eta[s] ~ dnorm(0,tau.p) # effet aleatoire site sur detection
    }#c'est tordu s ici = vecsite[m] plus haut
    sigma.p ~ dunif(0, 5)                
    tau.p <- pow(sigma.p, -2)
    sigma2.p <- pow(sigma.p, 2)
    

    # Define the multinomial likelihood
    for (m in 1:nbmarrays){ # boucle sur site*sexe
    for (t in 1:(n.occasions-1)){
    marr[t,1:n.occasions,m] ~ dmulti(pr[m,t, ], rel[t,m])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[m,t] <- 1-p[m,t]                # Probability of non-recapture
    pr[m,t,t] <- phi[m,t]*p[m,t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[m,t,j] <- prod(phi[m,t:j])*prod(q[m,t:(j-1)])*p[m,j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr[m,t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr[m,t,n.occasions] <- 1-sum(pr[m,t,1:(n.occasions-1)])
    } #t
    } #m
    }
    ",fill = TRUE)
sink()

Minput<-Marray[,,1:20]
nbmarrays<-dim(Minput)[3]
nsitesInput<-nbmarrays/2
lsiteInput<-rep(c(1:10),2)

#avec jeu de donn?es partiel:
jags.data <- list(marr = Minput, n.occasions = dim(Minput[,,1])[2], 
                  rel = apply(Minput,3,rowSums),
                  nbmarrays=nbmarrays,
                  vecsex=vecsp[1:nbmarrays],
                  nsites=nsitesInput,
                  slopeLsex=rep(c(-0.5,0.5),10),vecsite=lsiteInput)

#quand on utilise tout le jeu de donnees:
nbmarrays<-dim(Marray)[3]
jags.data <- list(marr = Marray, n.occasions = dim(Marray[,,1])[2], 
                  rel = apply(Marray,3,rowSums),
                  nbmarrays=nbmarrays,
                  vecsex=vecsp[,1],
                  nsites=nsites,
                  vecsite=lsite,slopeLsex=slopeLsex[1:nbmarrays])

# Initial values
inits <- function(){list(mean.phi = runif(1, 0, 1),
                        mean.p = runif(1, 0, 1),
                        beta_phi_sex=rnorm(1,0,sqrt(1/0.1)),
                        beta_p_sex=rnorm(1,0,sqrt(1/0.1)),
                        sigma.p=runif(n = 1,0,5))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p","beta_phi_sex","beta_p_sex","sigma.p")

# MCMC settings
ni <- 200
nt <- 3
nb <- 50
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs <- jags(jags.data, inits, parameters, "cjs-mnl.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(cjs, digits = 3) 
traceplot(cjs)

