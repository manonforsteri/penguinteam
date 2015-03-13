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


#Si on suppose que dans I, m�me proportion de sexes que dans datas identifi�es
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


# on charge le package qui permet d appeler jags depuis R
library(R2jags)

# on specifie le modele en syntaxe BUGS 
# ou, par exemple, RE(time) signifie effet aleatoire (RE) du temps (time)
sink("cjs-msme.jags")
cat("
model{
#-- effets sur les parametres
	for (t in 1:(n.occasions-1)){
		logit(phi[t]) <-  mu.phi +  beta*t
		for (s in 1:nsitessex){
                  p[s,t] <- mu.p  
                            + eta[lsite[s]] 
                            + slopeLsex[s]*beta_p_Sex
      }
	}

	for (s in 1:nsites){
		eta[s] ~ dnorm(0,tau.p) # effet aleatoire site sur detection
	}


#-- priors

mean.phi ~ dunif(0, 1)             
mu.phi <- log(mean.phi / (1-mean.phi)) 

mean.p ~ dunif(0, 1)	           
mu.p <- log(mean.p / (1-mean.p)) 
sigma.p ~ dunif(0, 5)                
tau.p <- pow(sigma.p, -2)
sigma2.p <- pow(sigma.p, 2)

beta ~ dnorm (0, 0.0001)
beta_p_Sex ~ dnorm (0, 0.0001)
#-- vraisemblance multinomiale
	for (s in 1:nsitessex){ # boucle sur site X s
		for (t in 1:(n.occasions-1)){ # boucle sur temps
			marr[t,1:n.occasions,s] ~ dmulti(pr[s,t,], rel[t,s])
		}
		# Define the cell probabilities of the m-array:
		# Main diagonal
		for (t in 1:(n.occasions-1)){
			q[s,t] <- 1-p[s,t]
			pr[s,t,t] <- phi[t]*p[s,t]
			# Above main diagonal
			for (j in (t+1):(n.occasions-1)){
				pr[s,t,j] <- prod(phi[t:j])*prod(q[s,t:(j-1)])*p[s,j]
			} #j
			# Below main diagonal
			for (j in 1:(t-1)){
				pr[s,t,j]<-0
			} #j
		} #t
		# Last column: probability of non-recapture
		for (t in 1:(n.occasions-1)){
			pr[s,t,n.occasions] <- 1-sum(pr[s,t,1:(n.occasions-1)])
		} # t
	} # s
}
",fill = TRUE)
sink()

# met les donnees dans une liste
#nsites[k] ; nb sites par espèce k
#vecsp[s,k] ; pour chaque espèce k, on a les 2 indices qui donnent ou sont les sites, puis la dedans on peut sélectionner le site s
#rel[t,s,k] ; nb individus relâches au temps t, sur site s pour espèce k
jags.data <- list(marr = Marray, n.occasions = dim(Marray[,,1])[2]
                  , rel = apply(Marray,3,rowSums), nsites = lgr, nsitessex=length(slopeLsex),lsite=lsite,slopeLsex=slopeLsex)

# genere des valeurs initiales
inits <- function(){list(mean.phi = runif(1, 0, 1)
                         #, sigma.phi = runif(1, 0, 5)
                         , mean.p = runif(1, 0, 1)
                         #,sigma.p = runif(1, 0, 5)
                         #,sigma.alpha = runif(1, 0, 5)
                          )}  

# parametres que l on souhaite estimer
parameters <- c("phi","mean.p", "mean.phi"
                #, "sigma2.phi"
                #,"sigma2.p"
                #,"sigma2.alpha"
                #,"sigma2.beta"
                )

# specifications MCMC
ni <- 10 # nb iterations
nt <- 1 # nb thining
nb <- 1 # nb burnin
nc <- 3 # nb chains

# on appelle jags depuis R ; le hic est qu on a un bug, benin, mais un chtit bug
# ce bug est du type : 
#--------
# Erreur dans jags.model(model.file, data = data, inits = init.values, n.chains = n.chains,  : 
# Error in node marr[XX,XX:XX,XX]
# Invalid parent values
#--------
# 2 options
# i) on relance jags jusqu a ce que ca tourne
#temp.result <- jags(jags.data, inits, parameters, "cjs-msme.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
# ii) on boucle tant que le bug est la (on utilise pour ca la fn try), des que ca marche, on stocke les res

tmp <- jags(jags.data, inits, parameters, "cjs-msme.jags"
            , n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb
            , working.directory = getwd())

 while(inherits(try(tmp <- jags(jags.data, inits, parameters, "cjs-msme.jags"
                               , n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb
                               , working.directory = getwd()), silent=TRUE), "try-error")) {}

temp.result <- tmp 


# affiche les resultats
print(temp.result, digits = 3)

write(temp.result, file="resultsmodel3.txt")
