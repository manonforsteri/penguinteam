#----------------------------------------------------------
#--- modele de cormack-jolly-seber multi-site multi-espece
#--- note : on ignore le sexe
#----------------------------------------------------------
#rm=list(ls())
#setwd('C:/Users/Manon Ghislain L/Documents/Scripts/JAGS')
setwd('C:/Users/Thimothee Admin/Documents/GitHub/penguinteam/')
#setwd('C:/Users/Timoth�e/Documents/GitHub/penguinteam/')

#---- charge les packages dont on aura besoin
library(R2jags) # pour appeler jags depuis R
library(FSA) # pour calculer les m-array (statistiques reduites)
library(abind) # pour combiner des arrays (generalise r/cbind)

#---- lit les donnees et les explore
dat <- read.table('datadefinitives.txt',header=T)
dat<-
head(dat)
#tail(dat)
names(dat)
dim(dat)
attach(dat)
sum(S.) # nb total d individus
levels(X.COV.site) # modalites var site
levels(X.COV.sexe) # modalites var sexe
levels(X.COV.espece) # modalites var espece
nsites = length(levels(X.COV.site)) # nb sites 
nsp = length(levels(X.COV.espece)) # nb especes
sexes<-levels(X.COV.sexe) # modalites var sexe

#---- met les donnees dans un format utilisable pour l analyse bayesienne
# pre-allocation memoire 
lgr = rep(NA,nsp) # nb site par espece*sexe
lsex = vector()
lsite = vector()
N = NULL # liste des marray par espece x site

uni<-1#unique marray index
for (j in 1:nsp)# boucle sur les especes
  { 
  	k = levels(X.COV.espece)[j] # recupere nom espece courante
  	M = list() # pre-alloue memoire pour les marray de espece courante
  	ind = 1
    tempsp<-subset(dat, X.COV.espece == k )
    for (sex in 1:length(sexes))
      {
        focsex <- sexes[sex]# recupere nom sex courant
        tempsex<-subset(tempsp,X.COV.sexe == focsex)
        for (s in 1:nsites)# boucle sur les sites
          { 
        		i = levels(X.COV.site)[s] # recupere nom site courant
        		temp <- subset(tempsex, X.COV.site == i) # selectionne le jeu de donnees pour site i espece k
        		if(nrow(temp)==0)
              {
          			print(paste(k, focsex, i)) # si ce site ne contient pas d obs pour espece k, affiche l id de ce site
          			next() # et passe au site suivant
        			} 
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
      }#end for (sex in 1:length(sexes))
  	lgr[j] <- length(M) # stocke le nb de sites*sexes de l espece courante
  	N <- c(N,M) # concatene les sites de l espece courante avec ceux des especes precedentes
  }# end for (j in 1:nsp)

length(N) # affiche nb de sites*espece*sex total

# transforme la liste en array qui sera lu par JAGS
Marray=abind(N, along = 3)
dim(Marray) # 25 25  318 // 318 sites avec pour chacun un marray de dim 25x25


save.image(file = "MarrayCreated.RData")#pour pas avoir � refaire tourner la boucle ci-dessus

#############################################################################################
load(file = "MarrayCreated.RData")


# cree un vecteur avec les intervalles qui definissent les sites pour une espece donnee
vecsp = matrix(NA,nrow=max(lgr),ncol=nsp)
somcum = cumsum(lgr)
temp = 1:somcum[1]
vecsp[1:length(temp),1] = temp
for (i in 2:nsp){
	temp <- (somcum[i-1]+1):somcum[i]
	vecsp[1:length(temp),i] <- temp
}
lsex
slopeLsex<-(as.numeric(as.factor(lsex))-2)/2 # -0.5 pour femelles, 0.5 pour males, et 0 pour indetermines (suppose un melange proche de 50%)

# on charge le package qui permet d appeler jags depuis R
library(R2jags)

# on specifie le modele en syntaxe BUGS 
# ou, par exemple, RE(time) signifie effet aleatoire (RE) du temps (time)
sink("cjs-msme.jags")
cat("
model{

#-- effets sur les parametres
for (k in 1:nsp){
	for (t in 1:(n.occasions-1)){
		logit(phi[k,t]) <-  mu.phi +
                         ( epsilon[t] 
                         *alpha[k] )
		for (s in 1:nsites[k]){
			p[k,vecsp[s,k],t] <- mu.p + beta[k] +slopeLsex[vecsp[s,k]]*beta_p_Sex
                          + eta[vecsp[s,k]] # le nb de sites differe par espece - LE CAUCHEMAR
		}
	}
}

for (k in 1:nsp){
	for (s in 1:nsites[k]){
		eta[vecsp[s,k]] ~ dnorm(0,tau.p) # effet aleatoire site sur detection
	}
}

for (k in 1:nsp){	
	alpha[k] ~ dunif(0,1) # effet fixe espece sur survie
	beta[k] ~ dnorm(0,tau.beta) # effet aleatoire espece sur detection
}

for (t in 1:(n.occasions-1)){
	epsilon[t] ~ dunif(0, 1) # effet aleatoire temps sur survie
}

beta_p_Sex ~ dunif(-10,10)#Sex difference in detection probability

#-- priors
#sigma.alpha ~ dunif(0, 5)                
#tau.alpha <- pow(sigma.alpha, -2)
#sigma2.alpha <- pow(sigma.alpha, 2)

sigma.beta ~ dunif(0, 5)                
tau.beta <- pow(sigma.beta, -2)
sigma2.beta <- pow(sigma.beta, 2)

mean.phi ~ dunif(0, 1)             
mu.phi <- log(mean.phi / (1-mean.phi)) 
sigma.phi ~ dunif(0, 5)                
tau.phi <- pow(sigma.phi, -2)
sigma2.phi <- pow(sigma.phi, 2)

mean.p ~ dunif(0, 1)	           
mu.p <- log(mean.p / (1-mean.p)) 
sigma.p ~ dunif(0, 5)                
tau.p <- pow(sigma.p, -2)
sigma2.p <- pow(sigma.p, 2)

#-- vraisemblance multinomiale espece x site
for (k in 1:nsp){ # boucle sur espece
	for (s in 1:nsites[k]){ # boucle sur site
		for (t in 1:(n.occasions-1)){ # boucle sur temps
			marr[t,1:n.occasions,vecsp[s,k]] ~ dmulti(pr[k,vecsp[s,k],t,], rel[t,vecsp[s,k]])
		}
		# Define the cell probabilities of the m-array:
		# Main diagonal
		for (t in 1:(n.occasions-1)){
			q[k,vecsp[s,k],t] <- 1-p[k,vecsp[s,k],t]
			pr[k,vecsp[s,k],t,t] <- phi[k,t]*p[k,vecsp[s,k],t]
			# Above main diagonal
			for (j in (t+1):(n.occasions-1)){
				pr[k,vecsp[s,k],t,j] <- prod(phi[k,t:j])*prod(q[k,vecsp[s,k],t:(j-1)])*p[k,vecsp[s,k],j]
			} #j
			# Below main diagonal
			for (j in 1:(t-1)){
				pr[k,vecsp[s,k],t,j]<-0
			} #j
		} #t
		# Last column: probability of non-recapture
		for (t in 1:(n.occasions-1)){
			pr[k,vecsp[s,k],t,n.occasions] <- 1-sum(pr[k,vecsp[s,k],t,1:(n.occasions-1)])
		} # t
	} # s
} # k
}
",fill = TRUE)
sink()

# met les donnees dans une liste
#nsites[k] ; nb sites par espèce k
#vecsp[s,k] ; pour chaque espèce k, on a les 2 indices qui donnent ou sont les sites, puis la dedans on peut sélectionner le site s
#rel[t,s,k] ; nb individus relâches au temps t, sur site s pour espèce k
jags.data <- list(marr = Marray, n.occasions = dim(Marray[,,1])[2]
                  , rel = apply(Marray,3,rowSums), nsites = lgr, vecsp = vecsp,nsp=nsp, slopeLsex=slopeLsex)

# genere des valeurs initiales
inits <- function(){list(mean.phi = runif(1, 0, 1)
                         , sigma.phi = runif(1, 0, 5)
                         , mean.p = runif(1, 0, 1),sigma.p = runif(1, 0, 5)
                         #,sigma.alpha = runif(1, 0, 5)
                          )}  

# parametres que l on souhaite estimer
parameters <- c("phi","mean.p", "mean.phi"
                , "sigma2.phi"
                ,"sigma2.p"
                #,"sigma2.alpha"
                ,"sigma2.beta",
                "beta_p_Sex"
                )

# specifications MCMC
ni <- 500 # nb iterations
nt <- 1 # nb thining
nb <- 100 # nb burnin
nc <- 1 # nb chains

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
while(inherits(try(tmp <- jags(jags.data, inits, parameters, "cjs-msme.jags"
                               , n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb
                               , working.directory = getwd()), silent=TRUE), "try-error")) {}
temp.result <- tmp 

# affiche les resultats
print(temp.result, digits = 3)


