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

#avec jeu de données partiel:
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

