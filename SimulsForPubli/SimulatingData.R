setwd(dir = "~/Documents/GitHub/penguinteam/")
library(R2jags)

###### Simulation model phi(sp+t) p(.) #####
n.occasions <- 10                   # Number of capture occasions
meanphi<-0#logit intercept survival (i.e 0.5 on survival scale)
effect_Occasions<-c(0,rnorm(n = n.occasions-2,mean = 0,sd = 0.5))#time effect

n.sp<-7
effect_Sp<-c(0,rnorm(n = n.sp-1,mean=0,sd=1))


marked <- round(100*(1:n.sp)^-2)   # Annual number of newly marked individuals
p <- rep(0.4, n.occasions-1)#constant recapture probability

# Define matrices with survival and recapture probabilities
PHI <- matrix(NA, ncol = n.occasions-1, nrow = sum(marked))
sp_ind<-sample(x = 1:n.sp,size = sum(marked),replace = TRUE)#to which species does this individual belong?
for (i in 1:sum(marked))
  {
    for (t in 1:(n.occasions-1))
      {
        PHI[i,t]<-1/(1+exp(-( meanphi + effect_Occasions[t] + effect_Sp[sp_ind[i]] ) ) )# additive effects only
      }
  }

P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])  # Define a vector with the occasion of marking
  
  for (i in 1:sum(marked)){# Fill the CH matrix
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break  	# If dead, move to next individual 
      
      rp <- rbinom(1, 1, P[i,t-1])# Bernoulli trial: is individual recaptured? 
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

# Execute function
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}


cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}


#try in JAGS
sink("cjs-phi(sp+t)p(.).jags")
cat("
    model {
    # Priors and constraints
    for (i in 1:nind){
      for (t in f[i]:(n.occasions-1)){
        logit(phi[i,t]) <- mu.phi + beta.t[t] + beta.sp[ sp_ind[i] ]
        logit(p[i,t]) <- mu.p
      } #t
    } #i
    
    meanphi ~ dunif(0, 1)              # Priors for group-specific survival

    mu.phi<- log(phi.mean/(1-phi.mean))
    phi.mean ~ dunif(0,1)

    beta.t[1]<-0
    for (t in 2:(n.occasions-1)){
      beta.t[t] ~ dnorm(0,0.01)
    }

    beta.sp[1]<-0
    for (sp in 2:n.sp){
      beta.sp[sp] ~ dnorm(0,0.01)
    }

    mu.p <- log(p.mean/(1-p.mean))
    p.mean ~ dunif(0,1)
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
    } #i
    }
    ",fill = TRUE)
sink()


sink("cjs-phi(spXt)p(.).jags")
cat("
    model {
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu.phi + beta.t[t] + beta.sp[ sp_ind[i] ] + beta.tXsp[ t,sp_ind[i] ]
    logit(p[i,t]) <- mu.p
    } #t
    } #i
    
    meanphi ~ dunif(0, 1)              # Priors for group-specific survival
    
    mu.phi<- log(phi.mean/(1-phi.mean))
    phi.mean ~ dunif(0,1)
    
    for (t in 1:(n.occasions-1)){
    beta.t[t] ~ dnorm(0,tau.t)
    }
    tau.t <- pow(sigma.t, -2)
    sigma.t ~ dunif(0,10)
    
    for (sp in 1:n.sp){
    beta.sp[sp] ~ dnorm(0,tau.sp)
    for (t in 1:(n.occasions-1)){
    beta.tXsp[t,sp] ~ dnorm(0,tau.tXsp)
    }
    }
    tau.sp <- pow(sigma.sp, -2)
    sigma.sp ~ dunif(0,10)
    tau.tXsp <- pow(sigma.tXsp, -2)
    sigma.tXsp ~ dunif(0,10)

    mu.p <- log(p.mean/(1-p.mean))
    p.mean ~ dunif(0,1)
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

sink("cjs-phi(spXt_no_t)p(.).jags")
cat("
    model {
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
    logit(phi[i,t]) <- mu.phi + beta.sp[ sp_ind[i] ] + beta.tXsp[ t,sp_ind[i] ]
    logit(p[i,t]) <- mu.p
    } #t
    } #i
    
    meanphi ~ dunif(0, 1)              # Priors for group-specific survival
    
    mu.phi<- log(phi.mean/(1-phi.mean))
    phi.mean ~ dunif(0,1)

    for (sp in 1:n.sp){
    beta.sp[sp] ~ dnorm(0,tau.sp)
    for (t in 1:(n.occasions-1)){
    beta.tXsp[t,sp] ~ dnorm(0,tau.tXsp)
    }
    }
    tau.sp <- pow(sigma.sp, -2)
    sigma.sp ~ dunif(0,10)
    tau.tXsp <- pow(sigma.tXsp, -2)
    sigma.tXsp ~ dunif(0,10)
    
    mu.p <- log(p.mean/(1-p.mean))
    p.mean ~ dunif(0,1)
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
    # State process
    z[i,t] ~ dbern(mu1[i,t])
    mu1[i,t] <- phi[i,t-1] * z[i,t-1]
    # Observation process
    y[i,t] ~ dbern(mu2[i,t])
    mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), sp_ind=sp_ind, n.sp=length(unique(sp_ind)) )

#### Initial values####
#interaction init

beta.tXsp_ini<-function()
{beta.tXsp<-matrix(NA,nrow = n.occasions-1,ncol = n.sp)
for (i in 2:(n.occasions-1))
  {
    for (j in 2:n.sp)
      {
        beta.tXsp[i,j]<-rnorm(n = 1,mean = 0,sd = 1)
      }
  }
return(beta.tXsp)
}

inits <- function(){list(z = cjs.init.z(CH, f), phi.mean = runif(1, 0, 1), p.mean = runif(1, 0, 1),
                         beta.t=c(NA,rnorm(n = n.occasions-2,mean = 0,sd=1 ) ),
                         beta.sp=c(NA,rnorm(n = n.sp-1,mean = 0,sd=1 ) ),
                         beta.tXsp= beta.tXsp_ini() )}  

# Parameters monitored
parameters <- c("phi.mean", "p.mean","beta.t","beta.sp","beta.tXsp", "sigma.sp", "sigma.t", "sigma.tXsp")

# MCMC settings
ni <- 11000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R 
#cjs.additif <- jags(jags.data, inits, parameters, "cjs-phi(sp+t)p(.).jags",
#                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
#print(cjs.additif, digits = 3)

cjs.interactionS <- jags(jags.data, inits, parameters, "cjs-phi(spXt)p(.).jags",
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(cjs.interactionS, digits = 3)

# Summarize posteriors

predictions <- cjs.interactionS$BUGSoutput$mean$beta.tXsp+ +matrix(cjs.interactionS$BUGSoutput$mean$beta.t, nrow = 9, ncol = n.sp, byrow = FALSE)

matrix(cjs.interactionS$BUGSoutput$mean$beta.sp, nrow = 9, ncol = n.sp, byrow = TRUE)

plot(predictions[,1], type="l", ylim=c(-0.5,0.5), xlab="year", main="With year RE")
for(i in 2:5)
{
  lines(predictions[,i])
}

cjs.interaction$BUGSoutput$mean$sigma.t/(cjs.interaction$BUGSoutput$mean$sigma.tXsp+cjs.interaction$BUGSoutput$mean$sigma.t)


parameters <- c("phi.mean", "p.mean","beta.sp","beta.tXsp", "sigma.sp", "sigma.tXsp")
cjs.interaction_no_t <- jags(jags.data, inits, parameters, "cjs-phi(spXt_no_t)p(.).jags",
                             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(cjs.interaction_no_t, digits = 3)

predictions2 <- cjs.interaction_no_t$BUGSoutput$mean$beta.tXsp
matrix(cjs.interaction_no_t$BUGSoutput$mean$beta.sp, nrow = 9, ncol = n.sp, byrow = TRUE)

plot(predictions2[,1], type="l", ylim=c(-0.5,0.5), main = "Without year RE", xlab="year")
for (i in 2:n.sp)
{
  lines(predictions2[,i])
}

var(apply(predictions2, MARGIN = 1, mean))/var(as.vector(predictions2)-apply(predictions2, 2, mean))


var(apply(predictions, MARGIN = 1, mean))/var(as.vector(predictions)-apply(predictions, 2, mean))

var(apply(predictions2, 2, mean))
