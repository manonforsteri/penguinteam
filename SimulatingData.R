setwd(dir = "C:/Users/Thimothee Admin/Documents/GitHub/penguinteam/")

###### Simulation model phi(sp+t) p(.) #####
n.occasions <- 6                   # Number of capture occasions
effect_Occasions<-rnorm(n = n.occasions-1,mean = 0,sd = 0.5)#time effect
meanphi<-0.5
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
p <- rep(0.4, n.occasions-1)#constant recapture probability

# Define matrices with survival and recapture probabilities
PHI <- matrix(1/(1+exp(-(meanphi+effect_Occasions))), ncol = n.occasions-1, nrow = sum(marked))
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
    logit(phi[i,t]) <- mu.phi
    logit(p[i,t]) <- mu.p
    } #t
    } #i
    
    meanphi ~ dunif(0, 1)              # Priors for group-specific survival

    mu.phi<- log(phi.mean/(1-phi.mean))
    phi.mean ~ dunif(0,1)

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
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), phi.mean = runif(1, 0, 1), p.mean = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("phi.mean", "p.mean")

# MCMC settings
ni <- 1000
nt <- 3
nb <- 500
nc <- 3

# Call JAGS from R (BRT 2 min)
cjs.group <- jags(jags.data, inits, parameters, "cjs-phi(sp+t)p(.).jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.group, digits = 3)
