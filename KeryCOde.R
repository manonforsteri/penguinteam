# Define parameter values
n.occasions <- 6                   # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.65, n.occasions-1)
p <- rep(0.4, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))
# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)){
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i]==n.occasions) next
    for (t in (mark.occ[i]+1):n.occasions){
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i,t-1])
      if (sur==0) break  	# If dead, move to next individual 
      # Bernoulli trial: is individual recaptured? 
      rp <- rbinom(1, 1, P[i,t-1])
      if (rp==1) CH[i,t] <- 1
    } #t
  } #i
  return(CH)
}

# Execute function
CH <- simul.cjs(PHI, P, marked)

# 7.10. Fitting the CJS to data in the m-array format: the multinomial likelihood
# 7.10.1. Introduction
# Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

# 7.10.2. Time-dependent models
# Specify model in BUGS language
sink("cjs-mnl.jags")
cat("
    model {
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
    phi[t] ~ dunif(0, 1)         # Priors for survival
    p[t] ~ dunif(0, 1)           # Priors for recapture
    }
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
    marr[t,1:n.occasions] ~ dmulti(pr[t, ], rel[t])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
    q[t] <- 1-p[t]                # Probability of non-recapture
    pr[t,t] <- phi[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
    pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr[t,j] <- 0
    } #j
    } #t
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
    pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } #t
    
    }
    ",fill = TRUE)
sink()

# Create the m-array from the capture-histories
marr <- marray(CH)

# Bundle data
jags.data <- list(marr = marr, n.occasions = dim(marr)[2], rel = rowSums(marr))

# Initial values
inits <- function(){list(phi = runif(dim(marr)[2]-1, 0, 1), p = runif(dim(marr)[2]-1, 0, 1))}  

# Parameters monitored
parameters <- c("phi", "p")

# MCMC settings
ni <- 10000
nt <- 3
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs <- jags(jags.data, inits, parameters, "cjs-mnl.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

print(cjs, digits = 3) 
