library(R2jags)
library(Hmisc)
Sys.setenv(MAKEFLAGS = "-j8")  

sink("modelB")
cat("
    
    model{
    
    #Model
    for(i in 1:I){
    for(t in (f[i]):(Kmax-1)) {
    logit(phi[i,t])<-  phi.temps[t] + phi.trans[x[i,t]] + phi.sex[sex[i]] + phi.species[species[i]] + random.phi.site[site[i]]
    lp[i,t]<-  p.temps[t] + p.sex[sex[i]] + random.p.site[site[i]] + b1*covhetero[i] + p.species[species[i]]
    p[i,t]<-   A[i,t]/(1+exp(-lp[i,t]))                         #inverse logit
    
    } #t
    }#i
    
    #####Priors#####
    
    
    #for survival parameters
    for(t in 1:(Kmax-1)){
    phi.temps[t] ~ dnorm(0,0.01)I(-10,10)
    }
    
    phi.trans[1]<-0
    phi.trans[2] ~ dnorm(0,0.01)I(-10,10)
    
    phi.sex[1] <- 0
    phi.sex[2] ~ dnorm(0,0.01)I(-10,10)
    phi.sex[3] ~ dnorm(0,0.01)I(-10,10)
    
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
    
    p.sex[1]<-0
    p.sex[2] ~ dnorm(0,0.01)I(-10,10)
    p.sex[3] ~ dnorm(0,0.01)I(-10,10)
    
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
    p.temps[t] ~ dnorm(0,0.01)I(-10,10)
    }
    
    
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
    phi.temps=rnorm(Kmax-1, mean = 0, sd = 0.01),
    phi.trans=c(NA,rnorm(1)),
    phi.sex=c(NA,rnorm(nsex-1, mean = 0, sd = 0.01)),
    phi.species=c(NA,rnorm(nspecies-1)),
    p.sex=c(NA,rnorm(nsex-1, mean = 0, sd = 0.01)),
    p.temps=rnorm(Kmax-1, mean = 0, sd = 0.01),
    b1=rnorm(1, mean = 0, sd = 0.01),
    p.species=c(NA,rnorm(nspecies-1, mean = 0, sd = 0.01))
  )}

modelB<- jags(data=data, inits=inits, parameters=parameters, "modelB", n.thin=10, n.chains=3, n.burnin=500, n.iter=1000)
# modB<-autojags(modB) 
save(modelB, file="C:/Users/MGGCBRUNOY/Desktop/nouveau JAGS/modelB")


