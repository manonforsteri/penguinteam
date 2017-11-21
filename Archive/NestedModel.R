# data simulation

nl1 <- 40
level1 <- rnorm(nl1, mean = 0, sd = 1.5)

nl2 <- 40

reff <- matrix(data = NA, nrow = nl1, ncol = nl2)
for (i in 1:nl1)
{
  reff[i,] <- rnorm(nl2, mean = level1[i], sd =0.5)
}


ndata <- nl1*nl2*30

sl1 <- sample(x=1:nl1, size = ndata, replace = TRUE)
sl2 <- sample(x=1:nl2, size = ndata, replace = TRUE)

y = vector(length=ndata)
for (i in 1:ndata)
  {
    y[i] = 1+reff[sl1[i],sl2[i]]+rnorm(n = 1,0,0.1)
  }
ydf = data.frame(y=y, sl1 = sl1, sl2=(sl1-1)*nl2+sl2)

library(lme4)
lmer(formula = y ~ 1 + (1|sl1/sl2), data=ydf)
lmer(formula = y ~ 1 + (1|sl1)+(1|sl2), data=ydf)
#it's the same to nest or not nest explicitly

setwd("/home/timothee/Documents/GitHub/penguinteam/")
#############################
#### Now in JAGS
sink("model1nest")
cat("
    
    model{
#prior  
  int ~ dnorm(0,10)
  tau.resid <- pow(sigma.resid, -2)
  sigma.resid ~ dunif(0,10)

  for (j in 1:nl1)
    {
      random.l1[j] ~ dnorm(0,tau.sl1)
    }
  tau.sl1 <- pow(sigma.sl1, -2)
  sigma.sl1 ~ dunif(0,10)

  for (h in 1:nl2)
    {
      random.l2[h] ~ dnorm(0,tau.sl2)
    }
  tau.sl2 <- pow(sigma.sl2, -2)
  sigma.sl2 ~ dunif(0,10)

#likelihood
  for (i in 1:ndata)
    {
      y[i] ~ dnorm(ye[i], tau.resid)
      ye[i] <- int + random.l1[sl1[i]] + random.l2[sl2[i]]
    }

    }
",fill=TRUE)
sink()

library(R2jags)
parameters<-c("int", "sigma.resid", "sigma.sl1", "sigma.sl2")
data = list(y = ydf$y, ndata=nrow(ydf), sl1 = ydf$sl1, nl1 = length(unique(ydf$sl1)),  sl2=ydf$sl2,nl2 = length(unique(ydf$sl2)))
inits = function(){list(int=rnorm(n = 1,mean = 0,sd = 10),
sigma.resid=runif(1,0,10), sigma.sl1=dunif(1,0,10), sigma.sl2=dunif(1,0,10))}

jags1<- jags(data=data, inits=inits, parameters=parameters, "model1nest", 
                    n.thin=10, n.chains=3, n.burnin=100, n.iter=1100)
print(jags1)
