library(rjags)

setwd("C:/Users/Manon Ghislain L/Documents/Scripts/MCMC Titimounet")

CH<-data.matrix(read.table(file="PHYCOL.txt",header=F,colClasses="double"))
