library(rwc)
library(rgr)
source("rwc.fcns.r")


#############################################
## Load covariate list and Fst Matrix (D)

matrix.list<-readRDS (file="covariates.rds")
names(matrix.list)<-c("intercept", "main.cov", "maindist.cov", "edgelength.cov", "direction.cov","barrier.cov","dry.cov", "order.cov", "area.cov",
                       "link.cov","elev.cov","slope.cov","roadXing.cov")

D<-readRDS(file="D.rds")


##########################################
## Model fitting using MCMC##
## Aditional details can be found in supporting documentation for the rwc package (Hanks 2018)
## model is D~GenWish(df, Sigma), Sigma from Hanks (2017) or Hanks and Hooten (2013)
##########################################

n.nodes=c(103,72:74,76:102,104,105) # Create a vector containing the column numbers of observed nodes within the larger covariate matrices

fit=mcmc.wish(D, ## genetic distance matrix 
               df=12, ## df = number of microsats
               matrix.list, ## list of covariate matrices 
               obs.idx=n.nodes, ## numeric vector of the sampling locations.
               ## the length of obs.idx = nrow(D), and the i-th entry of
               ## obs.idx is the index of the node in the stream network
               ## corresponding to the i-th row of D
               model="SAR", ## use "SAR" for Hanks (2017) method, which allows for asymmetric flow
               beta.start=rep(0.15,length(matrix.list)), ## starting values for gene flow parameters. The model may be sensitive to the starting value, in which case chains will break. 
               n.mcmc=50000,## number of MCMC iterations
               adapt.max=25000, #specifying the last iteration at which the covariance matrix of the proposal distribution will be adapted; usually half the value of n.mcmc
               adapt.int=500, #Interval at which the covariance matrix of the proposal distribution is adapted
               print.iter=TRUE,## set to FALSE if you don't want to see output at each iteration
               output.trace.plot=TRUE) ## set to FALSE if you don't want the output to be saved to a pdf in the working directoy. 

## plot trace plots for betas, calculate DIC, and obtain parameter estimates
matplot(fit$beta,type="l")
str(fit)
apply(fit$beta,2,mean)
apply(fit$beta,2,quantile,c(0.025,0.975))








