### Q2 Biomass, flowering, and survival for Elymus scribneri

## Packages
library(R2jags);library(car);library(tidyverse);library(lubridate)

# load in data
dat <- read.csv("data_lives_here/MasterHerbExpData.csv")
head(dat)

# make above ground (cage) and below ground (gopher) treatments factors
dat$cage <- as.factor(dat$cage)
dat$gop <- as.factor(dat$gopher)

# create a group for individuals
dat$ind <- as.factor(paste(dat$peak, dat$site, dat$plant ,dat$species, sep=""))

# create a group for site
dat$sit <- as.factor(paste(dat$peak, dat$site, sep=""))

# get day of the year for repeat measures analysis
dat$fixdat <- mdy(dat$date)
dat$doy <- yday(dat$fixdat)

# so that everything is estimated from the  minimum date in the data set
dat$doy <- dat$doy-min(dat$doy)

# remove the pre-transplant ("PreTrans") data
dat <- dat[!(dat$samp_period=="PreTrans"),]

#data for each species-- here, ELSC or Elymus scribneri
esdat <- dat[dat$species=="ELSC",]

# use allometric equation to get predicted biomass
# predictions can be found allometric equation script
esdat$allolivebio <- 0.0251527*esdat$tiller+0.0014635*esdat$tiller^2
 

### biomass models for Elymus scribneri

biodat <- esdat[!is.na(esdat$allolivebio),]

# random effects
sit <- factor(biodat$sit)
ind <- factor(biodat$ind)

# dependent variable
bio <- log(biodat$allolivebio)
N <- as.numeric(length(bio))

# independent variables
doy <- biodat$doy
biomat <- model.matrix(~site*cage*gopher, biodat)
head(biomat)

# set up data and parameters to use
jags.data <- list("bio", "N", "biomat", "sit", "ind", "doy")
jags.param <- c("a","bt", "prec1","prec2", "prec3","sig1","sig2", "sig3",
                "rss", "rss.new")

# the model
biomass.mod <- function(){
  for (j in 1:9){d[j]~dnorm(0, prec1)}# site random effects
  for (j in 1:177){en[j]~dnorm(0, prec2)} # individual random effect
  for (i in 1:N){
    bio[i]~ dnorm(mu[i], prec3) 
    mu[i] <- inprod(a, biomat[i,])+bt*doy[i]+d[sit[i]]+en[ind[i]]
    
    res[i] <- (bio[i]-mu[i])^2 # posterior predictive cheeck
    bio.new[i]~ dnorm(mu[i], prec3)
    res.new[i] <- (bio.new[i]-mu[i])^2
  }
  #priors
  for(j in 1:12){a[j]~ dnorm(0,1.0E-6)}
  bt~dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001) 
  prec2~dgamma(0.001,0.001) 
  prec3~dgamma(0.001,0.001) 
  sig1 <- 1/sqrt(prec1)
  sig2 <- 1/sqrt(prec2)
  sig3 <- 1/sqrt(prec3)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
biojags <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                         n.iter=50000 ,model.file=biomass.mod, n.thin=5, n.chains=3)

biojags # results

#posterior predictive check
bio.draws <- biojags$BUGSoutput$sims.list

plot(bio.draws$rss,bio.draws$rss.new, main="Biomass",
     xlab="SSQ observed",ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(bio.draws$rss>bio.draws$rss.new) 


### flower number models for ES

# only include 17 and 18 when plants flowered.  
# Flowers from a year were summed to "a" or august samp_period
flowdat <- esdat %>% filter(samp_period %in% c("a17", "a18")) 
flowdat <- flowdat[!is.na(flowdat$flower),] # take out NAs
flowdat <- flowdat[flowdat$flower>0,] # only analyze individuals that flowered

# group random effects
sit <- factor(flowdat$sit)
ind <- factor(flowdat$ind)

# dependent variables
flow <- flowdat$flower
N <- as.numeric(length(flow))

# independent variables
sp <- factor(flowdat$samp_period)
flowmat <- model.matrix(~site*cage*gopher, flowdat)

# set up data and parameters to use
jags.data <- list("flow", "N", "flowmat", "sit", "ind", "sp")
jags.param <- c("a", "bt", "prec1", "prec2","sig1", "sig2","rss", "rss.new", "r")

# the model
flower.mod <- function(){
  for (j in 1:9){d[j]~dnorm(0, prec1)}# site random effects
  for (j in 1:107){en[j]~dnorm(0, prec2)} # individual
  for (i in 1:N){
    flow[i]~ dnegbin(p[i], r) 
    log(mu[i]) <- inprod(a, flowmat[i,])+bt*sp[i]+d[sit[i]]+en[ind[i]]
    p[i] <- r/(r + mu[i])  
    
    res[i] <- (flow[i]-p[i])^2 # posterior predictive check
    flow.new[i]~ dnegbin(p[i], r)
    res.new[i] <- (flow.new[i]-p[i])^2
  }
  #priors
  r~ dunif(0,50)
  for(j in 1:12){a[j]~ dnorm(0,1.0E-6)}
  bt ~ dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001) 
  prec2~dgamma(0.001,0.001) 
  sig1 <- 1/sqrt(prec1)
  sig2 <- 1/sqrt(prec2)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
flowjags <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                          n.iter=500000 ,model.file=flower.mod, n.thin=5, n.chains=3)

flowjags # results

#posterior predictive check
flow.draws <- flowjags$BUGSoutput$sims.list

plot(flow.draws$rss,flow.draws$rss.new, main="Flower",
     xlab="SSQ observed",ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(flow.draws$rss>flow.draws$rss.new) 


### survival analysis
# survival asses until end of experiment- use biomass data
datsurv <- read.csv("data_lives_here/MamExpBioSurvData.csv")

# pull out Poa alpina data
essurv <- datsurv[datsurv$species=="ELSC",]
head(essurv)

# dependent variable
surv <- essurv$surv
N <- as.numeric(length(surv))

# independent variables
survmat <- model.matrix(~site*cage*gopher, essurv)

# set up data and parameters to use
jags.data <- list("surv", "N", "survmat")
jags.param <- c("a",  "rss", "rss.new", "r")

# the model
survival.mod <- function(){
  for (i in 1:N){
    surv[i]~ dbern(mu[i]) 
    logit(mu[i]) <- inprod(a, survmat[i,])
    
    res[i] <- (surv[i]-mu[i])^2 # posterior predictive check
    surv.new[i]~ dbern(mu[i])
    res.new[i] <- (surv.new[i]-mu[i])^2
  }
  #priors
  for(j in 1:12){a[j]~ dnorm(0,1.0E-6)}
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
survjags <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                          n.iter=50000 ,model.file=survival.mod, n.thin=5, n.chains=3)

survjags # results

#posterior predictive check
surv.draws <- survjags$BUGSoutput$sims.list

plot(surv.draws$rss,surv.draws$rss.new, main="Survival",
     xlab="SSQ observed",ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(surv.draws$rss>surv.draws$rss.new) 

