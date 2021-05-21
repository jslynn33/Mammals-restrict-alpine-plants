### Q1 Analysis of damage on Festuca brachyphylla (FB)


## Packages
library(R2jags);library(car);library(tidyverse);library(lubridate)

# load in data
dat <- read.csv("data_lives_here/MasterHerbExpData.csv")
head(dat)

# make above ground (cage) and below ground (gopher) treatments factors
dat$cage <- as.factor(dat$cage)
dat$gop <- as.factor(dat$gopher)

# create column of damage presences and absence
dat$presherb <- 1*(dat$damage>0) 

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

#data for each species-- here, FEBR or Festuca brachyphylla
fbdat <- dat[dat$species=="FEBR",]

# set up model for herbivory presence
presdat <- fbdat[!is.na(fbdat$presherb),]

# set up random effects
ind <- factor(presdat$ind)
sit <- factor(presdat$sit)

# Dependent variables
pres <- presdat$presherb
N <- as.numeric(length(pres))

# independent variables
doy <- presdat$doy
presmat <- model.matrix(~site*cage*gopher, data=presdat)
head(presmat)

# set up data and parameters to use
jags.data <- list("pres", "N", "doy", "presmat", "sit", "ind")
jags.param <- c("a", "bt", "prec1", "prec2","sig1", "sig2","rss", "rss.new")

# the model
dampres.mod <- function(){
  for (j in 1:9){d[j]~dnorm(0, prec1)}# site random effects
  for (j in 1:177){en[j]~dnorm(0, prec2)} # individual random effect
  for (i in 1:N){
    pres[i]~dbern(mu[i])
    logit(mu[i]) <- inprod(a, presmat[i,]) + bt*doy[i]+d[sit[i]]+en[ind[i]]
    
    res[i] <- (pres[i]-mu[i])^2 # for posterior predictive checks
    pres.new[i]~dbern(mu[i])
    res.new[i] <- (pres.new[i]-mu[i])^2
  }
  #priors
  for(j in 1:12){a[j]~ dnorm(0,1.0E-6)}
  bt~dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001) 
  prec2~dgamma(0.001,0.001) 
  sig1 <- 1/sqrt(prec1)
  sig2 <- 1/sqrt(prec2)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
presdam <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                         n.iter=50000 ,model.file=dampres.mod, n.thin=5, n.chains=3)

presdam # results

#posterior predictive check
presdam.draws <- presdam$BUGSoutput$sims.list

plot(presdam.draws$rss,presdam.draws$rss.new, main="Herbivory Presence",
     xlab="SSQ observed",ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(presdam.draws$rss>presdam.draws$rss.new) 


### all damage classes

# subset data to only data with greater than zero damage
damdat <- fbdat[fbdat$damage>0,]
damdat <- damdat[!is.na(damdat$damage),]
damdat <- damdat[!(damdat$DamCat=="gopher"),] # exclude gopher disturbed (only two cases of this)

# set up random effects
ind <- factor(damdat$ind)
sit <- factor(damdat$sit)

# logit transform damage as dependent variable
dam <- logit(damdat$damage/100, adjust=0.01)
N <- as.numeric(length(dam))

# independent variables
doy <- damdat$doy
dammat <- model.matrix(~site*cage*gop, data=damdat)

# set up data and parameters to use
jags.data <- list("dam", "N", "dammat", "sit", "ind", "doy")
jags.param <- c("a", "bt", "prec1", "prec2", "prec3", "sig1", "sig2", "sig3",
                "rss", "rss.new")

# the model
damall.mod <- function(){
  for (j in 1:9){d[j]~dnorm(0, prec1)}# site random effects
  for (j in 1:79){en[j]~dnorm(0, prec2)} # individual random effect
  for (i in 1:N){
    dam[i]~ dnorm(mu[i], prec3) 
    mu[i] <- inprod(a, dammat[i,])+bt*doy[i] +d[sit[i]]+en[ind[i]]
    
    res[i] <- (dam[i]-mu[i])^2 # posterior predictive checks
    dam.new[i]~ dnorm(mu[i], prec3)
    res.new[i] <- (dam.new[i]-mu[i])^2
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
damall <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=50000 ,model.file=damall.mod, n.thin=5, n.chains=3)

damall # the results

#posterior predictive check
damall.draws <- damall$BUGSoutput$sims.list

plot(damall.draws$rss,damall.draws$rss.new, main="Herbivory all damage",
     xlab="SSQ observed",ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(damall.draws$rss>damall.draws$rss.new) 


#### model for just insect damage
insdamdat <- damdat[damdat$DamCat=="insect",]

# set up random effects
ind <- factor(insdamdat$ind)
sit <- factor(insdamdat$sit)

# dependent variable
dam <- logit(insdamdat$damage/100, adjust=0.01)
N <- as.numeric(length(dam))

# independent variables
doy <- insdamdat$doy
insdammat <- model.matrix(~site*cage*gopher, data=insdamdat)

# set up data and parameters to use
jags.data <- list("dam", "N", "insdammat", "sit", "ind", "doy")
jags.param <- c("a", "bt", "prec1", "prec2", "prec3","sig1", "sig2", "sig3",
                "rss", "rss.new")

# the model
damins.mod <- function(){
  for (j in 1:9){d[j]~dnorm(0, prec1)}# site random effects
  for (j in 1:62){en[j]~dnorm(0, prec2)} # individual random effects
  for (i in 1:N){
    dam[i]~ dnorm(mu[i], prec3) 
    mu[i] <- inprod(a, insdammat[i,])+bt*doy[i] +d[sit[i]]+en[ind[i]]
    
    res[i] <- (dam[i]-mu[i])^2 # posterior predictive check
    dam.new[i]~ dnorm(mu[i], prec3)
    res.new[i] <- (dam.new[i]-mu[i])^2
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
damins <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=50000 ,model.file=damins.mod, n.thin=5, n.chains=3)

damins # results

#posterior predictive check
damins.draws <- damins$BUGSoutput$sims.list

plot(damins.draws$rss,damins.draws$rss.new, main="Herbivory insect damage",
     xlab="SSQ observed",ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(damins.draws$rss>damins.draws$rss.new) 


#### model for just mammal damage

mamdamdat <- damdat[damdat$DamCat=="mammal",]

# set up random effects
ind <- factor(mamdamdat$ind)
sit <- factor(mamdamdat$sit)

# dependent variable
dam <- logit(mamdamdat$damage/100, adjust=0.01)
N <- as.numeric(length(dam))

# independent variables
doy <- mamdamdat$doy
mamdammat <- model.matrix(~site, data=mamdamdat)
head(mamdammat)

# set up data and parameters to use
jags.data <- list("dam", "N", "mamdammat", "sit", "ind", "doy")
jags.param <- c("a", "bt", "prec1", "prec2", "prec3","sig1", "sig2", 
                "sig3","rss", "rss.new")

# the model 
dammam.mod <- function(){
  for (j in 1:7){d[j]~dnorm(0, prec1)}# site random effects
  for (j in 1:27){en[j]~dnorm(0, prec2)} # individual random effect
  for (i in 1:N){
    dam[i]~ dnorm(mu[i], prec3) 
    mu[i] <- inprod(a, mamdammat[i,])+bt*doy[i] +d[sit[i]]+en[ind[i]]
    
    res[i] <- (dam[i]-mu[i])^2 # poseterior predictive check
    dam.new[i]~ dnorm(mu[i], prec3)
    res.new[i] <- (dam.new[i]-mu[i])^2
  }
  #priors
  for(j in 1:3){a[j]~ dnorm(0,1.0E-6)}
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
dammam <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=50000 ,model.file=dammam.mod, n.thin=5, n.chains=3)

dammam # resutls

#posterior predictive check
dammam.draws <- dammam$BUGSoutput$sims.list

plot(dammam.draws$rss,dammam.draws$rss.new, main="Herbivory mammal damage",
     xlab="SSQ observed",ylab="SSQ simulated")
abline(0,1)

# Bayesian p-values
mean(dammam.draws$rss>dammam.draws$rss.new) 

