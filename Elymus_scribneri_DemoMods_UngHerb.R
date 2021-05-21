### Demographic analyses Poa alpina


# Packages
library(R2jags); library(car);library(tidyverse)

# load in demographic data
demodat <- read.csv("data_lives_here/MasterDemoDataset.csv")
head(demodat)

# make year and plot number a factor
demodat$year <- as.factor(demodat$year)
demodat$plotnumb <- as.factor(demodat$plotnumb)

# pull out Poa alpina (PA) data
esdat <- demodat[demodat$species=="ELSC",]

### vital rate calculations


# Parameterize the growth model for ES

# remove NA's in ESdata
esdat2 <- esdat[!is.na(esdat$tillernumb1 & esdat$tillernumb2),]


tillernumb1 <- esdat2$tillernumb1 # tiller number at in year t
tillernumb2 <- esdat2$tillernumb2 # tiller number at in year t+1

year <- factor(esdat2$year) # for random effect
N <- as.numeric(length(esdat2$tillernumb1))# get length of dataset

# set up data and parameters to use
jags.data <- list("N","tillernumb1","tillernumb2", "year")
jags.param <- c( "a","bESg","r","rss", "rss.new", "prec1")

# the model
esgrowthmod <- function(){
  for(j in 1:3){aar[j]~dnorm(0, prec1)} # year random variables
  for (i in 1:N){
    tillernumb2[i]~dnegbin(p[i],r)
    log(mu[i]) <- a+ bESg*log(tillernumb1[i])+aar[year[i]]
    p[i] <- r/(r+mu[i])
    
    res[i] <- (tillernumb2[i]-mu[i])^2 # posterior predictive check
    tiller.new[i]~dnegbin(p[i],r)
    res.new[i] <- (tiller.new[i]-mu[i])^2
  }
  # priors
  r~ dunif(0,50)
  a~ dnorm(0,1.0E-6)
  bESg~dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  # derived parameters
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
esgrowth<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                n.iter=100000 ,model.file=esgrowthmod, n.thin=5, n.chains=3)

esgrowth # results

# post predictive check
esgrowth.paramlist <- esgrowth$BUGSoutput$sims.list

plot(esgrowth.paramlist$rss,esgrowth.paramlist$rss.new, main="A. growth",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(esgrowth.paramlist$rss>esgrowth.paramlist$rss.new)

## take the 95% cred. intevals of each
bESg <- sort(esgrowth.paramlist$bESg) [750:29250]
aESg <- sort(esgrowth.paramlist$a) [750:29250]
rESg <- sort(esgrowth.paramlist$r) [750:29250]

esgrowth_post <- as.data.frame(cbind(bESg,aESg,rESg))
write.csv(esgrowth_post, "esgrowth_post.csv")


### survival vital rates for ES

esdat3 <- esdat[!is.na(esdat$survto2),] #subset out individuals still dead

# independent variables
tillernumb1 <- esdat3$tillernumb1
year <- factor(esdat3$year) # year random effect

# dependent variables 
survto2 <- esdat3$survto2
N <- as.numeric(length(esdat3$survto2))

# set up data and parameters to use
jags.data <- list("N","tillernumb1","survto2", "year")
jags.param <- c( "a","bESs", "rss", "rss.new", "prec1")

# the model 
essurvmod <- function(){
  for(j in 1:3){y[j]~dnorm(0, prec1)} # year random effect
  for (i in 1:N){
    survto2[i]~dbern(mu[i])
    logit(mu[i]) <- a + bESs*log(tillernumb1[i])+y[year[i]]
    
    res[i] <- (survto2[i]-mu[i])^2 # posterior predictive check
    surv.new[i]~dbern(mu[i])
    res.new[i] <- (surv.new[i]-mu[i])^2
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bESs~ dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
essurv<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
              n.iter=200000 ,model.file=essurvmod, n.thin=5, n.chains=3)

essurv # results

#posterior predictive check
essurv.paramlist <- essurv$BUGSoutput$sims.list

plot(essurv.paramlist$rss,essurv.paramlist$rss.new, main="B. survival",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(essurv.paramlist$rss>essurv.paramlist$rss.new)

# get 95 % cred of vital rate params
bESs <- sort(essurv.paramlist$bESs) [1500:58500]
aESs <- sort(essurv.paramlist$a) [1500:58500]

essurv_post <- as.data.frame(cbind(bESs,aESs))
write.csv(essurv_post, "essurv_post.csv")


## probability of flowering vital rates

# create a 0,1 column for flowering
demoflow <- demodat[!is.na(demodat$flowernumb),]
did_flower <- 0
for(i in 1:length(demoflow$flowernumb)){
  did_flower[i]<- if(demoflow$flowernumb[i] > 0) 1 else 0
}
demoflow$did_flower<- did_flower


## get ES data
esdat4 <- demoflow[demoflow$species=="ELSC",]

# independent variables
tillernumb1 <-esdat4$tillernumb1 
year <- esdat4$year # random year effect

# dependent variable
did_flow <- esdat4$did_flower
N <- as.numeric(length(esdat4$did_flower))

# set up data and parameters to use
jags.data <- list("N","tillernumb1","did_flow", "year")
jags.param <- c( "a","bESdf", "rss","rss.new", "y", "prec1")

# the model
esdidflowmod <- function(){
  for(j in 1:4){y[j]~dnorm(0, prec1)} # random effect of year
  for (i in 1:N){
    did_flow[i]~dbern(mu[i])
    logit(mu[i]) <- a + bESdf*log(tillernumb1[i])+y[year[i]]
    
    res[i] <- (did_flow[i]-mu[i])^2 # posterior predictive check
    flow.new[i]~dbern(mu[i])
    res.new[i] <- (flow.new[i]-mu[i])^2
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bESdf~ dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
esdidflow<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                 n.iter=50000 ,model.file=esdidflowmod, n.thin=5, n.chains=3)

esdidflow # results

# posterior predictive check
esdidflow.paramlist <- esdidflow$BUGSoutput$sims.list

plot(esdidflow.paramlist$rss,esdidflow.paramlist$rss.new, main="C. did flower?",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(esdidflow.paramlist$rss>esdidflow.paramlist$rss.new)

# pull out 95% cred interval
bESdf <- sort(esdidflow.paramlist$bESdf) [375:14625]
aESdf <- sort(esdidflow.paramlist$a)[375:14625]

# Save
esdidflow_post <- as.data.frame(cbind(bESdf,aESdf))
write.csv(esdidflow_post, "esdidflow_post.csv")


### if flowering, how many flowers predicted by size

# get the proper data
esdat5 <- esdat[!esdat$flowernumb==0 ,]
esdat5 <- esdat5[!is.na(esdat5$flowernumb),]

# independent variables
year <- factor(esdat5$year) # year random effect
tillernumb1 <- esdat5$tillernumb1

# dependent variables
flowernumb <- log(esdat5$flowernumb)
N <- as.numeric(length(esdat5$flowernumb))

# set up data and parameters to use
jags.data <- list("N","flowernumb","tillernumb1", "year")
jags.param <- c( "a","bESfn","prec2","rss","rss.new", "prec1")

# the model
esflownumbmod <- function(){
  for(j in 1:4){y[j]~dnorm(0, prec1)} # year random effect 
  for (i in 1:N){
    flowernumb[i]~dnorm(mu[i],prec2)
    mu[i] <- a+ bESfn*log(tillernumb1[i])+y[year[i]]
    
    res[i] <- (flowernumb[i]-mu[i])^2 # posterior predictive check
    flnu.new[i]~dnorm(mu[i],prec2)
    res.new[i] <- (flnu.new[i]-mu[i])^2
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bESfn~dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  prec2~dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
esflownumb <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                  n.iter=50000 ,model.file=esflownumbmod, n.thin=5, n.chains=3)

esflownumb # results

# posterior check
esflownumb.paramlist <- esflownumb$BUGSoutput$sims.list

plot(esflownumb.paramlist$rss,esflownumb.paramlist$rss.new, main="D. flower number",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesina p-value
mean(esflownumb.paramlist$rss>esflownumb.paramlist$rss.new)

# 95% cred of vital rate params
bESfn <- sort(esflownumb.paramlist$bESfn) [375:14625]
aESfn <- sort(esflownumb.paramlist$a)[375:14625]

# save them
esflownumb_post <- as.data.frame(cbind(bESfn,aESfn))
write.csv(esflownumb_post, "esflownumb_post.csv")


# recruitment as a function of seed production
# this data is on github
recdat <- read.csv("data_lives_here/AlpDemo_recruit.csv")
esrec <- recdat[recdat$species=="ELSC",]
esdat6 <- esrec[!is.na(esrec$recruits_y2) ,]

# dependent variables
recruits_y2 <- esdat6$recruits_y2
N <- as.numeric(length(esdat6$recruits_y2))

# for binomial, number of trials (seeds produced)
seed <- round(esdat6$seed_numb)

# set up data and parameters to use
jags.data <- list("N","recruits_y2","seed")
jags.param <- c( "aESr", "rss", "rss.new")

# the model
esrecruit <- function(){
  for (i in 1:N){
    recruits_y2[i]~dbinom(mu[i],seed[i])
    logit(mu[i]) <- aESr 
    
    res[i] <- (recruits_y2[i]-mu[i])^2
    rec.new[i]~dbinom(mu[i],seed[i])
    res.new[i] <- (rec.new[i]-mu[i])^2
  }
  # prior
  aESr~ dnorm(0,1.0E-6)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
esrec<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
             n.iter=2000000 ,model.file=esrecruit, n.thin=5, n.chains=3)

esrec # results

# posterior check
esrec.paramlist <- esrec$BUGSoutput$sims.list

plot(esrec.paramlist$rss,esrec.paramlist$rss.new, main="E. recruitment",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

mean(esrec.paramlist$rss>esrec.paramlist$rss.new)

length(esrec.paramlist$aESr)

aESr <- sort(esrec.paramlist$aESr) [15000:585000]
hist(aESr)

esrec_post <- as.data.frame(cbind(aESr))
write.csv(esrec_post, "esrec_post.csv")


#########################################################################################
## Define functions
invlogit<-function(x){a <- exp(x)/(1+exp(x))
a[is.nan(a)]=1
return(a)
}

length2<-function(x){length(x[!is.na(x)]) }

## Assemble vital rate coefficient vectors

# load in posteriors from vital rate analyses
g_post <- read.csv("esgrowth_post.csv")
s_post <- read.csv("essurv_post.csv")
df_post <- read.csv("esdidflow_post.csv")
fn_post <- read.csv("esflownumb_post.csv")
r_post <- read.csv("esrec_post.csv")

# assemble posterior data frames
# naming goes "parameterSPECIESvitalrate"- e.g., "a" for intercept, "PA" for Poa alpina, and "s" for survival
# draw randomly 1000 times from 95% posterior
aESs <- sample(s_post$aESs, 1000)
bESs <- sample(s_post$bESs, 1000)
aESg <- sample(g_post$aESg, 1000)
bESg <- sample(g_post$bESg, 1000)
rESg <- sample(g_post$rESg, 1000)
aESdf <- sample(df_post$aESdf, 1000)
bESdf <- sample(df_post$bESdf, 1000)
aESfn <- sample(fn_post$aESfn, 1000)
bESfn <- sample(fn_post$bESfn, 1000)
aESr <- sample(r_post$aESr,1000)
min_siz <- rep(1, 1000)
max_size <- rep(90, 1000) # max six in experiment is 48!
seed_flower <- 11.04635449

ESdem_post <- as.data.frame(cbind(aESs,bESs,aESg, bESg,rESg,aESdf,bESdf,aESfn,bESfn, 
                                  aESr, min_siz,max_size, seed_flower))


# Above is the a base control representing core site with no cages 
# now use treatment effects to alter intercepts according to treatment effects in main analysis. 

# load in means of biomass, survival, and flowering. 
s_treat <- read.csv("data_lives_here/surv_postMeans95_ungherb.csv")
g_treat <- read.csv("data_lives_here/Bio_postMeans95_ungherb.csv")
fn_treat <- read.csv("data_lives_here/flow_postMeans95_ungherb.csv")

# filter to species and get effect sizes using high controls as reference
s_treat <- s_treat %>% filter(species=="Elymus scribneri") 
s_treat$eff <- s_treat$post/s_treat$post[10]-1 # s_treat$post[10] needs to be core control

g_treat <- g_treat %>% filter(species=="Elymus scribneri") 
g_treat$eff <- (g_treat$post/ g_treat$post[10])-1 # g_treat$post[10] needs to be core control

fn_treat <- fn_treat %>% filter(species=="Elymus scribneri") 
fn_treat$eff <- fn_treat$post/fn_treat$post[10]-1 # fn_treat$post[10] needs to be core control


# begin multiplying through the effect sizes
# only alter intercepts where we have treatment effects
# check that eff correspons to desired treatments

# core above
aESs_ha <- aESs+(abs(mean(aESs))*s_treat$eff[1])
aESg_ha <- aESg+abs(mean(aESg))*g_treat$eff[1]
aESfn_ha <- aESfn+abs(mean(aESfn))*fn_treat$eff[1]

ESdem_samp_ha <- as.data.frame(cbind(aESs_ha,bESs,aESg_ha, bESg,rESg,aESdf,bESdf,aESfn_ha,bESfn, 
                                     aESr, min_siz, max_size,seed_flower))


# core below
aESs_hb <- aESs+(abs(mean(aESs))*s_treat$eff[4])
aESg_hb <- aESg+abs(mean(aESg))*g_treat$eff[4]
aESfn_hb <- aESfn+abs(mean(aESfn))*fn_treat$eff[4]

ESdem_samp_hb <- as.data.frame(cbind(aESs_hb,bESs,aESg_hb, bESg,rESg,aESdf,bESdf,aESfn_hb,bESfn, 
                                     aESr, min_siz, max_size,seed_flower))

# core both exclosure
aESs_hab <- aESs+(abs(mean(aESs))*s_treat$eff[7])
aESg_hab <- aESg+abs(mean(aESg))*g_treat$eff[7]
aESfn_hab <- aESfn+abs(mean(aESfn))*fn_treat$eff[7]

ESdem_samp_hab <- as.data.frame(cbind(aESs_hab,bESs,aESg_hab, bESg,rESg,aESdf,bESdf,
                                      aESfn_hab,bESfn,aESr, min_siz, max_size,seed_flower))

# limit control
aESs_mc <- aESs+(abs(mean(aESs))*s_treat$eff[11])
aESg_mc <- aESg+abs(mean(aESg))*g_treat$eff[11]
aESfn_mc <- aESfn+abs(mean(aESfn))*fn_treat$eff[11]

ESdem_samp_mc <- as.data.frame(cbind(aESs_mc,bESs,aESg_mc, bESg,rESg,aESdf,bESdf,
                                     aESfn_mc,bESfn,aESr, min_siz, max_size,seed_flower))

# limit above
aESs_ma <- aESs+(abs(mean(aESs))*s_treat$eff[2])
aESg_ma <- aESg+abs(mean(aESg))*g_treat$eff[2]
aESfn_ma <- aESfn+abs(mean(aESfn))*fn_treat$eff[2]

ESdem_samp_ma <- as.data.frame(cbind(aESs_ma,bESs,aESg_ma, bESg,rESg,aESdf,bESdf,
                                     aESfn_ma,bESfn,aESr, min_siz, max_size,seed_flower))

# limit below
aESs_mb <- aESs+(abs(mean(aESs))*s_treat$eff[5])
aESg_mb <- aESg+abs(mean(aESg))*g_treat$eff[5]
aESfn_mb <- aESfn+abs(mean(aESfn))*fn_treat$eff[5]

ESdem_samp_mb <- as.data.frame(cbind(aESs_mb,bESs,aESg_mb, bESg,rESg,aESdf,bESdf,
                                     aESfn_mb,bESfn,aESr, min_siz, max_size,seed_flower))

# limit both
aESs_mab <- aESs+(abs(mean(aESs))*s_treat$eff[8])
aESg_mab <- aESg+abs(mean(aESg))*g_treat$eff[8]
aESfn_mab <- aESfn+abs(mean(aESfn))*fn_treat$eff[8]

ESdem_samp_mab <- as.data.frame(cbind(aESs_mab,bESs,aESg_mab, bESg,rESg,aESdf,bESdf,
                                      aESfn_mab,bESfn,aESr, min_siz, max_size,seed_flower))

# novel control
aESs_lc <- aESs+(abs(mean(aESs))*s_treat$eff[12])
aESg_lc <- aESg+abs(mean(aESg))*g_treat$eff[12]
aESfn_lc <- aESfn+abs(mean(aESfn))*fn_treat$eff[12]

ESdem_samp_lc <- as.data.frame(cbind(aESs_lc,bESs,aESg_lc, bESg,rESg,aESdf,bESdf,
                                     aESfn_lc,bESfn,aESr, min_siz, max_size,seed_flower))

# novel above
aESs_la <- aESs+(abs(mean(aESs))*s_treat$eff[3])
aESg_la <- aESg+abs(mean(aESg))*g_treat$eff[3]
aESfn_la <- aESfn+abs(mean(aESfn))*fn_treat$eff[3]

ESdem_samp_la <- as.data.frame(cbind(aESs_la,bESs,aESg_la, bESg,rESg,aESdf,bESdf,
                                     aESfn_la,bESfn,aESr, min_siz, max_size,seed_flower))

# novel below
aESs_lb <- aESs+(abs(mean(aESs))*s_treat$eff[6])
aESg_lb <- aESg+abs(mean(aESg))*g_treat$eff[6]
aESfn_lb <- aESfn+abs(mean(aESfn))*fn_treat$eff[6]

ESdem_samp_lb <- as.data.frame(cbind(aESs_lb,bESs,aESg_lb, bESg,rESg,aESdf,bESdf,
                                     aESfn_lb,bESfn,aESr, min_siz, max_size,seed_flower))

# novel both
aESs_lab <- aESs+(abs(mean(aESs))*s_treat$eff[9])
aESg_lab <- aESg+abs(mean(aESg))*g_treat$eff[9]
aESfn_lab <- aESfn+abs(mean(aESfn))*fn_treat$eff[9]

ESdem_samp_lab <- as.data.frame(cbind(aESs_lab,bESs,aESg_lab, bESg,rESg,aESdf,bESdf,
                                      aESfn_lab,bESfn,aESr, min_siz, max_size,seed_flower))


#########################################################################################
## Define vital rate functions
## These models correspond to statistical models fit by above (coefficients read in above)

#SURVIVAL AT SIZE X.
sx<-function(x,params){
  surv.mean<-params[1] + params[2]*log(x) 
  return(invlogit(surv.mean))
}

#PROBABILITY OF GROWTH FROM SIZE X TO Y
#This function truncates the density asscociation with x==0 and x>x.max
gxy<-function(x,y,params){
  grow.mean<-params[3] + params[4]*log(x) 
  grow<-dnbinom(x=y,mu=exp(grow.mean),size=params[5],log=F)
  truncLower<-dnbinom(x=0,mu=exp(grow.mean),size=params[5],log=F)
  truncUpper<-sum(dnbinom(x=(params[12]+1):7000,mu=exp(grow.mean),size=params[5],log=F))
  return(grow/(1-(truncLower+truncUpper)))
}

#SURVIVAL*GROWTH
pxy<-function(x,y,params){
  sx(x,params) * gxy(x,y,params)
}

# PROBABILITY OF FLOWERING
Pfx<-function(x,params){
  flow.mean<-params[6] + params[7]*log(x)  
  return(invlogit(flow.mean))
}

#NUMBER OF flowers
Nfx<-function(x,params){
  stalks.mean<-params[8] + params[9]*log(x) 
  return(exp(stalks.mean))
}


#SEED GERMINATION
Germ<-function(params){
  germ.mean<-params[10] 
  return(invlogit(germ.mean))
}

Fertx<-function(x,params){
  seedlings<-Pfx(x,params)*Nfx(x,params)*params[13]*Germ(params)
  return(seedlings)
}


################################################################################
## Put it all together
bigmatrix<-function(params){   
  
  matdim<-params[12]         ## bigmatrix dimension is max size 
  y<-1:params[12]
  
  # Fertility matrix
  Fmat<-matrix(0,matdim,matdim)
  Fmat[1,]<-Fertx(x=y,params=params) 
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,matdim,matdim)
  Tmat[]<-t(outer(y,y,pxy,params=params)) 
  
  # Put it all together
  MPMmat<-Tmat+Fmat #sum the Tmat & Fmat to get the whole matrix
  #this is the discrete approximation to the continuous kernel
  
  return(list(MPMmat=MPMmat, Fmat=Fmat,Tmat=Tmat))
}


#######################################################################
# calculate the lambdas

# core site controls
esdem_hc<- as.matrix(ESdem_post) # full posterior draws
esdem_hc_m <- apply(esdem_hc, 2, median) # median of each parameter after draws

lambda_hc <- 0
for(i in 1:1000){
  lambda_hc[i] <- Re(eigen(bigmatrix(esdem_hc[i,])$MPMmat)$value[1])
}

lambda_hc_m <-Re(eigen(bigmatrix(esdem_hc_m)$MPMmat)$value[1])

# core above
esdem_ha <- as.matrix(ESdem_samp_ha)
esdem_ha_m <- apply(esdem_ha, 2, median)

lambda_ha <- 0
for(i in 1:1000){
  lambda_ha[i] <- Re(eigen(bigmatrix(esdem_ha[i,])$MPMmat)$value[1])
}

lambda_ha_m <-Re(eigen(bigmatrix(esdem_ha_m)$MPMmat)$value[1])

# core below
esdem_hb <- as.matrix(ESdem_samp_hb)
esdem_hb_m <- apply(esdem_hb, 2, median)

lambda_hb <- 0
for(i in 1:1000){
  lambda_hb[i] <- Re(eigen(bigmatrix(esdem_hb[i,])$MPMmat)$value[1])
}

lambda_hb_m <-Re(eigen(bigmatrix(esdem_hb_m)$MPMmat)$value[1])

# core both
esdem_hab <- as.matrix(ESdem_samp_hab)
esdem_hab_m <- apply(esdem_hab, 2, median)

lambda_hab <- 0
for(i in 1:1000){
  lambda_hab[i] <- Re(eigen(bigmatrix(esdem_hab[i,])$MPMmat)$value[1])
}

lambda_hab_m <-Re(eigen(bigmatrix(esdem_hab_m)$MPMmat)$value[1])

# limit control
esdem_mc <- as.matrix(ESdem_samp_mc)
esdem_mc_m <- apply(esdem_mc, 2, median)

lambda_mc <- 0
for(i in 1:1000){
  lambda_mc[i] <- Re(eigen(bigmatrix(esdem_mc[i,])$MPMmat)$value[1])
}

lambda_mc_m <-Re(eigen(bigmatrix(esdem_mc_m)$MPMmat)$value[1])

# limit above
esdem_ma <- as.matrix(ESdem_samp_ma)
esdem_ma_m <- apply(esdem_ma, 2, median)

lambda_ma <- 0
for(i in 1:1000){
  lambda_ma[i] <- Re(eigen(bigmatrix(esdem_ma[i,])$MPMmat)$value[1])
}

lambda_ma_m <-Re(eigen(bigmatrix(esdem_ma_m)$MPMmat)$value[1])

# limit below
esdem_mb <- as.matrix(ESdem_samp_mb)
esdem_mb_m <- apply(esdem_mb, 2, median)

lambda_mb <- 0
for(i in 1:1000){
  lambda_mb[i] <- Re(eigen(bigmatrix(esdem_mb[i,])$MPMmat)$value[1])
}

lambda_mb_m <-Re(eigen(bigmatrix(esdem_mb_m)$MPMmat)$value[1])

# limit both
esdem_mab <- as.matrix(ESdem_samp_mab)
esdem_mab_m <- apply(esdem_mab, 2, median)

lambda_mab <- 0
for(i in 1:1000){
  lambda_mab[i] <- Re(eigen(bigmatrix(esdem_mab[i,])$MPMmat)$value[1])
}

lambda_mab_m <-Re(eigen(bigmatrix(esdem_mab_m)$MPMmat)$value[1])

# novel control
esdem_lc <- as.matrix(ESdem_samp_lc)
esdem_lc_m <- apply(esdem_lc, 2, median)

lambda_lc <- 0
for(i in 1:1000){
  lambda_lc[i] <- Re(eigen(bigmatrix(esdem_lc[i,])$MPMmat)$value[1])
}

lambda_lc_m <-Re(eigen(bigmatrix(esdem_lc_m)$MPMmat)$value[1])

# novel above
esdem_la <- as.matrix(ESdem_samp_la)
esdem_la_m <- apply(esdem_la, 2, median)

lambda_la <- 0
for(i in 1:1000){
  lambda_la[i] <- Re(eigen(bigmatrix(esdem_la[i,])$MPMmat)$value[1])
}

lambda_la_m <-Re(eigen(bigmatrix(esdem_la_m)$MPMmat)$value[1])

# novel below
esdem_lb <- as.matrix(ESdem_samp_lb)
esdem_lb_m <- apply(esdem_lb, 2, median)

lambda_lb <- 0
for(i in 1:1000){
  lambda_lb[i] <- Re(eigen(bigmatrix(esdem_lb[i,])$MPMmat)$value[1])
}

lambda_lb_m <-Re(eigen(bigmatrix(esdem_lb_m)$MPMmat)$value[1])

# novel both
esdem_lab <- as.matrix(ESdem_samp_lab)
esdem_lab_m <- apply(esdem_lab, 2, median)

lambda_lab <- 0
for(i in 1:1000){
  lambda_lab[i] <- Re(eigen(bigmatrix(esdem_lab[i,])$MPMmat)$value[1])
}

lambda_lab_m <-Re(eigen(bigmatrix(esdem_lab_m)$MPMmat)$value[1])


# assemble posterior and median lambda datasets
eslambda<- as.data.frame(cbind(lambda_hc,lambda_ha,lambda_hb,lambda_hab,
                               lambda_lc,lambda_la,lambda_lb,lambda_lab,
                               lambda_mc,lambda_ma,lambda_mb,lambda_mab))
eslambda_m <- as.data.frame(cbind(lambda_hc_m,lambda_ha_m,lambda_hb_m,lambda_hab_m,
                                  lambda_lc_m,lambda_la_m,lambda_lb_m,lambda_lab_m,
                                  lambda_mc_m,lambda_ma_m,lambda_mb_m,lambda_mab_m))

write.csv(eslambda, "eslambda_nat_withExpTreat.csv")
write.csv(eslambda_m, "eslambda_nat_withExpTreat_median.csv")


#################################################################################
## LTRE calculations

# create midpoints between reference core/control and other treatment combinations
ha_midpoint <- (esdem_ha_m+esdem_hc_m)/2
hb_midpoint <- (esdem_hb_m+esdem_hc_m)/2
hab_midpoint <- (esdem_hab_m+esdem_hc_m)/2
mc_midpoint <- (esdem_mc_m+esdem_hc_m)/2
ma_midpoint <- (esdem_ma_m+esdem_hc_m)/2
mb_midpoint <- (esdem_mb_m+esdem_hc_m)/2
mab_midpoint <- (esdem_mab_m+esdem_hc_m)/2
lc_midpoint <- (esdem_lc_m+esdem_hc_m)/2
la_midpoint <- (esdem_la_m+esdem_hc_m)/2
lb_midpoint <- (esdem_lb_m+esdem_hc_m)/2
lab_midpoint <- (esdem_lab_m+esdem_hc_m)/2

# Midpoint lambda
midpoint.lam_ha<-Re(eigen(bigmatrix(ha_midpoint)$MPMmat)$values[1])
midpoint.lam_hb<-Re(eigen(bigmatrix(hb_midpoint)$MPMmat)$values[1])
midpoint.lam_hab<-Re(eigen(bigmatrix(hab_midpoint)$MPMmat)$values[1])
midpoint.lam_mc<-Re(eigen(bigmatrix(mc_midpoint)$MPMmat)$values[1])
midpoint.lam_ma<-Re(eigen(bigmatrix(ma_midpoint)$MPMmat)$values[1])
midpoint.lam_mb<-Re(eigen(bigmatrix(mb_midpoint)$MPMmat)$values[1])
midpoint.lam_mab<-Re(eigen(bigmatrix(mab_midpoint)$MPMmat)$values[1])
midpoint.lam_lc<-Re(eigen(bigmatrix(lc_midpoint)$MPMmat)$values[1])
midpoint.lam_la<-Re(eigen(bigmatrix(la_midpoint)$MPMmat)$values[1])
midpoint.lam_lb<-Re(eigen(bigmatrix(lb_midpoint)$MPMmat)$values[1])
midpoint.lam_lab<-Re(eigen(bigmatrix(lab_midpoint)$MPMmat)$values[1])


# LTRE sensitivities prellocation
LTRE.sens_ha <-c();LTRE.sens_hb <-c();LTRE.sens_hab <-c()
LTRE.sens_mc <-c();LTRE.sens_ma <-c();LTRE.sens_mb <-c();LTRE.sens_mab <-c()
LTRE.sens_lc <-c();LTRE.sens_la <-c();LTRE.sens_lb <-c();LTRE.sens_lab <-c()

# perturbation
perturbation<-0.001

# run the LTREs
for(j in 1:10){
  LTRE.mid<-ha_midpoint
  LTRE.mid[j]<-ha_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_ha[j]<-(LTRE.lam-midpoint.lam_ha)/perturbation
}
for(j in 1:10){
  LTRE.mid<-hb_midpoint
  LTRE.mid[j]<-hb_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_hb[j]<-(LTRE.lam-midpoint.lam_hb)/perturbation
}
for(j in 1:10){
  LTRE.mid<-hab_midpoint
  LTRE.mid[j]<-hab_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_hab[j]<-(LTRE.lam-midpoint.lam_hab)/perturbation
}
for(j in 1:10){
  LTRE.mid<-lc_midpoint
  LTRE.mid[j]<-lc_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_lc[j]<-(LTRE.lam-midpoint.lam_lc)/perturbation
}
for(j in 1:10){
  LTRE.mid<-la_midpoint
  LTRE.mid[j]<-la_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_la[j]<-(LTRE.lam-midpoint.lam_la)/perturbation
}
for(j in 1:10){
  LTRE.mid<-lb_midpoint
  LTRE.mid[j]<-lb_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_lb[j]<-(LTRE.lam-midpoint.lam_lb)/perturbation
}
for(j in 1:10){
  LTRE.mid<-lab_midpoint
  LTRE.mid[j]<-lab_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_lab[j]<-(LTRE.lam-midpoint.lam_lab)/perturbation
}
for(j in 1:10){
  LTRE.mid<-mc_midpoint
  LTRE.mid[j]<-mc_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_mc[j]<-(LTRE.lam-midpoint.lam_mc)/perturbation
}
for(j in 1:10){
  LTRE.mid<-ma_midpoint
  LTRE.mid[j]<-ma_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_ma[j]<-(LTRE.lam-midpoint.lam_ma)/perturbation
}
for(j in 1:10){
  LTRE.mid<-mb_midpoint
  LTRE.mid[j]<-mb_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_mb[j]<-(LTRE.lam-midpoint.lam_mb)/perturbation
}
for(j in 1:10){
  LTRE.mid<-mab_midpoint
  LTRE.mid[j]<-mab_midpoint[j]+perturbation
  LTRE.lam<-Re(eigen(bigmatrix(LTRE.mid)$MPMmat)$values[1])
  LTRE.sens_mab[j]<-(LTRE.lam-midpoint.lam_mab)/perturbation
}

# lambda differences
delta.lam.ha <- lambda_ha_m-lambda_hc_m
delta.lam.hb <- lambda_hb_m-lambda_hc_m
delta.lam.hab <- lambda_hab_m-lambda_hc_m
delta.lam.lc <- lambda_lc_m-lambda_hc_m
delta.lam.la <- lambda_la_m-lambda_hc_m
delta.lam.lb <- lambda_lb_m-lambda_hc_m
delta.lam.lab <- lambda_lab_m-lambda_hc_m
delta.lam.mc <- lambda_mc_m-lambda_hc_m
delta.lam.ma <- lambda_ma_m-lambda_hc_m
delta.lam.mb <- lambda_mb_m-lambda_hc_m
delta.lam.mab <- lambda_mab_m-lambda_hc_m

# sum the sensitivities
sum.ha <- sum(na.omit(c(LTRE.sens_ha[1:10]*(esdem_ha_m[1:10]-esdem_hc_m[1:10]))))
sum.hb <- sum(na.omit(c(LTRE.sens_hb[1:10]*(esdem_hb_m[1:10]-esdem_hc_m[1:10]))))
sum.hab <- sum(na.omit(c(LTRE.sens_hab[1:10]*(esdem_hab_m[1:10]-esdem_hc_m[1:10]))))
sum.lc <- sum(na.omit(c(LTRE.sens_lc[1:10]*(esdem_lc_m[1:10]-esdem_hc_m[1:10]))))
sum.la <- sum(na.omit(c(LTRE.sens_la[1:10]*(esdem_la_m[1:10]-esdem_hc_m[1:10]))))
sum.lb <- sum(na.omit(c(LTRE.sens_lb[1:10]*(esdem_lb_m[1:10]-esdem_hc_m[1:10]))))
sum.lab <- sum(na.omit(c(LTRE.sens_lab[1:10]*(esdem_lab_m[1:10]-esdem_hc_m[1:10]))))
sum.mc <- sum(na.omit(c(LTRE.sens_mc[1:10]*(esdem_mc_m[1:10]-esdem_hc_m[1:10]))))
sum.ma <- sum(na.omit(c(LTRE.sens_ma[1:10]*(esdem_ma_m[1:10]-esdem_hc_m[1:10]))))
sum.mb <- sum(na.omit(c(LTRE.sens_mb[1:10]*(esdem_mb_m[1:10]-esdem_hc_m[1:10]))))
sum.mab <- sum(na.omit(c(LTRE.sens_mab[1:10]*(esdem_mab_m[1:10]-esdem_hc_m[1:10]))))

# calculate the LTREs
ltre.ha <- (na.omit(c(LTRE.sens_ha[1:10]*(esdem_ha_m[1:10]-esdem_hc_m[1:10]))))
ltre.hb <- (na.omit(c(LTRE.sens_hb[1:10]*(esdem_hb_m[1:10]-esdem_hc_m[1:10]))))
ltre.hab <- (na.omit(c(LTRE.sens_hab[1:10]*(esdem_hab_m[1:10]-esdem_hc_m[1:10]))))
ltre.lc <- (na.omit(c(LTRE.sens_lc[1:10]*(esdem_lc_m[1:10]-esdem_hc_m[1:10]))))
ltre.la <- (na.omit(c(LTRE.sens_la[1:10]*(esdem_la_m[1:10]-esdem_hc_m[1:10]))))
ltre.lb <- (na.omit(c(LTRE.sens_lb[1:10]*(esdem_lb_m[1:10]-esdem_hc_m[1:10]))))
ltre.lab <- (na.omit(c(LTRE.sens_lab[1:10]*(esdem_lab_m[1:10]-esdem_hc_m[1:10]))))
ltre.mc <- (na.omit(c(LTRE.sens_mc[1:10]*(esdem_mc_m[1:10]-esdem_hc_m[1:10]))))
ltre.ma <- (na.omit(c(LTRE.sens_ma[1:10]*(esdem_ma_m[1:10]-esdem_hc_m[1:10]))))
ltre.mb <- (na.omit(c(LTRE.sens_mb[1:10]*(esdem_mb_m[1:10]-esdem_hc_m[1:10]))))
ltre.mab <- (na.omit(c(LTRE.sens_mab[1:10]*(esdem_mab_m[1:10]-esdem_hc_m[1:10]))))

# put them together
LTREresults <- as.data.frame(rbind(ltre.ha, ltre.hb, ltre.hab,ltre.mc, 
                                   ltre.ma,ltre.mb, ltre.mab,ltre.lc,
                                   ltre.la, ltre.lb, ltre.lab))

# add names
names(LTREresults) <- c("Survival", "b_surv", "Growth", "b_grow", "r_grow",
                        "a_df", "b_df", "Flower number", "b_fn", "a_rec")

# row names
LTREresults$model <- c("LTRE.sens_ha", "LTRE.sens_hb"," LTRE.sens_hab","LTRE.sens_mc", 
                       "LTRE.sens_ma","LTRE.sens_mb", "LTRE.sens_mab","LTRE.sens_lc",
                       "LTRE.sens_la", "LTRE.sens_lb", "LTRE.sens_lab")

# the LTRE sums
LTREresults$LTREsums <- c(sum.ha, sum.hb, sum.hab, sum.mc, sum.ma, sum.mb,sum.mab,
                          sum.lc, sum.la, sum.lb, sum.lab)

# add delta lambda
LTREresults$delt.lam <- c(delta.lam.ha, delta.lam.hb,delta.lam.hab, delta.lam.mc, delta.lam.ma,
                          delta.lam.mb,delta.lam.mab, delta.lam.lc, delta.lam.la, delta.lam.lb,
                          delta.lam.lab)

# add in site factors
LTREresults$site <- c("core", "core","core", "limit","limit","limit","limit","novel",
                      "novel","novel","novel")

# repeat the species
LTREresults$species <- rep("Elymus scribneri", 11)

# add in exclosure treatments
LTREresults$treat <- c("above", "below", "both", "control", "above", "below", "both",
                       "control", "above", "below", "both")

## save  the LTRE results data
write.csv(LTREresults, "LTREres_es_ungherb.csv")

