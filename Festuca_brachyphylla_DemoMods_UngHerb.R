### Demographic analyses Festuca brachyphylla


# Packages
library(R2jags); library(car);library(tidyverse)

# load in demographic data
demodat <- read.csv("data_lives_here/MasterDemoDataset.csv")
head(demodat)

# make year and plot number a factor
demodat$year <- as.factor(demodat$year)
demodat$plotnumb <- as.factor(demodat$plotnumb)

# pull out Festuca brachyphylla (FEBR/FB)
fbdat <- demodat[demodat$species=="FEBR",]

### vital rate calculations

# Parameterize the growth model for FB

# remove NA's in FB data
fbdat2 <- fbdat[!is.na(fbdat$tillernumb1 & fbdat$tillernumb2),]

tillernumb1 <- fbdat2$tillernumb1 # tiller number at in year t
tillernumb2 <- fbdat2$tillernumb2# tiller number at at year t+1

year <- factor(fbdat2$year) # for random effects
N <- as.numeric(length(fbdat2$tillernumb1)) # get length of dataset

# set up data and parameters to use
jags.data <- list("N","tillernumb1","tillernumb2", "year")
jags.param <- c( "a","bFBg","r", "rss", "rss.new", "prec1")

# the model
fbgrowthmod <- function(){
  for(j in 1:3){aar[j]~dnorm(0, prec1)} # year random effect
  for (i in 1:N){
    tillernumb2[i]~dnegbin(p[i],r)
    log(mu[i]) <- a+ bFBg*log(tillernumb1[i])+aar[year[i]]
    p[i] <- r/(r+mu[i])
    
    res[i] <- (tillernumb2[i]-mu[i])^2 # posterior predictive check
    tiller.new[i]~dnegbin(p[i],r)
    res.new[i] <- (tiller.new[i]-mu[i])^2
  }
  # priors
  r~ dunif(0,50)
  a~ dnorm(0,1.0E-6)
  bFBg~dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  # derived parameters
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
fbgrowth<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                n.iter=50000 ,model.file=fbgrowthmod, n.thin=5, n.chains=3)

fbgrowth # results

# posterior predictive check
fbgrowth.paramlist <- fbgrowth$BUGSoutput$sims.list

plot(fbgrowth.paramlist$rss,fbgrowth.paramlist$rss.new, main="A. growth",
     xlab="SSQ observed", ylab="SSQ simulated") ## looks good
abline(0,1)

# Bayesian p-value
mean(fbgrowth.paramlist$rss>fbgrowth.paramlist$rss.new)

## take the 95% cred. intevals of each
bFBg <- sort(fbgrowth.paramlist$bFBg) [375:14625]
aFBg <- sort(fbgrowth.paramlist$a)[375:14625]
rFBg <- sort(fbgrowth.paramlist$r)[375:14625]

fbgrowth_post <- as.data.frame(cbind(bFBg,aFBg,rFBg))
write.csv(fbgrowth_post, "fbgrowth_post.csv")


### survival vital rates for FB 
fbdat3 <- fbdat[!is.na(fbdat$survto2 & fbdat$tillernumb1),] #subset out individuals still dead

# dependent variables
survto2 <- fbdat3$survto2
N <- as.numeric(length(fbdat3$survto2))

# independent variables
tillernumb1 <- fbdat3$tillernumb1
year <- factor(fbdat3$year) # year random effect

# set up data and parameters to use
jags.data <- list("N","tillernumb1","survto2", "year")
jags.param <- c( "a","bFBs","rss", "rss.new", "prec1")

# the model 
fbsurvmod <- function(){
  for(j in 1:3){aar[j]~dnorm(0, prec1)} # year random effect
  for (i in 1:N){
    survto2[i]~dbern(mu[i])
    logit(mu[i]) <- a + bFBs*log(tillernumb1[i])+aar[year[i]]
    
    res[i] <- pow((survto2[i]-mu[i]),2) # posterior predictive checks
    surv.new[i]~dbern(mu[i])
    res.new[i] <- pow((surv.new[i]-mu[i]),2)
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bFBs~ dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through jags
fbsurv<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
              n.iter=100000 ,model.file=fbsurvmod, n.thin=5, n.chains=3)

fbsurv # results

#posterior predictive check
fbsurv.paramlist <- fbsurv$BUGSoutput$sims.list

plot(fbsurv.paramlist$rss,fbsurv.paramlist$rss.new,main="B. survival",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(fbsurv.paramlist$rss>fbsurv.paramlist$rss.new)

## pull out 95% cred
bFBs <- sort(fbsurv.paramlist$bFBs) [750:29250]
aFBs <- sort(fbsurv.paramlist$a) [750:29250]

bsurv_post <- as.data.frame(cbind(bFBs,aFBs))
write.csv(fbsurv_post, "fbsurv_post.csv")


## probabliltiy of flowering vital rates

# create a 0,1 column for flowering
demoflow <- demodat[!is.na(demodat$flowernumb),]
did_flower <- 0
for(i in 1:length(demoflow$flowernumb)){
  did_flower[i]<- if(demoflow$flowernumb[i] > 0) 1 else 0
}
demoflow$did_flower<- did_flower

## FB sort
fbdat4 <- demoflow[demoflow$species=="FEBR",]
fbdat4 <- fbdat4[!is.na(fbdat4$tillernumb1),]

# independent variable- tiller number is year of flowering
tillernumb1 <- fbdat4$tillernumb1
year <- factor(fbdat4$year) # random effect

# dependent variables 
did_flow <- fbdat4$did_flower
N <- as.numeric(length(fbdat4$did_flower))

# set up data and parameters to use
jags.data <- list("N","tillernumb1","did_flow", "year")
jags.param <- c( "a","bFBdf", "rss","rss.new", "prec1")

# the model 
fbdidflowmod <- function(){
  for(j in 1:4){aar[j]~dnorm(0, prec1)}  # random effect of years
  for (i in 1:N){
    did_flow[i]~dbern(mu[i])
    logit(mu[i]) <- a + bFBdf*log(tillernumb1[i])+aar[year[i]]
    
    res[i] <- (did_flow[i]-mu[i])^2 # posterior predictive check
    flow.new[i]~dbern(mu[i])
    res.new[i] <- (flow.new[i]-mu[i])^2
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bFBdf~ dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
fbdidflow<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                 n.iter=50000 ,model.file=fbdidflowmod, n.thin=5, n.chains=3)

fbdidflow # results

# posterior predictive check
fbdidflow.paramlist <- fbdidflow$BUGSoutput$sims.list

plot(fbdidflow.paramlist$rss,fbdidflow.paramlist$rss.new,main="C. did flower?",
     xlab="SSQ observed", ylab="SSQ simulated")## good
abline(0,1)

# Bayesian p-value
mean(fbdidflow.paramlist$rss>fbdidflow.paramlist$rss.new)

# 95% cred of vital rate params
bFBdf <- sort(fbdidflow.paramlist$bFBdf) [375:14625]
aFBdf <- sort(fbdidflow.paramlist$a)[375:14625]

fbdidflow_post <- as.data.frame(cbind(bFBdf,aFBdf))
write.csv(fbdidflow_post, "fbdidflow_post.csv")


### if flowering, how many flowers predicted by size
# get the proper data
fbdat5 <- fbdat[!fbdat$flowernumb==0 ,]
fbdat5 <- fbdat5[!is.na(fbdat5$flowernumb & fbdat5$tillernumb1),]

# independent variable
tillernumb1 <- fbdat5$tillernumb1
year <- factor(fbdat5$year) # year random effect

# dependent variable
flowernumb <- log(fbdat5$flowernumb)
N <- as.numeric(length(fbdat5$flowernumb))

# set up data and parameters to use
jags.data <- list("N","flowernumb","tillernumb1", "year")
jags.param <- c( "a","bFBfn", "rss","rss.new","prec2", "prec1")

# the model
fbflownumbmod <- function(){
  for(j in 1:4){aar[j]~dnorm(0,prec1)}# year random effect
  for (i in 1:N){
    flowernumb[i]~dnorm(mu[i],prec2)
    mu[i] <- a+ bFBfn*log(tillernumb1[i])+aar[year[i]]
    
    res[i] <- (flowernumb[i]-mu[i])^2# posterior predictive check
    flnu.new[i]~dnorm(mu[i],prec2)
    res.new[i] <- (flnu.new[i]-mu[i])^2
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bFBfn~dnorm(0,1.0E-6)
  prec1~dgamma(0.001, 0.001)
  prec2~dgamma(0.001, 0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
fbflownumb<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                  n.iter=50000 ,model.file=fbflownumbmod, n.thin=5, n.chains=3)

fbflownumb # results

# posterior check
fbflownumb.paramlist <- fbflownumb$BUGSoutput$sims.list

plot(fbflownumb.paramlist$rss,fbflownumb.paramlist$rss.new, main="D. flower number",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

#Bayesian p-value
mean(fbflownumb.paramlist$rss>fbflownumb.paramlist$rss.new)

# 95% cred of vital rate params
bFBfn <- sort(fbflownumb.paramlist$bFBfn) [375:14625]
aFBfn <- sort(fbflownumb.paramlist$a)[375:14625]

fbflownumb_post <- as.data.frame(cbind(bFBfn,aFBfn))
write.csv(fbflownumb_post, "fbflownumb_post.csv")


# recruitment as a function of seed production
# this data is on github
recdat <- read.csv("data_lives_here/AlpDemo_recruit.csv")
fbrecr <- recdat[recdat$species=="FEBR",]

# FB recruitment
fbdat6 <- fbrecr[! is.na(fbrecr$recruits_y2) ,]

# dependent variables
recruits_y2 <- fbdat6$recruits_y2
N <- as.numeric(length(fbdat6$recruits_y2))

# for binomial, number of trials (seeds produced)
seed <- round(fbdat6$seed_numb)

# set up data and parameters to use
jags.data <- list("N","recruits_y2","seed")
jags.param <- c( "aFBr", "rss","rss.new")

fbrecruit <- function(){
  for (i in 1:N){
    recruits_y2[i]~dbinom(mu[i],seed[i])
    logit(mu[i]) <- aFBr 
    
    res[i] <- (recruits_y2[i]-mu[i])^2# for posterior predictive check
    rec.new[i]~dbinom(mu[i],seed[i])
    res.new[i] <- (rec.new[i]-mu[i])^2
  }
  # priors
  aFBr~ dnorm(0,1.0E-6)
  rss <- sum(res[])
  rss.new <- sum(res.new[]) 
}

# run through jags
fbrec<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
             n.iter=4000000 ,model.file=fbrecruit, n.thin=5, n.chains=3)

fbrec # results

# posterior check
fbrec.paramlist <- fbrec$BUGSoutput$sims.list

plot(fbrec.paramlist$rss,fbrec.paramlist$rss.new, main="D. recruitment",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(fbrec.paramlist$rss>fbrec.paramlist$rss.new)

# 95% cred of recruit parm
aFBr <- sort(fbrec.paramlist$aFBr)[30000:1170000]

brec_post <- as.data.frame(cbind(aFBr))
write.csv(fbrec_post, "fbrec_post.csv")

# end vital rate set up
#########################################################################################
## Define functions
invlogit<-function(x){a <- exp(x)/(1+exp(x))
a[is.nan(a)]=1
return(a)
}

length2<-function(x){length(x[!is.na(x)]) }

## Assemble vital rate coefficient vectors

# load in posteriors from vital rate analyses
g_post <- read.csv("fbgrowth_post.csv")
s_post <- read.csv("fbsurv_post.csv")
df_post <- read.csv("fbdidflow_post.csv")
fn_post <- read.csv("fbflownumb_post.csv")
r_post <- read.csv("fbrec_post.csv")

# assemble posterior data frames
# naming goes "parameterSPECIESvitalrate"- e.g., "a" for intercept, "PA" for Poa alpina, and "s" for survival
# draw randomly 1000 times from 95% posterior

aFBs <- sample(s_post$aFBs, 1000)
bFBs <- sample(s_post$bFBs, 1000)
aFBg <- sample(g_post$aFBg, 1000)
bFBg <- sample(g_post$bFBg, 1000)
rFBg <- sample(g_post$rFBg, 1000)
aFBdf <- sample(df_post$aFBdf, 1000)
bFBdf <- sample(df_post$bFBdf, 1000)
aFBfn <- sample(fn_post$aFBfn, 1000)
bFBfn <- sample(fn_post$bFBfn, 1000)
aFBr <- sample(r_post$aFBr,1000)
min_siz <- rep(1, 1000)
max_size <- rep(133, 1000)# max of experimental data 133
seed_flower <- 9.835847763

FBdem_samp <- as.data.frame(cbind(aFBs,bFBs,aFBg, bFBg,rFBg,aFBdf,bFBdf,aFBfn,bFBfn, 
                                  aFBr, min_siz, max_size,seed_flower))

# Above is the a base control representing core site with no cages 
# now use treatment effects to alter intercepts according to treatment effects in main analysis. 

# load in means of biomass, survival, and flowering. 
s_treat <- read.csv("data_lives_here/surv_postMeans95_ungherb.csv")
g_treat <- read.csv("data_lives_here/Bio_postMeans95_ungherb.csv")
fn_treat <- read.csv("data_lives_here/flow_postMeans95_ungherb.csv")

# filter to species and get effect sizes using high controls as 
s_treat <- s_treat %>% filter(species=="Festuca brachyphylla") 
s_treat$eff <- s_treat$post/s_treat$post[10]-1 # s_treat$post[10] needs to be core control

g_treat <- g_treat %>% filter(species=="Festuca brachyphylla") 
g_treat$eff <- (g_treat$post/ g_treat$post[10])-1 # g_treat$post[10] needs to be core control

fn_treat <- fn_treat %>% filter(species=="Festuca brachyphylla") 
fn_treat$eff <- fn_treat$post/fn_treat$post[10]-1 # fn_treat$post[10] needs to be core control

# begin multiplying through the effect sizes
# only alter intercepts where we have treatment effects
# check that eff corresponds to desired treatments

# core above
aFBs_ha <- aFBs+(abs(mean(aFBs))*s_treat$eff[1])
aFBg_ha <- aFBg+abs(mean(aFBg))*g_treat$eff[1]
aFBfn_ha <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[1]

FBdem_samp_ha <- as.data.frame(cbind(aFBs_ha,bFBs,aFBg_ha, bFBg,rFBg,aFBdf,bFBdf,
                                     aFBfn_ha,bFBfn,aFBr, min_siz, max_size,seed_flower))

# core below
aFBs_hb <- aFBs+abs(mean(aFBs))*s_treat$eff[4]
aFBg_hb <- aFBg+abs(mean(aFBg))*g_treat$eff[4]
aFBfn_hb <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[4]

FBdem_samp_hb <- as.data.frame(cbind(aFBs_hb,bFBs,aFBg_hb, bFBg,rFBg,aFBdf,bFBdf,aFBfn_hb,bFBfn, 
                                     aFBr, min_siz, max_size,seed_flower))

# core both
aFBs_hab <- aFBs+abs(mean(aFBs))*s_treat$eff[7]
aFBg_hab <- aFBg+abs(mean(aFBg))*g_treat$eff[7]
aFBfn_hab <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[7]

FBdem_samp_hab <- as.data.frame(cbind(aFBs_hab, bFBs,aFBg_hab, bFBg,rFBg,aFBdf,bFBdf,
                                      aFBfn_hab,bFBfn, aFBr, min_siz, max_size,seed_flower))

# limit  control
aFBs_mc <- aFBs+abs(mean(aFBs))*s_treat$eff[11]
aFBg_mc <- aFBg+abs(mean(aFBg))*g_treat$eff[11]
aFBfn_mc <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[11]

FBdem_samp_mc <- as.data.frame(cbind(aFBs_mc,bFBs,aFBg_mc, bFBg,rFBg,aFBdf,bFBdf,aFBfn_mc,bFBfn, 
                                     aFBr, min_siz, max_size,seed_flower))

# limit above
aFBs_ma <- aFBs+abs(mean(aFBs))*s_treat$eff[2]
aFBg_ma <- aFBg+abs(mean(aFBg))*g_treat$eff[2]
aFBfn_ma <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[2]

FBdem_samp_ma <- as.data.frame(cbind(aFBs_ma,bFBs,aFBg_ma, bFBg,rFBg,aFBdf,bFBdf,aFBfn_ma,bFBfn, 
                                     aFBr, min_siz, max_size,seed_flower))

# limit below
aFBs_mb <- aFBs+abs(mean(aFBs))*s_treat$eff[5]
aFBg_mb <- aFBg+abs(mean(aFBg))*g_treat$eff[5]
aFBfn_mb <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[5]

FBdem_samp_mb <- as.data.frame(cbind(aFBs_mb,bFBs,aFBg_mb, bFBg,rFBg,aFBdf,bFBdf,aFBfn_mb,bFBfn, 
                                     aFBr, min_siz, max_size,seed_flower))

#limit both
aFBs_mab <- aFBs+abs(mean(aFBs))*s_treat$eff[8]
aFBg_mab <- aFBg+abs(mean(aFBg))*g_treat$eff[8]
aFBfn_mab <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[8]

FBdem_samp_mab <- as.data.frame(cbind(aFBs_mab,bFBs  ,aFBg_mab, bFBg,rFBg,aFBdf,bFBdf,
                                      aFBfn_mab,bFBfn,aFBr, min_siz, max_size,seed_flower))

# novel control
aFBs_lc <- aFBs+abs(mean(aFBs))*s_treat$eff[12]
aFBg_lc <- aFBg+abs(mean(aFBg))*g_treat$eff[12]
aFBfn_lc <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[12]

FBdem_samp_lc <- as.data.frame(cbind(aFBs_lc,bFBs,aFBg_lc, bFBg,rFBg,aFBdf,bFBdf,aFBfn_lc,bFBfn, 
                                     aFBr, min_siz, max_size,seed_flower))

# novel above
aFBs_la <- aFBs+abs(mean(aFBs))*s_treat$eff[3]
aFBg_la <- aFBg+abs(mean(aFBg))*g_treat$eff[3]
aFBfn_la <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[3]

FBdem_samp_la <- as.data.frame(cbind(aFBs_la,bFBs,aFBg_la, bFBg,rFBg,aFBdf,bFBdf,aFBfn_la,bFBfn, 
                                     aFBr, min_siz, max_size,seed_flower))

# novel below
aFBs_lb <- aFBs+abs(mean(aFBs))*s_treat$eff[6]
aFBg_lb <- aFBg+abs(mean(aFBg))*g_treat$eff[6]
aFBfn_lb <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[6]

FBdem_samp_lb <- as.data.frame(cbind(aFBs_lb,bFBs,aFBg_lb, bFBg,rFBg,aFBdf,bFBdf,aFBfn_lb,bFBfn, 
                                     aFBr, min_siz, max_size,seed_flower))

# novel both
aFBs_lab <- aFBs+abs(mean(aFBs))*s_treat$eff[9]
aFBg_lab <- aFBg+abs(mean(aFBg))*g_treat$eff[9]
aFBfn_lab <- aFBfn+abs(mean(aFBfn))*fn_treat$eff[9]

FBdem_samp_lab <- as.data.frame(cbind(aFBs_lab,bFBs,aFBg_lab, bFBg,rFBg,aFBdf,bFBdf,
                                      aFBfn_lab,bFBfn,aFBr, min_siz, max_size,seed_flower))

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

# FERTILITY function
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
fbdem_hc<- as.matrix(FBdem_samp) # full posterior draws
fbdem_hc_m <- apply(fbdem_hc, 2, median) # median of each parameter after draws

lambda_hc <- 0
for(i in 1:1000){
  lambda_hc[i] <-Re(eigen(bigmatrix(fbdem_hc[i,])$MPMmat)$value[1])
}

lambda_hc_m <-Re(eigen(bigmatrix(fbdem_hc_m)$MPMmat)$value[1])

# core above
fbdem_ha <- as.matrix(FBdem_samp_ha)
fbdem_ha_m <- apply(fbdem_ha, 2, median)

lambda_ha <- 0
for(i in 1:1000){
  lambda_ha[i] <- Re(eigen(bigmatrix(fbdem_ha[i,])$MPMmat)$value[1])
}

lambda_ha_m <-Re(eigen(bigmatrix(fbdem_ha_m)$MPMmat)$value[1])

# core below
fbdem_hb <- as.matrix(FBdem_samp_hb)
fbdem_hb_m <- apply(fbdem_hb, 2, median)

lambda_hb <- 0
for(i in 1:1000){
  lambda_hb[i] <-Re(eigen(bigmatrix(fbdem_hb[i,])$MPMmat)$value[1])
}

lambda_hb_m <-Re(eigen(bigmatrix(fbdem_hb_m)$MPMmat)$value[1])

# core both
fbdem_hab<- as.matrix(FBdem_samp_hab)
fbdem_hab_m <- apply(fbdem_hab, 2, median)

lambda_hab <- 0
for(i in 1:1000){
  lambda_hab[i] <-Re(eigen(bigmatrix(fbdem_hab[i,])$MPMmat)$value[1])
}

lambda_hab_m <-Re(eigen(bigmatrix(fbdem_hab_m)$MPMmat)$value[1])

# limit control
fbdem_mc<- as.matrix(FBdem_samp_mc)
fbdem_mc_m <- apply(fbdem_mc, 2, median)

lambda_mc <- 0
for(i in 1:1000){
  lambda_mc[i] <-Re(eigen(bigmatrix(fbdem_mc[i,])$MPMmat)$value[1])
}

lambda_mc_m <-Re(eigen(bigmatrix(fbdem_mc_m)$MPMmat)$value[1])

# limit above
fbdem_ma<- as.matrix(FBdem_samp_ma)
fbdem_ma_m <- apply(fbdem_ma, 2, median)

lambda_ma <- 0
for(i in 1:1000){
  lambda_ma[i] <-Re(eigen(bigmatrix(fbdem_ma[i,])$MPMmat)$value[1])
}

lambda_ma_m <-Re(eigen(bigmatrix(fbdem_ma_m)$MPMmat)$value[1])

# limit below
fbdem_mb<- as.matrix(FBdem_samp_mb)
fbdem_mb_m <- apply(fbdem_mb, 2, median)

lambda_mb <- 0
for(i in 1:1000){
  lambda_mb[i] <-Re(eigen(bigmatrix(fbdem_mb[i,])$MPMmat)$value[1])
}

lambda_mb_m <-Re(eigen(bigmatrix(fbdem_mb_m)$MPMmat)$value[1])

# limit both
fbdem_mab <- as.matrix(FBdem_samp_mab)
fbdem_mab_m <- apply(fbdem_mab, 2, median)

lambda_mab <- 0
for(i in 1:1000){
  lambda_mab[i] <-Re(eigen(bigmatrix(fbdem_mab[i,])$MPMmat)$value[1])
}

lambda_mab_m <-Re(eigen(bigmatrix(fbdem_mab_m)$MPMmat)$value[1])

# novel control
fbdem_lc<- as.matrix(FBdem_samp_lc)
fbdem_lc_m <- apply(fbdem_lc, 2, median)

lambda_lc <- 0
for(i in 1:1000){
  lambda_lc[i] <-Re(eigen(bigmatrix(fbdem_lc[i,])$MPMmat)$value[1])
}

lambda_lc_m <-Re(eigen(bigmatrix(fbdem_lc_m)$MPMmat)$value[1])

# novel above
fbdem_la<- as.matrix(FBdem_samp_la)
fbdem_la_m <- apply(fbdem_la, 2, median)

lambda_la <- 0
for(i in 1:1000){
  lambda_la[i] <-Re(eigen(bigmatrix(fbdem_la[i,])$MPMmat)$value[1])
}

lambda_la_m <-Re(eigen(bigmatrix(fbdem_la_m)$MPMmat)$value[1])

# novel below
fbdem_lb<- as.matrix(FBdem_samp_lb)
fbdem_lb_m <- apply(fbdem_lb, 2, median)

lambda_lb <- 0
for(i in 1:1000){
  lambda_lb[i] <-Re(eigen(bigmatrix(fbdem_lb[i,])$MPMmat)$value[1])
}

lambda_lb_m <-Re(eigen(bigmatrix(fbdem_lb_m)$MPMmat)$value[1])

# novel both
fbdem_lab <- as.matrix(FBdem_samp_lab)
fbdem_lab_m <- apply(fbdem_lab, 2, median)

lambda_lab<- 0
for(i in 1:1000){
  lambda_lab[i] <-Re(eigen(bigmatrix(fbdem_lab[i,])$MPMmat)$value[1])
}

lambda_lab_m <-Re(eigen(bigmatrix(fbdem_lab_m)$MPMmat)$value[1])


# assemble posterior and median lambda datasets
fblambda<- as.data.frame(cbind(lambda_hc,lambda_ha,lambda_hb,lambda_hab,
                               lambda_lc,lambda_la,lambda_lb,lambda_lab,
                               lambda_mc,lambda_ma,lambda_mb,lambda_mab))
fblambda_m <- as.data.frame(cbind(lambda_hc_m,lambda_ha_m,lambda_hb_m,lambda_hab_m,
                                  lambda_lc_m,lambda_la_m,lambda_lb_m,lambda_lab_m,
                                  lambda_mc_m,lambda_ma_m,lambda_mb_m,lambda_mab_m))

# save them
write.csv(fblambda, "fblambda_nat_withExpTreat.csv")
write.csv(fblambda_m, "fblambda_nat_withExpTreat_median.csv")


#################################################################################
## LTRE calculations

# create midpoints between reference core/control and other treatment combinations
ha_midpoint <- (fbdem_ha_m+fbdem_hc_m)/2
hb_midpoint <- (fbdem_hb_m+fbdem_hc_m)/2
hab_midpoint <- (fbdem_hab_m+fbdem_hc_m)/2
mc_midpoint <- (fbdem_mc_m+fbdem_hc_m)/2
ma_midpoint <- (fbdem_ma_m+fbdem_hc_m)/2
mb_midpoint <- (fbdem_mb_m+fbdem_hc_m)/2
mab_midpoint <- (fbdem_mab_m+fbdem_hc_m)/2
lc_midpoint <- (fbdem_lc_m+fbdem_hc_m)/2
la_midpoint <- (fbdem_la_m+fbdem_hc_m)/2
lb_midpoint <- (fbdem_lb_m+fbdem_hc_m)/2
lab_midpoint <- (fbdem_lab_m+fbdem_hc_m)/2

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
sum.ha <- sum(na.omit(c(LTRE.sens_ha[1:10]*(fbdem_ha_m[1:10]-fbdem_hc_m[1:10]))))
sum.hb <- sum(na.omit(c(LTRE.sens_hb[1:10]*(fbdem_hb_m[1:10]-fbdem_hc_m[1:10]))))
sum.hab <- sum(na.omit(c(LTRE.sens_hab[1:10]*(fbdem_hab_m[1:10]-fbdem_hc_m[1:10]))))
sum.lc <- sum(na.omit(c(LTRE.sens_lc[1:10]*(fbdem_lc_m[1:10]-fbdem_hc_m[1:10]))))
sum.la <- sum(na.omit(c(LTRE.sens_la[1:10]*(fbdem_la_m[1:10]-fbdem_hc_m[1:10]))))
sum.lb <- sum(na.omit(c(LTRE.sens_lb[1:10]*(fbdem_lb_m[1:10]-fbdem_hc_m[1:10]))))
sum.lab <- sum(na.omit(c(LTRE.sens_lab[1:10]*(fbdem_lab_m[1:10]-fbdem_hc_m[1:10]))))
sum.mc <- sum(na.omit(c(LTRE.sens_mc[1:10]*(fbdem_mc_m[1:10]-fbdem_hc_m[1:10]))))
sum.ma <- sum(na.omit(c(LTRE.sens_ma[1:10]*(fbdem_ma_m[1:10]-fbdem_hc_m[1:10]))))
sum.mb <- sum(na.omit(c(LTRE.sens_mb[1:10]*(fbdem_mb_m[1:10]-fbdem_hc_m[1:10]))))
sum.mab <- sum(na.omit(c(LTRE.sens_mab[1:10]*(fbdem_mab_m[1:10]-fbdem_hc_m[1:10]))))

# calculate the LTREs
ltre.ha <- (na.omit(c(LTRE.sens_ha[1:10]*(fbdem_ha_m[1:10]-fbdem_hc_m[1:10]))))
ltre.hb <- (na.omit(c(LTRE.sens_hb[1:10]*(fbdem_hb_m[1:10]-fbdem_hc_m[1:10]))))
ltre.hab <- (na.omit(c(LTRE.sens_hab[1:10]*(fbdem_hab_m[1:10]-fbdem_hc_m[1:10]))))
ltre.lc <- (na.omit(c(LTRE.sens_lc[1:10]*(fbdem_lc_m[1:10]-fbdem_hc_m[1:10]))))
ltre.la <- (na.omit(c(LTRE.sens_la[1:10]*(fbdem_la_m[1:10]-fbdem_hc_m[1:10]))))
ltre.lb <- (na.omit(c(LTRE.sens_lb[1:10]*(fbdem_lb_m[1:10]-fbdem_hc_m[1:10]))))
ltre.lab <- (na.omit(c(LTRE.sens_lab[1:10]*(fbdem_lab_m[1:10]-fbdem_hc_m[1:10]))))
ltre.mc <- (na.omit(c(LTRE.sens_mc[1:10]*(fbdem_mc_m[1:10]-fbdem_hc_m[1:10]))))
ltre.ma <- (na.omit(c(LTRE.sens_ma[1:10]*(fbdem_ma_m[1:10]-fbdem_hc_m[1:10]))))
ltre.mb <- (na.omit(c(LTRE.sens_mb[1:10]*(fbdem_mb_m[1:10]-fbdem_hc_m[1:10]))))
ltre.mab <- (na.omit(c(LTRE.sens_mab[1:10]*(fbdem_mab_m[1:10]-fbdem_hc_m[1:10]))))

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

# repeat the species names
LTREresults$species <- rep("Festuca brachyphylla", 11)

# add in exclosure treatments
LTREresults$treat <- c("above", "below", "both", "control", "above", "below", "both",
                       "control", "above", "below", "both")

# save  the LTRE results data
write.csv(LTREresults, "LTREres_fb_ungherb.csv")

