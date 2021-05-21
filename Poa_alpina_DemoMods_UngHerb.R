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
padat <- demodat[demodat$species=="POAL",]

### vital rate calculations

# Parameterize the growth model for PA

# remove NA's in PA data
padat2 <- padat[!is.na(padat$tillernumb1 & padat$tillernumb2),]
padat2 <- padat2[padat2$tillernumb2<30,] # exclude the one huge individual

tillernumb1 <- padat2$tillernumb1 # tiller number at in year t
tillernumb2 <- padat2$tillernumb2 # tiller number at at year t+1

year <- padat2$year # for random effects
N <- as.numeric(length(padat2$tillernumb1)) # get length of dataset

# set up data and parameters to use
jags.data <- list("N","tillernumb1","tillernumb2", "year")
jags.param <- c( "bPAg","r","a", "prec1", "rss", "rss.new")

# the model
pagrowthmod <- function(){
  for(i in 1:3){y[i]~dnorm(0, prec1)} # year random effect
  for (i in 1:N){
    tillernumb2[i]~dnegbin(p[i],r)
    log(mu[i]) <- a+bPAg*log(tillernumb1[i])+y[year[i]]
    p[i] <- r/(r+mu[i])
    
    res[i] <- pow((tillernumb2[i]-mu[i]),2) # posterior predictive check
    tiller.new[i]~dnegbin(p[i],r)
    res.new[i] <- pow((tiller.new[i]-mu[i]),2)
  }
  # priors
  r~ dunif(0,50)
  a~ dnorm(0,1.0E-6)
  bPAg~dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  # derived parameters
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
pagrowth<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                n.iter=500000 ,model.file=pagrowthmod, n.thin=5, n.chains=3)

pagrowth # results

# post predictive check
pagrowth.paramlist <- pagrowth$BUGSoutput$sims.list

plot(pagrowth.paramlist$rss, pagrowth.paramlist$rss.new, main="A. growth",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(pagrowth.paramlist$rss> pagrowth.paramlist$rss.new) 

# save 95% Cred
bPAg <- sort(pagrowth.paramlist$bPAg) [3750:146250]
aPAg <- sort(pagrowth.paramlist$a) [3750:146250]
rPAg <- sort(pagrowth.paramlist$r) [3750:146250]

pagrowth_post <- as.data.frame(cbind(bPAg,aPAg,rPAg))
write.csv(pagrowth_post, "pagrowth_post.csv")


### survival vital rates for PA

padat3 <- padat[!is.na(padat$survto2),] #subset out individuals still dead
padat3 <- padat3[padat3$tillernumb1<30,] 

# dependent variable 
survto2 <- padat3$survto2
N <- as.numeric(length(padat3$survto2))

# independent variable
tillernumb1 <- padat3$tillernumb1

# for random effect
year <- padat3$year

# set up data and parameters to use
jags.data <- list("N","tillernumb1","survto2", "year")
jags.param <- c( "a","bPAs", "rss", "rss.new", "y", "prec2")

# the model 
pasurvmod <- function(){
  for(j in 1:3){y[j]~dnorm(0, prec2)} # year random effect
  for (i in 1:N){
    survto2[i]~dbern(mu[i])
    logit(mu[i]) <- a + bPAs*log(tillernumb1[i])+y[year[i]]
    
    res[i] <- (survto2[i]-mu[i])^2 # posterior predictive checks
    surv.new[i]~dbern(mu[i])
    res.new[i] <- (surv.new[i]-mu[i])^2
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bPAs~ dnorm(0,1.0E-6)
  prec2~dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
pasurv<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
              n.iter=100000 ,model.file=pasurvmod, n.thin=5, n.chains=3)

pasurv # results


#posterior predictive check
pasurv.paramlist <- pasurv$BUGSoutput$sims.list

plot(pasurv.paramlist$rss,pasurv.paramlist$rss.new, main="B. survival",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(pasurv.paramlist$rss>pasurv.paramlist$rss.new) 

# get 95 % cred of vital rate params
bPAs <- sort(pasurv.paramlist$bPAs) [750:29250]
aPAs <- sort(pasurv.paramlist$a) [750:29250]

pasurv_post <- as.data.frame(cbind(bPAs,aPAs))
write.csv(pasurv_post, "pasurv_post.csv")


## probability of flowering vital rates

# create a 0,1 column for flowering
demoflow <- demodat[!is.na(demodat$flowernumb),]
did_flower <- 0
for(i in 1:length(demoflow$flowernumb)){
  did_flower[i]<- if(demoflow$flowernumb[i] > 0) 1 else 0
}
demoflow$did_flower<- did_flower

# sort to PA
padat4 <- demoflow[demoflow$species=="PA",]
padat4 <- padat4[padat4$tillernumb1<30,] 

# dependent variable
did_flow <- padat4$did_flow
N <- as.numeric(length(padat4$did_flower))

# random effect 
year <- padat4$year

# independent variable- tiller number is year of flowering
tillernumb1 <- padat4$tillernumb1

# set up data and parameters to use
jags.data <- list("N","tillernumb1","did_flow", "year")
jags.param <- c( "a","bPAdf","rss", "rss.new", "prec1", "y")

# the model
padidflowmod <- function(){
  for(j in 1:4){y[j]~dnorm(0,prec1)} # random effect of years
  for (i in 1:N){
    did_flow[i]~dbern(mu[i])
    logit(mu[i]) <- a + bPAdf*log(tillernumb1[i])+y[year[i]]
    
    res[i] <- (did_flow[i]-mu[i])^2 # posterior predictive check
    flow.new[i]~dbern(mu[i])
    res.new[i] <- (flow.new[i]-mu[i])^2
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bPAdf~ dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
padidflow <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                 n.iter=50000 ,model.file=padidflowmod, n.thin=5, n.chains=3)

padidflow # results

#posterior predictive check
padidflow.paramlist <- padidflow$BUGSoutput$sims.list

plot(padidflow.paramlist$rss,padidflow.paramlist$rss.new,
     main="C. did flower?", xlab="SSQ observed",
     ylab="SSQ simulated")
abline(0,1)

# bayesian p- value
mean(padidflow.paramlist$rss>padidflow.paramlist$rss.new) 

# take params from 95% cred
bPAdf <- sort(padidflow.paramlist$bPAdf) [375:14625]
aPAdf <- sort(padidflow.paramlist$a) [375:14625]

#save
padidflow_post <- as.data.frame(cbind(bPAdf,aPAdf))
write.csv(padidflow_post, "padidflow_post.csv")


### if flowering, how many flowers predicted by size
# get the proper data
padat5 <- padat[!padat$flowernumb==0 ,]
padat5 <- padat5[!is.na(padat5$flowernumb),]
padat5 <- padat5[padat5$tillernumb1<30,]

# independent variable
tillernumb1 <- padat5$tillernumb1

# random effect 
year <- padat5$year

# dependent variables
flowernumb <- log(padat5$flowernumb)
N <- as.numeric(length(padat5$flowernumb))

# set up data and parameters to use
jags.data <- list("N","flowernumb","tillernumb1", "year" )
jags.param <- c( "a","bPAfn","rss","rss.new","y", "prec1","prec2")

# the model 
paflownumbmod <- function(){
  for(j in 1:4){y[j]~dnorm(0,prec1)} # year random effect
  for (i in 1:N){
    flowernumb[i]~dnorm(mu[i],prec2)
    mu[i] <- a + bPAfn*log(tillernumb1[i])+y[year[i]]
    
    res[i] <- (flowernumb[i]-mu[i])^2 # posterior predictive check
    flnu.new[i]~dnorm(mu[i],prec2)
    res.new[i] <- (flnu.new[i]-mu[i])^2
  }
  # priors
  a~ dnorm(0,1.0E-6)
  bPAfn~dnorm(0,1.0E-6)
  prec1~dgamma(0.001,0.001)
  prec2~dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# run through JAGS
paflownumb<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                  n.iter=50000 ,model.file=paflownumbmod, n.thin=5, n.chains=3)

paflownumb# results

# posterior check predictive checks
paflownumb.paramlist <- paflownumb$BUGSoutput$sims.list

plot(paflownumb.paramlist$rss,paflownumb.paramlist$rss.new,main="D. flower number",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(paflownumb.paramlist$rss>paflownumb.paramlist$rss.new) 

# 95% cred of vital rate params
bPAfn <- sort(paflownumb.paramlist$bPAfn) [375:14625]
aPAfn <- sort(paflownumb.paramlist$a) [375:14625]

# save them 
paflownumb_post <- as.data.frame(cbind(bPAfn,aPAfn))
write.csv(paflownumb_post, "paflownumb_post.csv")


# recruitment as a function of seed production
# this data is on github
recdat <- read.csv("data_lives_here/AlpDemo_recruit.csv")
parecr <- recdat[recdat$species=="POAL",]

# PA recruitment
padat6 <- parecr[! is.na(parecr$recruits_y2) ,]

# dependent variable
recruits_y2 <- padat6$recruits_y2
N <- as.numeric(length(padat6$recruits_y2))

# for binomial, number of trials (seeds produced)
seed <- round(padat6$seed_numb)

# set up data and parameters to use
jags.data <- list("N","recruits_y2","seed")
jags.param <- c( "aPAr", "rss", "rss.new")

# the model
parecruit <- function(){
  for (i in 1:N){
    recruits_y2[i]~dbinom(mu[i],seed[i])
    logit(mu[i]) <- aPAr 
    
    res[i] <- (recruits_y2[i]-mu[i])^2 # for posterior predictive check
    rec.new[i]~dbinom(mu[i],seed[i])
    res.new[i] <- (rec.new[i]-mu[i])^2
  }
  # priors
  aPAr~ dnorm(0,1.0E-6)
  prec1~ dgamma(0.001,0.001)
  rss <- sum(res[])
  rss.new <- sum(res.new[]) 
}

# run through jags
parec<- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
             n.iter=4000000 ,model.file=parecruit, n.thin=5, n.chains=3)

parec # results

# posterior check
parec.paramlist <- parec$BUGSoutput$sims.list

plot(parec.paramlist$rss,parec.paramlist$rss.new, main="E. recruitment",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

# Bayesian p-value
mean(parec.paramlist$rss>parec.paramlist$rss.new) 

# 95% cred of recruit parm
aPAr <- sort(parec.paramlist$aPAr)[30000:1170000]

# save
parec_post <- as.data.frame(cbind(aPAr))
write.csv(parec_post, "parec_post.csv")

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
g_post <- read.csv("pagrowth_post.csv")
s_post <- read.csv("pasurv_post.csv")
df_post <- read.csv("padidflow_post.csv")
fn_post <- read.csv("paflownumb_post.csv")
r_post <- read.csv("parec_post.csv")

# assemble posterior data frames
# naming goes "parameterSPECIESvitalrate"- e.g., "a" for intercept, "PA" for Poa alpina, and "s" for survival
# draw randomly 1000 times from 95% posterior
aPAs <- sample(s_post$aPAs, 1000)
bPAs <- sample(s_post$bPAs, 1000)
aPAg <- sample(g_post$aPAg, 1000)
bPAg <- sample(g_post$bPAg, 1000)
rPAg <- sample(g_post$rPAg, 1000)
aPAdf <- sample(df_post$aPAdf, 1000)
bPAdf <- sample(df_post$bPAdf, 1000)
aPAfn <- sample(fn_post$aPAfn, 1000)
bPAfn <- sample(fn_post$bPAfn, 1000)
aPAr <- sample(r_post$aPAr,1000)
min_siz <- rep(1, 1000)
max_size <- rep(80, 1000)
seed_flower <- 26.29880952 # estimated number of seed per flower 

PAdem_samp <- as.data.frame(cbind(aPAs,bPAs,aPAg, bPAg,rPAg,aPAdf,bPAdf,aPAfn,bPAfn, 
                                  aPAr, min_siz, max_size,seed_flower))

# Above is the a base control representing core site with no cages 
# now use treatment effects to alter intercepts according to treatment effects in main analysis. 

# load in means of biomass, survival, and flowering. 
s_treat <- read.csv("data_lives_here/surv_postMeans95_ungherb.csv")
g_treat <- read.csv("data_lives_here/Bio_postMeans95_ungherb.csv")
fn_treat <- read.csv("data_lives_here/flow_postMeans95_ungherb.csv")

# filter to species and get effect sizes using high controls as reference
s_treat <- s_treat %>% filter(species=="Poa alpina") 
s_treat$eff <- s_treat$post/s_treat$post[10]-1 # s_treat$post[10] needs to be core control

g_treat <- g_treat %>% filter(species=="Poa alpina") 
g_treat$eff <- (g_treat$post/ g_treat$post[10])-1 # g_treat$post[10] needs to be core control

fn_treat <- fn_treat %>% filter(species=="Poa alpina") 
fn_treat$eff <- fn_treat$post/fn_treat$post[10]-1 # fn_treat$post[10] needs to be core control


# begin multiplying through the effect sizes
# only alter intercepts where we have treatment effects
# check that eff corresponds to desired treatments

# core above
aPAs_ha <- aPAs+(abs(mean(aPAs))*s_treat$eff[1])
aPAg_ha <- aPAg+abs(mean(aPAg))*g_treat$eff[1]
aPAfn_ha <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[1]

PAdem_samp_ha <- as.data.frame(cbind(aPAs_ha,bPAs,aPAg_ha, bPAg,rPAg,aPAdf,bPAdf,aPAfn_ha,bPAfn, 
                                  aPAr, min_siz, max_size,seed_flower))

# core below
aPAs_hb <- aPAs+abs(mean(aPAs))*s_treat$eff[4]
aPAg_hb <- aPAg+abs(mean(aPAg))*g_treat$eff[4]
aPAfn_hb <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[4]

PAdem_samp_hb <- as.data.frame(cbind(aPAs_hb,bPAs,aPAg_hb, bPAg,rPAg,aPAdf,bPAdf,aPAfn_hb,bPAfn, 
                                     aPAr, min_siz, max_size,seed_flower))

# core both
aPAs_hab <- aPAs+abs(mean(aPAs))*s_treat$eff[7]
aPAg_hab <- aPAg+abs(mean(aPAg))*g_treat$eff[7]
aPAfn_hab <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[7]

PAdem_samp_hab <- as.data.frame(cbind(aPAs_hab,bPAs,aPAg_hab, bPAg,rPAg,aPAdf,bPAdf,
                                      aPAfn_hab,bPAfn, aPAr, min_siz, max_size,seed_flower))

# limit control
aPAs_mc <- aPAs+abs(mean(aPAs))*s_treat$eff[11]
aPAg_mc <- aPAg+abs(mean(aPAg))*g_treat$eff[11]
aPAfn_mc <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[11]

PAdem_samp_mc <- as.data.frame(cbind(aPAs_mc,bPAs,aPAg_mc, bPAg,rPAg,aPAdf,bPAdf,aPAfn_mc,bPAfn, 
                                     aPAr, min_siz, max_size,seed_flower))

# limit above
aPAs_ma <- aPAs+abs(mean(aPAs))*s_treat$eff[2]
aPAg_ma <- aPAg+abs(mean(aPAg))*g_treat$eff[2]
aPAfn_ma <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[2]

PAdem_samp_ma <- as.data.frame(cbind(aPAs_ma,bPAs,aPAg_ma, bPAg,rPAg,aPAdf,bPAdf,aPAfn_ma,bPAfn, 
                                     aPAr, min_siz, max_size,seed_flower))

# limit below
aPAs_mb <- aPAs+abs(mean(aPAs))*s_treat$eff[5]
aPAg_mb <- aPAg+abs(mean(aPAg))*g_treat$eff[5]
aPAfn_mb <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[5]

PAdem_samp_mb <- as.data.frame(cbind(aPAs_mb,bPAs,aPAg_mb, bPAg,rPAg,aPAdf,bPAdf,aPAfn_mb,bPAfn, 
                                     aPAr, min_siz, max_size,seed_flower))

# limit both
aPAs_mab <- aPAs+abs(mean(aPAs))*s_treat$eff[8]
aPAg_mab <- aPAg+abs(mean(aPAg))*g_treat$eff[8]
aPAfn_mab <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[8]

PAdem_samp_mab <- as.data.frame(cbind(aPAs_mab,bPAs,aPAg_mab, bPAg,rPAg,aPAdf,bPAdf,
                                      aPAfn_mab,bPAfn,aPAr, min_siz, max_size,seed_flower))

# novelsite control
aPAs_lc <- aPAs+abs(mean(aPAs))*s_treat$eff[12]
aPAg_lc <- aPAg+abs(mean(aPAg))*g_treat$eff[12]
aPAfn_lc <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[12]

PAdem_samp_lc <- as.data.frame(cbind(aPAs_lc,bPAs,aPAg_lc, bPAg,rPAg,aPAdf,bPAdf,aPAfn_lc,bPAfn, 
                                     aPAr, min_siz, max_size,seed_flower))

# novel above
aPAs_la <- aPAs+abs(mean(aPAs))*s_treat$eff[3]
aPAg_la <- aPAg+abs(mean(aPAg))*g_treat$eff[3]
aPAfn_la <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[3]

PAdem_samp_la <- as.data.frame(cbind(aPAs_la,bPAs,aPAg_la, bPAg,rPAg,aPAdf,bPAdf,aPAfn_la,bPAfn, 
                                     aPAr, min_siz, max_size,seed_flower))

# novel below
aPAs_lb <- aPAs+abs(mean(aPAs))*s_treat$eff[6]
aPAg_lb <- aPAg+abs(mean(aPAg))*g_treat$eff[6]
aPAfn_lb <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[6]

PAdem_samp_lb <- as.data.frame(cbind(aPAs_lb,bPAs,aPAg_lb, bPAg,rPAg,aPAdf,bPAdf,aPAfn_lb,bPAfn, 
                                     aPAr, min_siz, max_size,seed_flower))

# novel both
aPAs_lab <- aPAs+abs(mean(aPAs))*s_treat$eff[9]
aPAg_lab <- aPAg+abs(mean(aPAg))*g_treat$eff[9]
aPAfn_lab <- aPAfn+abs(mean(aPAfn))*fn_treat$eff[9]

PAdem_samp_lab <- as.data.frame(cbind(aPAs_lab,bPAs,aPAg_lab, bPAg,rPAg,aPAdf,bPAdf,
                                      aPAfn_lab,bPAfn,aPAr, min_siz, max_size,seed_flower))


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
padem_hc <- as.matrix(PAdem_samp) # full posterior draws
padem_hc_m <- apply(padem_hc, 2, median) # median of each parameter after draws

lambda_hc <- 0
for(i in 1:1000){
  lambda_hc[i] <-Re(eigen(bigmatrix(padem_hc[i,])$MPMmat)$value[1])
}

hist(lambda_hc);mean(lambda_hc);max(lambda_hc);sd(lambda_hc)
lambda_hc_m <- Re(eigen(bigmatrix(padem_hc_m)$MPMmat)$value[1])

# core above
padem_ha <- as.matrix(PAdem_samp_ha)
padem_ha_m <- apply(padem_ha, 2, median) 

lambda_ha <- 0
for(i in 1:1000){
  lambda_ha[i] <-Re(eigen(bigmatrix(padem_ha[i,])$MPMmat)$value[1])
}

lambda_ha_m <-Re(eigen(bigmatrix(padem_ha_m)$MPMmat)$value[1])

# core below
padem_hb <- as.matrix(PAdem_samp_hb)
padem_hb_m <- apply(padem_hb, 2, median)

lambda_hb <- 0
for(i in 1:1000){
  lambda_hb[i] <-Re(eigen(bigmatrix(padem_hb[i,])$MPMmat)$value[1])
}

lambda_hb_m <-Re(eigen(bigmatrix(padem_hb_m)$MPMmat)$value[1])

# core both
padem_hab <- as.matrix(PAdem_samp_hab)
padem_hab_m <- apply(padem_hab, 2, median)

lambda_hab <- 0
for(i in 1:1000){
  lambda_hab[i] <-Re(eigen(bigmatrix(padem_hab[i,])$MPMmat)$value[1])
}

lambda_hab_m <-Re(eigen(bigmatrix(padem_hab_m)$MPMmat)$value[1])

# limit control
padem_mc <- as.matrix(PAdem_samp_mc)
padem_mc_m <- apply(padem_mc, 2, median)

lambda_mc <- 0
for(i in 1:1000){
  lambda_mc[i] <-Re(eigen(bigmatrix(padem_mc[i,])$MPMmat)$value[1])
}

lambda_mc_m <-Re(eigen(bigmatrix(padem_mc_m)$MPMmat)$value[1])

# limit above
padem_ma <- as.matrix(PAdem_samp_ma)
padem_ma_m <- apply(padem_ma, 2, median)

lambda_ma <- 0
for(i in 1:1000){
  lambda_ma[i] <-Re(eigen(bigmatrix(padem_ma[i,])$MPMmat)$value[1])
}

lambda_ma_m <-Re(eigen(bigmatrix(padem_ma_m)$MPMmat)$value[1])

# limit below
padem_mb <- as.matrix(PAdem_samp_mb)
padem_mb_m <- apply(padem_mb, 2, median)

lambda_mb <- 0
for(i in 1:1000){
  lambda_mb[i] <-Re(eigen(bigmatrix(padem_mb[i,])$MPMmat)$value[1])
}

lambda_mb_m <-Re(eigen(bigmatrix(padem_mb_m)$MPMmat)$value[1])

# limit both
padem_mab <- as.matrix(PAdem_samp_mab)
padem_mab_m <- apply(padem_mab, 2, median)

lambda_mab <- 0
for(i in 1:1000){
  lambda_mab[i] <-Re(eigen(bigmatrix(padem_mab[i,])$MPMmat)$value[1])
}

lambda_mab_m <-Re(eigen(bigmatrix(padem_mab_m)$MPMmat)$value[1])

# novel control
padem_lc <- as.matrix(PAdem_samp_lc)
padem_lc_m <- apply(padem_lc, 2, median)

lambda_lc <- 0
for(i in 1:1000){
  lambda_lc[i] <-Re(eigen(bigmatrix(padem_lc[i,])$MPMmat)$value[1])
}

lambda_lc_m <-Re(eigen(bigmatrix(padem_lc_m)$MPMmat)$value[1])

# novel above
padem_la <- as.matrix(PAdem_samp_la)
padem_la_m <- apply(padem_la, 2, median)

lambda_la <- 0
for(i in 1:1000){
  lambda_la[i] <-Re(eigen(bigmatrix(padem_la[i,])$MPMmat)$value[1])
}

lambda_la_m <-Re(eigen(bigmatrix(padem_la_m)$MPMmat)$value[1])

# novel below
padem_lb <- as.matrix(PAdem_samp_lb)
padem_lb_m <- apply(padem_lb, 2, median)

lambda_lb <- 0
for(i in 1:1000){
  lambda_lb[i] <-Re(eigen(bigmatrix(padem_lb[i,])$MPMmat)$value[1])
}

lambda_lb_m <-Re(eigen(bigmatrix(padem_lb_m)$MPMmat)$value[1])

# novel both
padem_lab <- as.matrix(PAdem_samp_lab)
padem_lab_m <- apply(padem_lab, 2, median)

lambda_lab <- 0
for(i in 1:1000){
  lambda_lab[i] <-Re(eigen(bigmatrix(padem_lab[i,])$MPMmat)$value[1])
}

lambda_lab_m <-Re(eigen(bigmatrix(padem_lab_m)$MPMmat)$value[1])


# assemble posterior and median lambda datasets
palambda<- as.data.frame(cbind(lambda_hc,lambda_ha,lambda_hb,lambda_hab,
                               lambda_lc,lambda_la,lambda_lb,lambda_lab,
                               lambda_mc,lambda_ma,lambda_mb,lambda_mab))

palambda_m <- as.data.frame(cbind(lambda_hc_m,lambda_ha_m,lambda_hb_m,lambda_hab_m,
                               lambda_lc_m,lambda_la_m,lambda_lb_m,lambda_lab_m,
                               lambda_mc_m,lambda_ma_m,lambda_mb_m,lambda_mab_m))

write.csv(palambda, "palambda_nat_withExpTreat.csv")
write.csv(palambda_m, "palambda_nat_withExpTreat_median.csv")

#################################################################################
## LTRE calculations

# create midpoints between reference core/control and other treatment combinations
ha_midpoint <- (padem_ha_m+padem_hc_m)/2
hb_midpoint <- (padem_hb_m+padem_hc_m)/2
hab_midpoint <- (padem_hab_m+padem_hc_m)/2
mc_midpoint <- (padem_mc_m+padem_hc_m)/2
ma_midpoint <- (padem_ma_m+padem_hc_m)/2
mb_midpoint <- (padem_mb_m+padem_hc_m)/2
mab_midpoint <- (padem_mab_m+padem_hc_m)/2
lc_midpoint <- (padem_lc_m+padem_hc_m)/2
la_midpoint <- (padem_la_m+padem_hc_m)/2
lb_midpoint <- (padem_lb_m+padem_hc_m)/2
lab_midpoint <- (padem_lab_m+padem_hc_m)/2

# Midpoint lambdas
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
sum.ha <- sum(na.omit(c(LTRE.sens_ha[1:10]*(padem_ha_m[1:10]-padem_hc_m[1:10]))))
sum.hb <- sum(na.omit(c(LTRE.sens_hb[1:10]*(padem_hb_m[1:10]-padem_hc_m[1:10]))))
sum.hab <- sum(na.omit(c(LTRE.sens_hab[1:10]*(padem_hab_m[1:10]-padem_hc_m[1:10]))))
sum.lc <- sum(na.omit(c(LTRE.sens_lc[1:10]*(padem_lc_m[1:10]-padem_hc_m[1:10]))))
sum.la <- sum(na.omit(c(LTRE.sens_la[1:10]*(padem_la_m[1:10]-padem_hc_m[1:10]))))
sum.lb <- sum(na.omit(c(LTRE.sens_lb[1:10]*(padem_lb_m[1:10]-padem_hc_m[1:10]))))
sum.lab <- sum(na.omit(c(LTRE.sens_lab[1:10]*(padem_lab_m[1:10]-padem_hc_m[1:10]))))
sum.mc <- sum(na.omit(c(LTRE.sens_mc[1:10]*(padem_mc_m[1:10]-padem_hc_m[1:10]))))
sum.ma <- sum(na.omit(c(LTRE.sens_ma[1:10]*(padem_ma_m[1:10]-padem_hc_m[1:10]))))
sum.mb <- sum(na.omit(c(LTRE.sens_mb[1:10]*(padem_mb_m[1:10]-padem_hc_m[1:10]))))
sum.mab <- sum(na.omit(c(LTRE.sens_mab[1:10]*(padem_mab_m[1:10]-padem_hc_m[1:10]))))

# calculate the LTREs
ltre.ha <- (na.omit(c(LTRE.sens_ha[1:10]*(padem_ha_m[1:10]-padem_hc_m[1:10]))))
ltre.hb <- (na.omit(c(LTRE.sens_hb[1:10]*(padem_hb_m[1:10]-padem_hc_m[1:10]))))
ltre.hab <- (na.omit(c(LTRE.sens_hab[1:10]*(padem_hab_m[1:10]-padem_hc_m[1:10]))))
ltre.lc <- (na.omit(c(LTRE.sens_lc[1:10]*(padem_lc_m[1:10]-padem_hc_m[1:10]))))
ltre.la <- (na.omit(c(LTRE.sens_la[1:10]*(padem_la_m[1:10]-padem_hc_m[1:10]))))
ltre.lb <- (na.omit(c(LTRE.sens_lb[1:10]*(padem_lb_m[1:10]-padem_hc_m[1:10]))))
ltre.lab <- (na.omit(c(LTRE.sens_lab[1:10]*(padem_lab_m[1:10]-padem_hc_m[1:10]))))
ltre.mc <- (na.omit(c(LTRE.sens_mc[1:10]*(padem_mc_m[1:10]-padem_hc_m[1:10]))))
ltre.ma <- (na.omit(c(LTRE.sens_ma[1:10]*(padem_ma_m[1:10]-padem_hc_m[1:10]))))
ltre.mb <- (na.omit(c(LTRE.sens_mb[1:10]*(padem_mb_m[1:10]-padem_hc_m[1:10]))))
ltre.mab <- (na.omit(c(LTRE.sens_mab[1:10]*(padem_mab_m[1:10]-padem_hc_m[1:10]))))

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
LTREresults$species <- rep("Poa alpina", 11)

# add in exclosure treatments 
LTREresults$treat <- c("above", "below", "both", "control", "above", "below", "both",
                       "control", "above", "below", "both")

# save  the LTRE results data
write.csv(LTREresults, "LTREres_pa_ungherb.csv")

