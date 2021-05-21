### Allometric equations for biomass

# load in data
dat <- read.csv("data_lives_here/MamExpBioSurvData.csv")
head(dat)

#data for each species
padat <- dat[dat$species=="POAL",]
fbdat <- dat[dat$species=="FEBR",]
esdat <- dat[dat$species=="ELSC",]

# Poa alpina allometric equations
paallo <- lm(live_bio~tiller+I(tiller^2), data=padat)
summary(paallo) # use parameters from here in Fitness mode script

# Festuca brachyphylla allometric equations
fballo <- lm(live_bio~tiller+I(tiller^2)-1, data=fbdat) # removed intercept given it predicted negative biomass
summary(fballo) # use parameters from here in Fitness mode script

# Elymus scribneri allometric equations
esallo <- lm(live_bio~tiller+I(tiller^2)-1, data=esdat)# removed intercept given it predicted negative biomass
summary(esallo) # use parameters from here in Fitness mode script