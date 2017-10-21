##################################################################
# Model 1 -- hierarchical model of SNARC effect
##################################################################


library(tidyverse)
library(R2jags)
library(polspline)

dat <- read_csv("data.csv")

# look at RTs
dat %>%
  filter(error==0 & RT<3000) %>%
  ggplot(aes(x=RT)) +
  geom_density(aes(fill=parity), alpha=0.5)

filtered = dat %>%
  filter(error==0 & RT<3000)
  
data <- dat %>%
  filter(error==0 & RT<3000) %>%
  group_by(subject,stimulus,targetHand) %>%
  summarize(medRT=median(RT)) %>%
  mutate(dRT=medRT-lag(medRT)) %>%
  filter(!is.na(dRT)) %>%
  select(subject,stimulus,dRT)


#########################################################
# build model

ns<-length(unique(data$subject))
dRT=matrix(data$dRT,byrow=TRUE,nrow=length(unique(data$subject)),ncol=length(unique(data$stimulus)))
num=c(1,2,8,9)
jagsData <- list("dRT", "num", "ns")

# collect samples from posterior distributions
set.seed(315)
samples <- jags(jagsData, 
                 inits=NULL, 
                 parameters.to.save=c("B", "B.var", "sigma"), 
                 model.file="modelSnarc.txt", 
                 n.chains=3, 
                 n.iter=1e5,
                 n.burnin=5000, 
                 n.thin=10, 
                 DIC=T)

# print summary of posterior samples
samples
traceplot(samples)
posterior.slope = samples$BUGSoutput$sims.matrix[,1]

# 95% HPDI
sampMCMC = as.mcmc(posterior.slope)
HPDI=HPDinterval(sampMCMC, prob=0.95)


# posterior mode
dens=density(posterior.slope)
postMode=dens$x[which.max(dens$y)]

# construct figure for posterior density of slope
op = par(cex.main=1.5, mar=c(5,6,4,5)+0.1, mgp=c(3.5, 1,0), cex.lab=1.5, font.lab=2, cex.axis=1.3, bty="n", las=1)
plot(dens, xlim=c(-22,3), ylim=c(0,0.25), main="", xlab="", ylab="", lwd=2, lty=1, axes=FALSE)
lines(x=c(postMode,postMode), y=c(0.205,0.22), lwd=2, col="grey")

axis(1)
axis(2)
par(las=0)
mtext("Group-level slope b", side=1, line=2.5, cex=1.5)
mtext("Density", side=2, line=3.5, cex=1.5)

text(x=postMode, y=0.235, paste("Posterior mode =",round(postMode,1), sep=" "), cex=1.2)
lines(x=c(HPDI[1],HPDI[1]), y=c(0,0.02), lwd=2, col="grey")
lines(x=c(HPDI[2],HPDI[2]), y=c(0,0.02), lwd=2, col="grey")
lines(x=c(HPDI[1],HPDI[2]), y=c(0.01,0.01), lwd=2, col="grey")
text(x=(HPDI[1]+HPDI[2])/2, y=0.05, "95% HPDI", cex=1.2)
text(x=(HPDI[1]+HPDI[2])/2, y=0.025, paste("[", round(HPDI[1],1),", ",round(HPDI[2],1),"]", sep=""), cex=1.2)


# Savage Dickey density ratio - compute Bayes factor for overall snarc effect
fit.posterior <- logspline(posterior.slope)
posterior <- dlogspline(0, fit.posterior) 
prior     <- dunif(0,-20,20)                 
BF01      <- posterior/prior
BF01    
1/BF01

