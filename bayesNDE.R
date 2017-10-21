library(tidyverse)
library(R2jags)
library(polspline) # for computing Savage-Dickey density ratio

RTs = read.csv("fazioRTs.csv",header=TRUE)
err = read.csv("fazioErrrors.csv",header=TRUE)

RTs = RTs[,1:41]
err=err[,1:41]

longRTs = RTs %>%
  gather(key="stimPair", value="rt",-Subject) %>%
  separate(stimPair, into=c("stim1","stim2")) %>%
  separate(stim1, into=c("x","stim1"), sep=1) %>%
  select(-x)

longErr = err %>%
  gather(key="stimPair", value="correct", -Subject) %>%
  separate(stimPair, into=c("stim1", "stim2")) %>%
  separate(stim1, into=c("x","stim1"), sep=1) %>%
  select(-x)

dataMerge = data.frame(longRTs,longErr$correct)
names(dataMerge)[5]="correct"

dat = dataMerge %>%
  mutate(stim1=as.numeric(stim1), stim2=as.numeric(stim2)) %>%
  mutate(distance=abs(stim1 - stim2)) %>%
  mutate(ratio = ifelse(stim1/stim2>1, stim1/stim2, stim2/stim1)) %>%
  mutate(distance=as.factor(distance)) %>%
  filter(correct==1 & rt<5000)

dat$ratioBin = dat$ratio
dat$ratioBin[dat$ratio<1.28] <- 1
dat$ratioBin[dat$ratio>1.28 & dat$ratio<1.43] <- 2
dat$ratioBin[dat$ratio>1.48 & dat$ratio<1.65] <- 3
dat$ratioBin[dat$ratio>2.42] <- 4

# look at overall distance effect
dat %>%
  mutate(ratioBin=as.factor(ratioBin)) %>%
  ggplot(aes(x=ratioBin, y=rt)) +
  geom_boxplot()

dat %>%
  mutate(ratioBin=as.factor(ratioBin)) %>%
  ggplot(aes(x=rt, group=ratioBin)) +
  geom_density(aes(fill=ratioBin), alpha=0.5)

# prepare for modeling
data <- dat %>%
  group_by(Subject,ratioBin) %>%
  summarize(medRT=median(rt)) %>%
  select(Subject,ratioBin,medRT)

###################################
# modeling

ns<-length(unique(data$Subject))
medRT=matrix(data$medRT,byrow=TRUE,nrow=length(unique(data$Subject)),ncol=length(unique(data$ratioBin)))
ratioBin=c(1,2,3,4)
jagsData <- list("medRT", "ratioBin", "ns")

# collect samples from posterior distributions
set.seed(154)
samples <- jags(jagsData, 
                 inits=NULL, 
                 parameters.to.save=c("a", "b", "B"), 
                 model.file="modelDistance.txt", 
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

# construct figure for posterior density
op = par(cex.main=1.5, mar=c(5,6,4,5)+0.1, mgp=c(3.5, 1,0), cex.lab=1.5, font.lab=2, cex.axis=1.3, bty="n", las=1)
plot(dens, xlim=c(-100,-30), ylim=c(0,0.06), main="", xlab="", ylab="", lwd=2, lty=1, axes=FALSE)
lines(x=c(postMode,postMode), y=c(0.054,0.05), lwd=2, col="grey")

axis(1)
axis(2)
par(las=0)
mtext("Group-level slope b", side=1, line=2.5, cex=1.5)
mtext("Density", side=2, line=3.5, cex=1.5)

text(x=postMode, y=0.058, paste("Posterior mode =",round(postMode,1)), cex=1.2)
lines(x=c(HPDI[1],HPDI[1]), y=c(0,0.005), lwd=2, col="grey")
lines(x=c(HPDI[2],HPDI[2]), y=c(0,0.005), lwd=2, col="grey")
lines(x=c(HPDI[1],HPDI[2]), y=c(0.0025,0.0025), lwd=2, col="grey")
text(x=(HPDI[1]+HPDI[2])/2, y=0.014, "95% HPDI", cex=1.2)
text(x=(HPDI[1]+HPDI[2])/2, y=0.008, paste("[",round(HPDI[1],1),", ",round(HPDI[2],1),"]",sep=""), cex=1.2)


# Savage Dickey density ratio - compute Bayes factor for overall NDE
fit.posterior <- logspline(posterior.slope)
posterior <- dlogspline(0, fit.posterior) 
prior     <- dunif(0,-100,100)                 
BF01      <- posterior/prior
BF01    
1/BF01


