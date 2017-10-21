# simulation of SNARC effect
library(polspline)

B = -10 # mean of SNARC effect
sigma.B = 1 # sd of SNARC effect

N=15

set.seed(134) # random seed
numbers = c(1,2,8,9)
ind.slopes = rnorm(N, mean=B, sd=sigma.B)
ind.intercepts = runif(N, min=-200, max=200)

dRT = matrix(numeric(N*4), nrow=N, ncol=4)

for (i in 1:N){
  for (j in 1:4){
    dRT[i,j] = ind.intercepts[i] + ind.slopes[i]*numbers[j] + rnorm(1,0,100)
  }
}

# fit individual regression slopes
slopes=numeric(N)
intercepts=numeric(N)
for (i in 1:N){
  model = lm(dRT[i,]~numbers)
  slopes[i] = model$coefficients[2]
  intercepts[i] = model$coefficients[1]
}



# fit hierarchical Bayesian model
ns<-N
num=numbers
jagsData <- list("dRT", "num", "ns")

# collect samples from posterior distributions
set.seed(132)
samples <- jags(jagsData, 
                inits=NULL, 
                parameters.to.save=c("B"), 
                model.file="modelSnarc.txt", 
                n.chains=3, 
                n.iter=1e5,
                n.burnin=5000, 
                n.thin=10, 
                DIC=T)

posterior.slope = samples$BUGSoutput$sims.matrix[,1]
dens=density(posterior.slope)

sampMCMC = as.mcmc(posterior.slope)
HPDI=HPDinterval(sampMCMC, prob=0.95)

######################################################
# construct figure for posterior density of slope and interval estimates

op = par(cex.main=1.5, mar=c(5,6,4,5)+0.1, mgp=c(3.5, 1,0), cex.lab=1.5, font.lab=2, cex.axis=1.3, bty="n", las=1)
plot(dens, xlim=c(min(slopes)-1,max(slopes)+5), ylim=c(0,0.3), main="", xlab="", ylab="", lwd=2, lty=1, axes=FALSE)

axis(1)
axis(2)
par(las=0)
mtext("Group-level slope b", side=1, line=2.5, cex=1.5)
mtext("Density", side=2, line=3.5, cex=1.5)

par(new=TRUE)
hist(slopes, breaks=20, probability=TRUE, xlim=c(min(slopes)-1,max(slopes)+5),ylim=c(0,0.3), main="", xlab="", ylab="", lty=1, lwd=2, axes=FALSE)

# plot 95% HPDI
lines(x=c(HPDI[1],HPDI[1]), y=c(0.24,0.26), lwd=2, col="grey")
lines(x=c(HPDI[2],HPDI[2]), y=c(0.24,0.26), lwd=2, col="grey")
lines(x=c(HPDI[1],HPDI[2]), y=c(0.25,0.25), lwd=2, col="grey")
text(x=-19, y=0.25, "95% HPDI", cex=1.2, adj=c(1,0.5))

# plot 95% CI
tTest = t.test(slopes, mu=0)
CIleft = tTest$conf.int[1]
CIright = tTest$conf.int[2]

lines(x=c(CIleft,CIleft), y=c(0.19,0.21), lwd=2, col="grey")
lines(x=c(CIright,CIright), y=c(0.19,0.21), lwd=2, col="grey")
lines(x=c(CIleft,CIright), y=c(0.20,0.20), lwd=2, col="grey")
text(x=-19, y=0.20, "95% CI", cex=1.2, adj=c(1,0.5))

lines(x=c(-10,-10), y=c(0,0.3), lwd=2, lty=3, col="grey")


########################
# hypothesis tests
t.test(slopes, mu=0)

fit.posterior <- logspline(posterior.slope)
posterior <- dlogspline(0, fit.posterior) 
prior     <- dunif(0,-20,20)                 
BF01      <- posterior/prior
BF01    
1/BF01
