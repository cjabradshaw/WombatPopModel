##########################################################################################################################################
## common wombat (Vombatus ursinus) demographic model
## 
## Corey Bradshaw
## corey.bradshaw@flinders.edu.au
## Flinders University, September 2021
##########################################################################################################################################

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## source
source("matrixOperators.r")


##############################
## VOMBATUS (ursinus) (VU)
## sources: Roger et al. 2011 Popul Ecol 53:215-227

# mass
VU.mass <- 25 # Vombatus ursinus (Saran et al. 2011 Pacific Conservation Biology 17:310-319)

## predicted rm (from Henneman 1983 Oecologia 56:104-108)
## log10rm = 0.6914 - 0.2622*log10m (mass in g)
VU.rm.pred <- 10^(0.6914 - (0.2622*log10(VU.mass*1000)))
VU.lm.pred <- exp(VU.rm.pred)

## theoretical population density for mammalian herbivores based on body size (Damuth 1981; Freeland 1990)
## log10D = 4.196 − 0.74*(log10m)
VU.D.pred <- (10^(4.196 - (0.74*log10(VU.mass*1000))))/2 # divided by 2 for females only
VU.D.pred # animals/km2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
VU.age.max <- round(10^(0.89 + (0.13*log10(VU.mass*1000))), 0)
VU.age.max <- 26 # reset based on McIlroy 2008

## age vector
VU.age.vec <- 0:VU.age.max

## fertility
## total fecundity from Allainé et al. 1987 (Oecologia)
## lnF = 2.719 - 0.211lnM (all mammals)
VU.F.pred <- exp(2.719 - (0.211*log(VU.mass*1000)))/2 # divided by 2 for females
VU.F.pred <- 1/2/2 # inter-birth interval of 2 years / 2 for daughters only (McIlroy 1995)

## age at primiparity
## lnalpha = 0.214 + 0.263*lnM (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
VU.alpha <- ceiling(exp(-1.34 + (0.214*log(VU.mass*1000))))
VU.alpha <- 2  # reset according to Roger et al. (2011)

## define m function with age
## pouch young per year
VU.pypy <- 0.5

## proportion of females breeding each year
VU.pfbr <- 0.84

## proportion young that are female
VU.sr <- 0.5 

## m vector
VU.m.vec <- c(0, 0.3*VU.F.pred, 0.9*VU.F.pred, rep(VU.F.pred, 24))
VU.m.sd.vec <- 0.05*VU.m.vec
plot(VU.age.vec, VU.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
VU.m.dat <- data.frame(VU.age.vec, VU.m.vec)
param.init <- c(0.2, 2, -3)
VU.fit.logp <- nls(VU.m.vec ~ a / (1+(VU.age.vec/b)^c), 
                   data = VU.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
VU.fit.logp.summ <- summary(VU.fit.logp)
plot(VU.age.vec, VU.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
VU.age.vec.cont <- seq(0,max(VU.age.vec),1)
VU.pred.p.m <- coef(VU.fit.logp)[1] / (1+(VU.age.vec.cont/coef(VU.fit.logp)[2])^coef(VU.fit.logp)[3])
VU.pred.p.mm <- ifelse(VU.pred.p.m > 1, 1, VU.pred.p.m)
VU.pred.p.mm[2] <- 0
lines(VU.age.vec.cont, VU.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -0.5; b.s <- -0.25
VU.s.tran <- ln.a.s + b.s*log(VU.mass*1000) + log(1)
VU.s.ad.yr <- exp(-exp(VU.s.tran))

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (1.05*VU.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 3.1 # rate of mortality decline (also known as bt)
a2 <- 1 - VU.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 0.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.02 # rate of mortality increase
longev <- VU.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
VU.Sx <- c(0.99*VU.s.ad.yr, 1 - qx)
plot(x, VU.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
VU.s.sd.vec <- 0.05*VU.Sx

## create matrix
VU.popmat <- matrix(data = 0, nrow=VU.age.max+1, ncol=VU.age.max+1)
diag(VU.popmat[2:(VU.age.max+1),]) <- VU.Sx[-(VU.age.max+1)]
VU.popmat[VU.age.max+1,VU.age.max+1] <- 0
VU.popmat[1,] <- VU.pred.p.mm
colnames(VU.popmat) <- c(0:VU.age.max)
rownames(VU.popmat) <- c(0:VU.age.max)
VU.popmat.orig <- VU.popmat ## save original matrix

## matrix properties
max.lambda(VU.popmat.orig) ## 1-yr lambda
VU.lm.pred
max.r(VU.popmat.orig) # rate of population change, 1-yr
VU.ssd <- stable.stage.dist(VU.popmat.orig) ## stable stage distribution
plot(VU.age.vec, VU.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(VU.popmat.orig, VU.age.max) # reproductive value
VU.gen.l <- G.val(VU.popmat.orig, VU.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
VU.pop.found <- round(area*VU.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
VU.init.vec <- VU.ssd * VU.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*VU.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

VU.tot.F <- sum(VU.popmat.orig[1,])
VU.popmat <- VU.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
VU.n.mat <- matrix(0, nrow=VU.age.max+1,ncol=(t+1))
VU.n.mat[,1] <- VU.init.vec

## set up projection loop
for (i in 1:t) {
  VU.n.mat[,i+1] <- VU.popmat %*% VU.n.mat[,i]
}

VU.n.pred <- colSums(VU.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(VU.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
VU.K.max <- 1*VU.pop.found
VU.K.vec <- c(1, VU.K.max/2, 0.75*VU.K.max, VU.K.max) 
VU.red.vec <- c(1,0.985,0.951,0.873)
plot(VU.K.vec, VU.red.vec,pch=19,type="b")
VU.Kred.dat <- data.frame(VU.K.vec, VU.red.vec)

# logistic power function a/(1+(x/b)^c)
VU.param.init <- c(1, 2*VU.K.max, 2)
VU.fit.lp <- nls(VU.red.vec ~ a/(1+(VU.K.vec/b)^c), 
                 data = VU.Kred.dat,
                 algorithm = "port",
                 start = c(a = VU.param.init[1], b = VU.param.init[2], c = VU.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
VU.fit.lp.summ <- summary(VU.fit.lp)
plot(VU.K.vec, VU.red.vec, pch=19,xlab="N",ylab="reduction factor")
VU.K.vec.cont <- seq(1,2*VU.pop.found,1)
VU.pred.lp.fx <- coef(VU.fit.lp)[1]/(1+(VU.K.vec.cont/coef(VU.fit.lp)[2])^coef(VU.fit.lp)[3])
lines(VU.K.vec.cont, VU.pred.lp.fx, lty=3,lwd=3,col="red")

VU.a.lp <- coef(VU.fit.lp)[1]
VU.b.lp <- coef(VU.fit.lp)[2]
VU.c.lp <- coef(VU.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
VU.n.mat <- matrix(0, nrow=VU.age.max+1, ncol=(t+1))
VU.n.mat[,1] <- VU.init.vec
VU.popmat <- VU.popmat.orig

## set up projection loop
for (i in 1:t) {
  VU.totN.i <- sum(VU.n.mat[,i])
  VU.pred.red <- as.numeric(VU.a.lp/(1+(VU.totN.i/VU.b.lp)^VU.c.lp))
  diag(VU.popmat[2:(VU.age.max+1),]) <- (VU.Sx[-(VU.age.max+1)])*VU.pred.red
  VU.popmat[VU.age.max+1,VU.age.max+1] <- 0
  VU.popmat[1,] <- VU.pred.p.mm
  VU.n.mat[,i+1] <- VU.popmat %*% VU.n.mat[,i]
}

VU.n.pred <- colSums(VU.n.mat)
plot(yrs, VU.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=VU.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 100
itdiv <- iter/10

VU.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
VU.s.arr <- VU.m.arr <- array(data=NA, dim=c(t+1, VU.age.max+1, iter))

for (e in 1:iter) {
  VU.popmat <- VU.popmat.orig
  
  VU.n.mat <- matrix(0, nrow=VU.age.max+1,ncol=(t+1))
  VU.n.mat[,1] <- VU.init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    VU.s.alpha <- estBetaParams(VU.Sx, VU.s.sd.vec^2)$alpha
    VU.s.beta <- estBetaParams(VU.Sx, VU.s.sd.vec^2)$beta
    VU.s.stoch <- rbeta(length(VU.s.alpha), VU.s.alpha, VU.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    VU.fert.stch <- rnorm(length(VU.popmat[,1]), VU.pred.p.mm, VU.m.sd.vec)
    VU.m.arr[i,,e] <- ifelse(VU.fert.stch < 0, 0, VU.fert.stch)
    
    VU.totN.i <- sum(VU.n.mat[,i], na.rm=T)
    VU.pred.red <- VU.a.lp/(1+(VU.totN.i/VU.b.lp)^VU.c.lp)
    
    diag(VU.popmat[2:(VU.age.max+1),]) <- (VU.s.stoch[-(VU.age.max+1)])*VU.pred.red
    VU.popmat[VU.age.max+1,VU.age.max+1] <- 0
    VU.popmat[1,] <- VU.m.arr[i,,e]
    VU.n.mat[,i+1] <- VU.popmat %*% VU.n.mat[,i]

    VU.s.arr[i,,e] <- VU.s.stoch * VU.pred.red
    
  } # end i loop
  
  VU.n.sums.mat[e,] <- ((as.vector(colSums(VU.n.mat))/VU.pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

VU.n.md <- apply(VU.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
VU.n.up <- apply(VU.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
VU.n.lo <- apply(VU.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,VU.n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(VU.n.lo),1.05*max(VU.n.up)))
lines(yrs,VU.n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,VU.n.up,lty=2,col="red",lwd=1.5)

VU.s.add <- VU.m.add  <- rep(0, VU.age.max+1)
for (m in 1:iter) {
  VU.s.add <- rbind(VU.s.add, VU.s.arr[ceiling(VU.gen.l):(t+1),,m])
  VU.m.add <- rbind(VU.m.add, VU.m.arr[ceiling(VU.gen.l):(t+1),,m])
}
VU.s.add <- VU.s.add[-1,]
VU.m.add <- VU.m.add[-1,]

VU.s.md <- apply(VU.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
VU.s.up <- apply(VU.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
VU.s.lo <- apply(VU.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(VU.age.vec,VU.s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(VU.s.lo),1.05*max(VU.s.up)))
lines(VU.age.vec,VU.s.lo,lty=2,col="red",lwd=1.5)
lines(VU.age.vec,VU.s.up,lty=2,col="red",lwd=1.5)

VU.m.md <- apply(VU.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
VU.m.up <- apply(VU.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
VU.m.lo <- apply(VU.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(VU.age.vec,VU.m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(VU.m.lo),1.05*max(VU.m.up)))
lines(VU.age.vec,VU.m.lo,lty=2,col="red",lwd=1.5)
lines(VU.age.vec,VU.m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))

