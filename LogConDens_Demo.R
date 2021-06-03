source('LogConDens.R')

### Simulation of data ###

# Sample size:
n <- 1000

# # Standard Gaussian distribution:
# x <- sort(rnorm(n))

# # Standard exponential distribution:
# x <- sort(-log(runif(n)))

# # Standard logistic distribution:
# u <- sort(runif(n))
# x <- log(u/(1-u))

# Gamma distribution:
shape <- 3
x <- sort(rgamma(n,shape))

# # Beta distribution:
# shape1 <- 2
# shape2 <- 3
# x <- sort(rbeta(n,shape1,shape2))

### Estimation of the log-concave density ###

# First a version with graphical displays of all
# intermediate steps:
# (If you are not using RStudio, please make sure that
#  the graphics window is visible beside the
#  console window, and the latter should be in the foreground):

res <- LogConDens0_old(x)
res <- LogConDens0(x)


# The same without displaying intermediate steps:

system.time(res1 <- LogConDens1_old(x))
c('Loc.search'=res1$nr.local.search,
  'Newton.steps'=res1$nr.Newton)
system.time(res2 <- LogConDens1(x))
c('Loc.search'=res2$nr.local.search,
  'Newton.steps'=res2$nr.Newton)

# The same without displaying intermediate steps
# and without any bookkeeping:

t1 <- system.time(res1 <- LogConDens_old(x))
t2 <- system.time(res2 <- LogConDens(x))
rbind(t1, t2, t1/t2)

# Analyzing the result:

# Plot of directional derivatives for additional
# knots (should all be <= 0 with equality at the
#  actual estimated knots):
# par(mai=c(0.8,0.9,0.1,0.1))
plot(res$x,res$DLtau,type='l',
	xlab=expression(italic(tau)),
	ylab=expression(italic(DL(phi,V[tau]))))
rug(x)
abline(h=0,col='red')
points(res$knots,rep(0,length(res$knots)),col='blue')

# Plot of estimated (black) and true (green) log-density:
xx <- LocalInterpolate.1A(res$x,7)
phixx <- LocalInterpolate.1A(res$phix,7)

# # Standard Gaussian truth:
# f0xx <- dnorm(xx)
# phi0xx <- log(f0xx)

# # Standard exponential truth:
# f0xx <- exp(-xx)
# phi0xx <- -xx

# # Standard logistic truth:
# f0xx <- 1/(exp(xx) + exp(-xx) + 2)
# phi0xx <- log(f0xx)

# Gamma truth:
f0xx <- dgamma(xx,shape)
phi0xx <- log(f0xx)

# # Beta truth:
# f0xx <- dbeta(xx,shape1,shape2)
# phi0xx <- log(f0xx)

plot(xx,phi0xx,type='l',lwd=3,col='green',
	xlab=expression(italic(x)),
	ylab=expression(italic(phi(x))),
	ylim=range(c(res$phi,phi0xx)))
lines(xx,phi0xx,col='darkgreen')
lines(res$x,res$phix,lwd=2)
points(res$knots,res$phi,pch=16)
rug(x)

# Plot of estimated (black) and true (green) density:
plot(xx,f0xx,type='l',lwd=3,col='green',
	xlab=expression(italic(x)),
	ylab=expression(italic(f(x))),
	ylim=c(0,max(c(exp(res$phi),f0xx))))
lines(xx,f0xx,col='darkgreen')
lines(xx,exp(phixx),lwd=2)
rug(x)


####################################
### Comparisons of running times ###
####################################

n <- 1000

x <- sort(rnorm(n))
# x <- sort(rbeta(n,shape1,shape2))
# u <- runif(n,-1,1)
# x <- sign(u)*log(abs(u))

t1 <- system.time(res1 <- LogConDens_old(x,m0=1))
t2 <- system.time(res2 <- LogConDens(x,m0=1))
t1
t2
t1/t2


# Impact of m0:
m0v <- c(1,2,4,8,16,32,64)
mcsim <- 1000
n <- 100000
RT <- matrix(Inf,mcsim,length(m0v))
dimnames(RT)[[2]] <- paste('m0=',m0v,sep='')
for (s in 1:mcsim){
  print(paste("Simulation",s,"..."),quote=FALSE)
  x <- sort(rnorm(n))
#  x <- sort(- log(runif(n)))
  for (j in 1:length(m0v)){
    RT[s,j] <- system.time(res2 <- LogConDens(x,m0=m0v[j]))[3]
  }
}
head(RT)
par(cex=1,mai=c(0.5,0.5,0.1,0.1))
boxplot(RT,lty=1,lwd=1.5,ylim=c(0,max(RT)))

# Gaussian distribution:
write.table(RT,quote=FALSE,sep="\t",
            file=paste("RT_m0_n",n,".txt",sep=''))
# Exponential distribution:
write.table(RT,quote=FALSE,sep="\t",
            file=paste("RT_m0_n",n,"_exp.txt",sep=''))


### Simulations for paper:

n <- 100000
mcsim <- 200

RT <- matrix(n,mcsim,4)
dimnames(RT)[[2]] <- c('n','old','new','old/new')

s <- 1

while (s <= mcsim){
  print(paste('Simulation',s,':'))
  x <- sort(rnorm(n))
  RT[s,2] <- system.time(res <- LogConDens_old(x,m0=4))[3]
  RT[s,3] <- system.time(res <- LogConDens(x,m0=4))[3]
  RT[s,4] <- RT[s,2]/RT[s,3]
  s <- s+1
}

head(RT)
boxplot(RT[,2:3],lty=1,lwd=1.5)
boxplot(RT[,4],lty=1,lwd=1.5)

write.table(RT,file="RT_ratios.txt",
            row.names=FALSE,sep="\t",
            quote=FALSE,append=TRUE)

# Having stored all results:

ds <- read.table("RT_ratios.txt",header=TRUE,sep="\t")
str(ds)

# Ratios of running times:

par(cex=1,mai=c(0.5,0.5,0.1,0.1))
boxplot(ds[,4] ~ ds[,1],ylim=c(0,max(RT[,4])),lty=1,lwd=1.5,
        names=paste("n =",c('100','200','500','10^3','10^4','10^5')))
abline(h=c(0,1),col='blue',lty=2)

lm(ds[,4] ~ factor(ds[,1]) - 1)$coefficients

# Running times of new algorithm:

par(cex=1,mai=c(0.5,0.5,0.1,0.1))
boxplot(log10(ds[,3]) ~ ds[,1],lty=1,lwd=1.5,
        names=paste("n =",c('100','200','500','10^3','10^4','10^5')))
abline(h=c(-1,0,1),col='blue',lty=2)

lm(ds[,3] ~ factor(ds[,1]) - 1)$coefficients
