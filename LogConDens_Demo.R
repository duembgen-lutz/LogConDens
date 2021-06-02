source('LogConDens.R')


### Example 1: Gaussian data ###

# Sample size:
n <- 1000

# Sample from standard Gaussian distribution:
x <- sort(rnorm(n))

## Estimation of the log-concave density

# First a version with graphical displays of all
# intermediate steps:
# (If you are not using RStudio, please make sure that
#  the graphics window is visible beside the
#  console window, and the latter should be in the foreground):

res <- LogConDens0(x)

# The same without displaying intermediate steps:

res <- LogConDens(x)
system.time(res <- LogConDens(x))

# Analyzing the result:

# Plot of directional derivatives for additional
# knots (should all be <= 0 with equality at the
#  actual estimated knots):

plot(res$x,res$DLtau,type='l',
	xlab=expression(italic(tau)),
	ylab=expression(italic(DL(phi,V[tau]))))
rug(x)
abline(h=0,col='red')
points(res$knots,rep(0,length(res$knots)),col='blue')

# Plot of estimated (black) and true (green) log-density:

xx <- LocalInterpolate.1A(res$x,7)
phixx <- LocalInterpolate.1A(res$phix,7)

# True density and log-density:
f0xx <- dnorm(xx)
phi0xx <- log(f0xx)

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


### Example 2: Exponential data ###

# Sample size:
n <- 1000

# Sample from tandard exponential distribution:
x <- sort(-log(runif(n)))

## Estimation of the log-concave density

# With graphical displays of all intermediate steps:

res <- LogConDens0(x)

# The same without displaying intermediate steps:

res2 <- LogConDens(x)
system.time(res2 <- LogConDens(x))

# Analyzing the result:

# Plot of directional derivatives for additional
# knots (should all be <= 0 with equality at the
#  actual estimated knots):

plot(res$x,res$DLtau,type='l',
	xlab=expression(italic(tau)),
	ylab=expression(italic(DL(phi,V[tau]))))
rug(x)
abline(h=0,col='red')
points(res$knots,rep(0,length(res$knots)),col='blue')

# Plot of estimated (black) and true (green) log-density:

xx <- LocalInterpolate.1A(res$x,7)
phixx <- LocalInterpolate.1A(res$phix,7)

# True density and log-density:
f0xx <- exp(-xx)
phi0xx <- -xx

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


### Example 2: Exponential data ###

# Sample size:
n <- 1000

# Sample from tandard exponential distribution:
x <- sort(-log(runif(n)))

### Estimation of the log-concave density ###

# With graphical displays of all intermediate steps:

res <- LogConDens0(x)

# The same without displaying intermediate steps:

res2 <- LogConDens(x)
system.time(res2 <- LogConDens(x))

# Analyzing the result:

# Plot of directional derivatives for additional
# knots (should all be <= 0 with equality at the
#  actual estimated knots):

plot(res$x,res$DLtau,type='l',
	xlab=expression(italic(tau)),
	ylab=expression(italic(DL(phi,V[tau]))))
rug(x)
abline(h=0,col='red')
points(res$knots,rep(0,length(res$knots)),col='blue')

# Plot of estimated (black) and true (green) log-density:

xx <- LocalInterpolate.1A(res$x,7)
phixx <- LocalInterpolate.1A(res$phix,7)

# True density and log-density:
f0xx <- exp(-xx)
phi0xx <- -xx

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


### Example 3: Logistic data ###

# Sample size:
n <- 1000

# Sample from standard logistic distribution:
u <- sort(runif(n))
x <- log(u/(1-u))

## Estimation of the log-concave density

# First a version with graphical displays of all
# intermediate steps:
# (If you are not using RStudio, please make sure that
#  the graphics window is visible beside the
#  console window, and the latter should be in the foreground):

res <- LogConDens0(x)

# The same without displaying intermediate steps:

res <- LogConDens(x)
system.time(res <- LogConDens(x))

# Analyzing the result:

# Plot of directional derivatives for additional
# knots (should all be <= 0 with equality at the
#  actual estimated knots):

plot(res$x,res$DLtau,type='l',
	xlab=expression(italic(tau)),
	ylab=expression(italic(DL(phi,V[tau]))))
rug(x)
abline(h=0,col='red')
points(res$knots,rep(0,length(res$knots)),col='blue')

# Plot of estimated (black) and true (green) log-density:

xx <- LocalInterpolate.1A(res$x,7)
phixx <- LocalInterpolate.1A(res$phix,7)

# True density and log-density:
f0xx <- 1/(exp(xx) + exp(-xx) + 2)
phi0xx <- log(f0xx)

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


### Example 4: Gamma data ###

# Sample size:
n <- 1000

# Sample from gamma distribution:
shape <- 3
x <- sort(rgamma(n,shape))

## Estimation of the log-concave density

# First a version with graphical displays of all
# intermediate steps:
# (If you are not using RStudio, please make sure that
#  the graphics window is visible beside the
#  console window, and the latter should be in the foreground):

res <- LogConDens0(x)

# The same without displaying intermediate steps:

res <- LogConDens(x)
system.time(res <- LogConDens(x))

# Analyzing the result:

# Plot of directional derivatives for additional
# knots (should all be <= 0 with equality at the
#  actual estimated knots):

plot(res$x,res$DLtau,type='l',
	xlab=expression(italic(tau)),
	ylab=expression(italic(DL(phi,V[tau]))))
rug(x)
abline(h=0,col='red')
points(res$knots,rep(0,length(res$knots)),col='blue')

# Plot of estimated (black) and true (green) log-density:

xx <- LocalInterpolate.1A(res$x,7)
phixx <- LocalInterpolate.1A(res$phix,7)

# True density and log-density:
f0xx <- dgamma(xx,shape)
phi0xx <- log(f0xx)

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


### Example 5: Beta data ###

# Sample size:
n <- 1000

# Sample from gamma distribution:
shape1 <- 2
shape2 <- 3
x <- sort(rbeta(n,shape1,shape2))

## Estimation of the log-concave density

# First a version with graphical displays of all
# intermediate steps:
# (If you are not using RStudio, please make sure that
#  the graphics window is visible beside the
#  console window, and the latter should be in the foreground):

res <- LogConDens0(x)

# The same without displaying intermediate steps:

res <- LogConDens(x)
system.time(res <- LogConDens(x))

# Analyzing the result:

# Plot of directional derivatives for additional
# knots (should all be <= 0 with equality at the
#  actual estimated knots):

plot(res$x,res$DLtau,type='l',
	xlab=expression(italic(tau)),
	ylab=expression(italic(DL(phi,V[tau]))))
rug(x)
abline(h=0,col='red')
points(res$knots,rep(0,length(res$knots)),col='blue')

# Plot of estimated (black) and true (green) log-density:

xx <- LocalInterpolate.1A(res$x,7)
phixx <- LocalInterpolate.1A(res$phix,7)

# True density and log-density:
f0xx <- dbeta(xx,shape1,shape2)
phi0xx <- log(f0xx)

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
