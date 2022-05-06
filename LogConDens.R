#############################################################
### Active Set Algorithm for the Computation of the NPMLE ###
### of a Log-Concave Probability Density on the Real Line ###
#############################################################
#
# Lutz Duembgen, Alexandre Moesching, Christof Straehl
# December 2020
#
# Source:
# L. Duembgen, A. Moesching, C. Straehl (2021).
# Active Set Algorithms for Estimating Shape-Constrained
# Density Ratios.
# Computational Statistics and Data Analysis 163.
# arXiv:1808.09340

### Main programs ----

# Input:
# X: Vector of observations.
# W: Vector of corresponding weights (with default 1).
# m0: Number of initially deactivated constraints plus 1
#     (= number of knots minus 1).
# delta0: Parameter for the algorithms' thresholds:
#    For local searches (i.e. Newton procedure with
#     shape constraints), the threshold is
#         delta_Newton = delta0 / sample size.
#     For checking global optimality, the threshold is
#         delta_Knot = emp. standard deviation * delta_Newton.
# gamma0: Parameter for regularisation of Hessian matrix
#         in the Newton steps.

LogConDens0 <- function(X, W=rep(1,length(X)),
	m0=4, delta0=10^(-7), gamma0=10^(-5))
# Computation of log-concave density estimator with
# initial Gaussian fit on a grid of m0+1 knots.
# Graphical display of all intermediate steps...
{
	# Preparations of raw data:
	n <- length(X)
	if (any(X[1:(n-1)] >= X[2:n])){
		x <- X
		w <- W/sum(W)
	}else{
		tmp <- LocalPrepareData(X,W)
		x <- tmp$x
		w <- tmp$w
		n <- tmp$n
	}
	# Initialize knots and start with
	# Gaussian MLE phi:
	d <- floor((n-1)/m0)
	knotsind <- 1 + d*(0:m0)
	knotsind[m0+1] <- n
	knots <- x[knotsind]
	ww <- LocalLinearSplines2.1A(knots,x,w)
	mu <- sum(w*x)
	sigma <- sqrt(sum(w*(x - mu)^2))
	phi <- log(dnorm(knots,mean=mu,sd=sigma))
	phi <- LocalNormalize.1A(knots,phi)
	LL <- sum(ww*phi)
	message <- readline("First plot:")
	plot(knots,phi,type='l',lwd=2,
		xlab='x',ylab='phi(x)',
		main=paste('Starting point: ',
			'LL = ',round(LL,digits=9),sep=''))
	points(knots,phi,pch=16)
	rug(x)
	
	# Bookkeeping:
	nr.local.search <- 0
	nr.Newton <- 0
	
	# Prepare for first local search:
	maxDLtau <- Inf
	isnewknot <- rep(FALSE,length(knots))
	delta1 <- sigma*delta0/n
	
	while (maxDLtau > delta1)
	{
		# Store current value of LL:
		LL.old <- LL
		
		# Start local search:
		nr.local.search <- nr.local.search + 1
		phi.old <- phi
		knots.old <- knots
		proposal <- LocalNewton.1A(knots,ww,phi,gamma0=gamma0)
		nr.Newton <- nr.Newton + 1
		nr.knots <- matrix(c(length(knots),sum(isnewknot)),2,1)
		dimnames(nr.knots)[[1]] <- c('total','new')
		iter <- 0
		while (proposal$dirderiv > delta0/n && iter < 100)
		{
			phi.new <- proposal$phi.new
			message <- readline("New proposal:")
			plot(knots.old,phi.old,type='l',col='cyan',
				xlab='x',ylab='phi(x), phi.new(x)',
				ylim=range(c(phi.old,phi,phi.new)),
				main=paste('LL = ',round(LL,digits=9),sep=''))
			abline(v=knots[isnewknot],col='green')
			lines(knots,phi,lwd=2)
			rug(x)
			lines(knots,phi.new,col='blue')
			points(knots,phi.new,col='blue')
			
			tmp <- LocalStepForward.1A(knots,phi,phi.new,ww,isnewknot)
			phi <- tmp$phi
			if (length(phi) < length(knots))
			{
				LL.ref <- -Inf
				# Update of knots, ww and isnewknot:
				knots     <- tmp$knots
				ww        <- tmp$w
				isnewknot <- tmp$isnewknot
				# Graphics:
				message <- readline("2nd step size correction:")
				lines(knots,phi,col='magenta')
				points(knots,phi,col='magenta')
			}else{
				LL.ref <- LL
			}
			# Update of LL:
			LL <- sum(ww*phi)
			if (LL > LL.ref)
			{
				proposal <- LocalNewton.1A(knots,ww,phi,
					gamma0=gamma0)
				nr.Newton <- nr.Newton + 1
				nr.knots <- cbind(nr.knots,
					c(length(knots),sum(isnewknot)))
			}else{
				proposal$dirderiv <- 0
			}
			iter <- iter+1
		}
		message <- readline("Local optimum:")
		plot(knots.old,phi.old,type='l',col='cyan',
			xlab='x',ylab='phi(x)',ylim=range(phi),
			main=paste('Local optimum: ',
				'LL = ',round(LL,digits=9),sep=''))
		lines(knots,phi,lwd=2)
		points(knots,phi,pch=16)
		rug(x)
		print(paste('Local search ',nr.local.search,
			', numbers of knots:',sep=''),quote=FALSE)
		print(nr.knots)
		# Check for global optimality:
		if (LL > LL.old){
			DLtau <- LocalOptimality.1A(x,w,knots,phi)
			maxDLtau <- max(DLtau)
			if (maxDLtau > delta1){
				# Add new knots:
				delta2 <- max(delta1,maxDLtau*0.0001)
				newpars <- LocalNewKnots(x,knots,phi,DLtau,delta2)
				knots <- newpars$knots
				phi <- newpars$phi
				isnewknot <- newpars$isnewknot
				# Graphics:
				message <- readline("New knots:")
				plot(knots,phi,type='l',lwd=2,
					xlab='x',ylab='phi(x)',ylim=range(phi),
					main=paste('LL = ',round(LL,digits=9),sep=''))
				abline(v=knots[isnewknot],col='green')
				points(knots[isnewknot],phi[isnewknot],col='green')
				rug(x)
				# Update of ww:
				ww <- LocalLinearSplines2.1A(knots,x,w)
			}
		}else{
			maxDLtau <- 0
		}
	}
	
	message <- readline("Final estimate:")
	plot(knots,phi,type='l',lwd=2,
		xlab=expression(italic(x)),
		ylab=expression(italic(phi(x))),
		ylim=range(phi),
		main=paste('Global optimum: ',
			'LL = ',round(LL,digits=9),sep=''))
	points(knots,phi,pch=16)
	rug(x)
	
	print(paste('Number of local searches:',nr.local.search),
		quote=FALSE)	
	print(paste('Number of Newton steps:',nr.Newton),
		quote=FALSE)
	DD <- LocalLinearSplines1.1A(knots,x)
	phix <- DD %*% phi
	
	return(list(knots=knots,phi=phi,x=x,w=w,
		phix=phix,LL=LL,DLtau=DLtau,
		nr.local.search=nr.local.search,
		nr.Newton=nr.Newton))
}


LogConDens1 <- function(X, W=rep(1,length(X)),
	m0=4, delta0=10^(-7), gamma0=10^(-5))
# Same as LogConDens0, without graphics but with
# counting of local searches and Newton steps.
{
	# Preparations of raw data:
	n <- length(X)
	if (sum(X[1:(n-1)] >= X[2:n]) == 0){
		x <- X
		w <- W/sum(W)
	}else{
		tmp <- LocalPrepareData(X,W)
		x <- tmp$x
		w <- tmp$w
		n <- tmp$n
	}
	# Initialize knots and start with
	# Gaussian MLE phi:
	d <- floor((n-1)/m0)
	knotsind <- 1 + d*(0:m0)
	knotsind[m0+1] <- n
	knots <- x[knotsind]
	ww <- LocalLinearSplines2.1A(knots,x,w)
	mu <- sum(w*x)
	sigma <- sqrt(sum(w*(x - mu)^2))
	phi <- log(dnorm(knots,mean=mu,sd=sigma))
	phi <- LocalNormalize.1A(knots,phi)
	LL <- sum(ww*phi)
	
	# Bookkeeping:
	nr.local.search <- 0
	nr.Newton <- 0
	
	# Prepare for first local search:
	maxDLtau <- Inf
	isnewknot <- rep(FALSE,length(knots))
	delta1 <- sigma*delta0/n
	
	while (maxDLtau > delta1)
	{
		# Store current value of LL:
		LL.old <- LL
		
		# New local search:
		nr.local.search <- nr.local.search + 1
		proposal <- LocalNewton.1A(knots,ww,phi,gamma0=gamma0)
		nr.Newton <- nr.Newton + 1
		iter <- 0
		while (proposal$dirderiv > delta0/n && iter < 100)
		{
			phi.new <- proposal$phi.new
			tmp <- LocalStepForward.1A(knots,phi,phi.new,
				ww,isnewknot)
			phi <- tmp$phi
			if (length(phi) < length(knots)){
				LL.ref <- -Inf
				# Update of knots, ww and isnewknot:
				knots     <- tmp$knots
				ww        <- tmp$w
				isnewknot <- tmp$isnewknot
			}else{
				LL.ref <- LL
			}
			# Update of LL:
			LL <- sum(ww*phi)
			if (LL > LL.ref){
				proposal <- LocalNewton.1A(knots,ww,phi,
					gamma0=gamma0)
				nr.Newton <- nr.Newton + 1
			}else{
				proposal$dirderiv <- 0
			}
			iter <- iter+1
		}
		if (LL > LL.old){
			# Check for global optimality:
			DLtau <- LocalOptimality.1A(x,w,knots,phi)
			maxDLtau <- max(DLtau)
			if (maxDLtau > delta1){
				# Add new knots:
				delta2 <- max(delta1,maxDLtau*0.0001)
				newpars <- LocalNewKnots(x,knots,phi,DLtau,delta2)
				knots <- newpars$knots
				phi <- newpars$phi
				isnewknot <- newpars$isnewknot
				ww <- LocalLinearSplines2.1A(knots,x,w)
			}
		}else{
			# Stop the whole algorithm:
			maxDLtau <- 0
		}
	}
	
	DD <- LocalLinearSplines1.1A(knots,x)
	phix <- DD %*% phi
	
	return(list(knots=knots,phi=phi,x=x,w=w,
		phix=phix,LL=LL,DLtau=DLtau,
		nr.local.search=nr.local.search,
		nr.Newton=nr.Newton))
}


LogConDens <- function(X, W=rep(1,length(X)),
	m0=4, delta0=10^(-7), gamma0=10^(-5))
# Computation of log-concave density estimator with
# initial Gaussian fit on a grid of m0+1 knots
# (m0-1 deactivated constraints).
# Pure computation, no graphical displays, no counting...
{
	# Preparations of raw data:
	n <- length(X)
	if (sum(X[1:(n-1)] >= X[2:n]) == 0){
		x <- X
		w <- W/sum(W)
	}else{
		tmp <- LocalPrepareData(X,W)
		x <- tmp$x
		w <- tmp$w
		n <- tmp$n
	}
	# Initialize knots and start with
	# Gaussian MLE phi:
	d <- floor((n-1)/m0)
	knotsind <- 1 + d*(0:m0)
	knotsind[m0+1] <- n
	knots <- x[knotsind]
	ww <- LocalLinearSplines2.1A(knots,x,w)
	mu <- sum(w*x)
	sigma <- sqrt(sum(w*(x - mu)^2))
	phi <- log(dnorm(knots,mean=mu,sd=sigma))
	phi <- LocalNormalize.1A(knots,phi)
	LL <- sum(ww*phi)
	
	# Prepare for first local search:
	maxDLtau <- Inf
	isnewknot <- rep(FALSE,length(knots))
	delta1 <- sigma*delta0/n
	
	while (maxDLtau > delta1){
		# Store current value of LL:
		LL.old <- LL
		
		# New local search:
		proposal <- LocalNewton.1A(knots,ww,phi,gamma0=gamma0)
		iter <- 0
		while (proposal$dirderiv > delta0/n && iter < 100){
			phi.new <- proposal$phi.new
			tmp <- LocalStepForward.1A(knots,phi,phi.new,ww,
				isnewknot)
			phi <- tmp$phi
			if (length(phi) < length(knots)){
				LL.ref <- -Inf
				# Update of knots, ww and isnewknot:
				knots     <- tmp$knots
				ww        <- tmp$w
				isnewknot <- tmp$isnewknot
			}else{
				LL.ref <- LL
			}
			# Update of LL:
			LL <- sum(ww*phi)
			if (LL > LL.ref){
				proposal <- LocalNewton.1A(knots,ww,phi,
					gamma0=gamma0)
			}else{
				proposal$dirderiv <- 0
			}
			iter <- iter+1
		}
		if (LL > LL.old){
			# Check for global optimality:
			DLtau <- LocalOptimality.1A(x,w,knots,phi)
			maxDLtau <- max(DLtau)
			if (maxDLtau > delta1){
				# Add new knots:
				delta2 <- max(delta1,maxDLtau*0.0001)
				newpars <- LocalNewKnots(x,knots,phi,DLtau,delta2)
				knots <- newpars$knots
				phi <- newpars$phi
				isnewknot <- newpars$isnewknot
				ww <- LocalLinearSplines2.1A(knots,x,w)
			}
		}else{
			# Stop whole algorithm:
			maxDLtau <- 0
		}
	}
	
	DD <- LocalLinearSplines1.1A(knots,x)
	phix <- DD %*% phi
	
	return(list(knots=knots,phi=phi,x=x,w=w,
		phix=phix,LL=LL,DLtau=DLtau))
}


### Auxiliary programs 4 ----
### Checking for new potential knots 

LocalNewKnots <- function(x,knots,phi,DLtau,delta=0)
# Deactivation of constraints: Between any two consecutive
# knots knots[j] and knots[j+1], the constraint at x[i] is
# deactivated, if the local maximum of DLtau[.] is attained
# at i and greater than delta.
{
	m <- length(knots)
	knotsnew <- rep(NA,2*m-1)
	knotsnew[2*(1:m) - 1] <- knots
	phinew <- rep(NA,2*m-1)
	phinew[2*(1:m) - 1] <- phi
	isnewknot <- rep(FALSE,2*m-1)
	for (j in 1:(m-1)){
		II <- which(x > knots[j] & x < knots[j+1])
		if (length(II) > 0 && max(DLtau[II]) > delta){
			k <- which.max(DLtau[II])
			i <- II[k]
			j2 <- 2*j
			knotsnew[j2] <- x[i]
			lambda <- c(knots[j+1] - x[i],x[i] - knots[j])/
				(knots[j+1] - knots[j])
			phinew[j2] <- lambda[1]*phi[j] + lambda[2]*phi[j+1]
			isnewknot[j2] <- TRUE
		}
	}
	tmp <- !is.na(knotsnew)
	knotsnew <- knotsnew[tmp]
	phinew <- phinew[tmp]
	isnewknot <- isnewknot[tmp]
	return(list(knots=knotsnew,phi=phinew,isnewknot=isnewknot))
}

LocalOptimality.1A <- function(x,w,knots,phi)
# Determine the directional derivatives
# DL(phi,V_{tau,phi}), tau in {x[1],x[2],...}
{
	DLtau <- rep(0,length(x))
	for (j in 1:(length(knots)-1)){
		itmp <- which(x > knots[j] & x < knots[j+1])
		mtmp <- length(itmp)
		if (mtmp == 1){
			dtmp <- knots[j+1] - knots[j]
			x01tmp <- x[itmp] - knots[j]
			x10tmp <- knots[j+1] - x[itmp]
			phitmp <- (x10tmp*phi[j] + x01tmp*phi[j+1])/dtmp
			wtmp <- w[itmp]
			DLtau[itmp] <- wtmp*x10tmp*x01tmp/dtmp
			j01tmp <- x01tmp*J10.1A(phitmp,phi[j])
			j10tmp <- x10tmp*J10.1A(phitmp,phi[j+1])
			DLtau[itmp] <- DLtau[itmp] -
				x01tmp*x10tmp*(j01tmp + j10tmp)/dtmp
		}
		if (mtmp > 1){
			dtmp <- knots[j+1] - knots[j]
			x01tmp <- x[itmp] - knots[j]
			x10tmp <- knots[j+1] - x[itmp]
			phitmp <- (x10tmp*phi[j] + x01tmp*phi[j+1])/dtmp
			wtmp <- w[itmp]
			DLtau[itmp] <- (x10tmp*cumsum(wtmp*x01tmp) +
				x01tmp*cumsum(c(0,(wtmp*x10tmp)[mtmp:2]))[mtmp:1])/dtmp
			j01tmp <- x01tmp*J10.1A(phitmp,rep(phi[j],mtmp))
			j10tmp <- x10tmp*J10.1A(phitmp,rep(phi[j+1],mtmp))
			DLtau[itmp] <- DLtau[itmp] -
				x01tmp*x10tmp*(j01tmp + j10tmp)/dtmp
		}
	}
	return(DLtau)
}


### Auxiliary programs 3 ----
### 2nd step size correction

LocalConcavity.1A <- function(knots,phi)
# Determine negative changes of slope,
# all of which should be > 0 if phi defines
# a strictly concave function on {knots[i] : ...}.
{
	m <- length(knots)
	if (m == 2){
		return(c(Inf,Inf))
	}
	dk <- knots[2:m] - knots[1:(m-1)]
	dphi <- phi[2:m] - phi[1:(m-1)]
	dphidk <- dphi/dk
	chofslope <- dphidk[1:(m-2)] - dphidk[2:(m-1)]
	return(c(Inf,chofslope,Inf))
}

LocalStepForward.1A <- function(knots,phi,phi.new,w,
	isnewknot=rep(FALSE,length(knots)))
# Replace phi with
#    (1 - t0)*phi + t0*phi.new ,
# where t0 is the largest number in [0,1] such that
# the new phi defines still a concave function.
# Then normalize phi additively to become a log-probability
# density.
# Return new tuple of true knots, the new phi
# and the updated weight vector w.
# In case of an input argument isnewknot, the procedure
# checks first whether phi.new is strictly concave at all
# new knots.
{
	m <- length(knots)
	# 1st easy case:
	if (m == 2){
		phi <- LocalNormalize.1A(knots,phi.new)
		return(list(phi=phi))
	}
	# Changes of slope:
	chofsl.new <- LocalConcavity.1A(knots,phi.new)
	# 2nd easy (and ideal) case:
	if (min(chofsl.new) > 0){
		phi <- LocalNormalize.1A(knots,phi.new)
		return(list(phi=phi))
	}
	# All other cases with min(chofsl.new) <= 0 lead to
	# removal of one knot / activation of one constraint:
	if (any(isnewknot) && min(chofsl.new[isnewknot]) <= 0){
		# At least one of the new knots has to be removed:
		tmpi <- which(chofsl.new <= 0 & isnewknot)
		j0 <- which.min(chofsl.new[tmpi])
		j0 <- tmpi[j0]
		phi <- phi[-j0]
		tmp <- (knots[j0:(j0+1)] - knots[(j0-1):j0])/
			(knots[j0+1] - knots[j0-1])
		w[j0-1] <- w[j0-1] + tmp[2]*w[j0]
		w[j0+1] <- w[j0+1] + tmp[1]*w[j0]
		knots <- knots[-j0]
		w <- w[-j0]
		isnewknot <- isnewknot[-j0]
		return(list(knots=knots,phi=phi,w=w,isnewknot=isnewknot))
	}
	# At each knot, phi or phi.new is strictly concave,
	# and we replace phi.new with a proper convex combination
	# of phi and phi.new such that phi.new is concave:
	tmpi  <- which(chofsl.new <= 0)
	chofsl <- LocalConcavity.1A(knots,phi)
	t0 <- chofsl[tmpi]/(chofsl[tmpi] - chofsl.new[tmpi])
	j0 <- which.min(t0)
	t0 <- t0[j0]
	j0 <- tmpi[j0]
	phi <- (1 - t0)*phi[-j0] + t0*phi.new[-j0]
	tmp <- (knots[j0:(j0+1)] - knots[(j0-1):j0])/
		(knots[j0+1] - knots[j0-1])
	w[j0-1] <- w[j0-1] + tmp[2]*w[j0]
	w[j0+1] <- w[j0+1] + tmp[1]*w[j0]
	knots <- knots[-j0]
	w <- w[-j0]
	isnewknot <- isnewknot[-j0]
	# Now make sure that phi is even strictly concave
	# at each knot:
	chofsl <- LocalConcavity.1A(knots,phi)
	while (min(chofsl) <= 0){
		j0  <- which.min(chofsl)
		phi <- phi[-j0]
		tmp <- (knots[j0:(j0+1)] - knots[(j0-1):j0])/
			(knots[j0+1] - knots[j0-1])
		w[j0-1] <- w[j0-1] + tmp[2]*w[j0]
		w[j0+1] <- w[j0+1] + tmp[1]*w[j0]
		knots <- knots[-j0]
		w <- w[-j0]
		isnewknot <- isnewknot[-j0]
	}
	phi <- LocalNormalize.1A(knots,phi)
	return(list(knots=knots,phi=phi,w=w,isnewknot=isnewknot))
}


### Auxiliary programs 2 ----
### Normalization, log-likelihood plus derivatives
### and Newton step with 1st step size correction

LocalNormalize.1A <- function(x,phi)
# Normalize phi such that it defines
# a log-probability density
{
	n <- length(x)
	dx <- x[2:n] - x[1:(n-1)]
	tmp00 <- J00.1A(phi[1:(n-1)],phi[2:n])
	integral <- sum(tmp00*dx)
	phi.new <- phi - log(integral)
	return(phi.new)
}

LocalLL.1A <- function(x,w,phi)
# Log-likelihood
{
	n <- length(x)
	dx <- x[2:n] - x[1:(n-1)]
	LL <- sum(w*phi)
	tmp00 <- J00.1A(phi[1:(n-1)],phi[2:n])
	LL <- LL + 1 - sum(tmp00*dx)
	return(LL)
}

LocalLL2.1A <- function(x,w,phi)
# Log-likelihood plus its gradient vector
# and minus Hessian matrix.
# Only for checking numerical accuracy...
{
	n <- length(x)
	dx <- x[2:n] - x[1:(n-1)]
	LL <- sum(w*phi)
	tmp00 <- J00.1A(phi[1:(n-1)],phi[2:n])
	LL <- LL + 1 - sum(tmp00*dx)
	GLL <- w
	tmp10 <- J10.1A(phi[1:(n-1)],phi[2:n])
	GLL[1:(n-1)] <- GLL[1:(n-1)] - dx*tmp10
	tmp01 <- J10.1A(phi[2:n],phi[1:(n-1)])
	GLL[2:n] <- GLL[2:n] - dx*tmp01
	MHLL <- matrix(0,n,n)
	tmp20 <- J20.1A(phi[1:(n-1)],phi[2:n])
	A20 <- cbind(1:(n-1),1:(n-1))
	MHLL[A20] <- dx*tmp20
	tmp02 <- J20.1A(phi[2:n],phi[1:(n-1)])
	A02 <- cbind(2:n,2:n)
	MHLL[A02] <- MHLL[A02] + dx*tmp02
	tmp11 <- J11.1A(phi[1:(n-1)],phi[2:n])
	A11a <- cbind(1:(n-1),2:n)
	A11b <- cbind(2:n,1:(n-1))
	MHLL[A11a] <- dx*tmp11
	MHLL[A11b] <- dx*tmp11
	return(list(LL=LL,GLL=GLL,MHLL=MHLL))
}

LocalNewton.1A <- function(x,w,phi,
	delta0=10^(-9),gamma0=10^(-5))
# Newton step with step size correction
{
	n <- length(x)
	dx <- x[2:n] - x[1:(n-1)]
	LL <- sum(w*phi)
	tmp00 <- J00.1A(phi[1:(n-1)],phi[2:n])
	LL <- LL + 1 - sum(tmp00*dx)
	GLL <- w
	tmp10 <- J10.1A(phi[1:(n-1)],phi[2:n])
	GLL[1:(n-1)] <- GLL[1:(n-1)] - dx*tmp10
	tmp01 <- J10.1A(phi[2:n],phi[1:(n-1)])
	GLL[2:n] <- GLL[2:n] - dx*tmp01
	MHLL <- matrix(0,n,n)
	tmp20 <- J20.1A(phi[1:(n-1)],phi[2:n])
	A20 <- cbind(1:(n-1),1:(n-1))
	MHLL[A20] <- dx*tmp20
	tmp02 <- J20.1A(phi[2:n],phi[1:(n-1)])
	A02 <- cbind(2:n,2:n)
	MHLL[A02] <- MHLL[A02] + dx*tmp02
	tmp11 <- J11.1A(phi[1:(n-1)],phi[2:n])
	A11a <- cbind(1:(n-1),2:n)
	A11b <- cbind(2:n,1:(n-1))
	MHLL[A11a] <- dx*tmp11
	MHLL[A11b] <- dx*tmp11
	gamma <- gamma0*mean(abs(diag(MHLL)))
	diag(MHLL) <- diag(MHLL) + gamma
	dphi <- qr.solve(MHLL,GLL)
	dirderiv <- sum(dphi * GLL)
	delta <- dirderiv
	phi.new <- phi + dphi
	LL.new <- LocalLL.1A(x,w,phi.new)
	while (LL.new < LL + delta/3 && delta > delta0/n){
		phi.new <- (phi + phi.new)/2
		delta <- delta/2
		LL.new <- LocalLL.1A(x,w,phi.new)
	}
	if (LL.new > LL){
		return(list(phi.new=phi.new,dirderiv=dirderiv))
	}else{
		return(list(phi.new=phi,dirderiv=0))
	}	
}


### Auxiliary programs 1 ----
### Linear splines and
### basic function J(.,.) plus 1st and 2nd
### partial derivatives thereof.

LocalLinearSplines1.1A <- function(knots,x=knots)
# Computes a design matrix Dx for the B-spline basis
# of linear splines with given knots and evaluated
# at the components of the vector x.
# We assume that knots has strictly increasing
# components
#   knots[1] < knots[2] < ... < knots[m] ,
# and that the basis functions f_1 and f_m are
# extrapolated linearly on (-Inf,knots[2]] and
# [knots[m-1],Inf), respectively.
{
	m <- length(knots)
	dk <- knots[2:m] - knots[1:(m-1)]
	n <- length(x)
	Dx <- matrix(0,n,m)
	for (j in 1:(m-1)){
		tmp <- (x > knots[j] & x <= knots[j+1])
		Dx[tmp,j] <- (knots[j+1] - x[tmp])/dk[j]
		Dx[tmp,j+1] <- (x[tmp] - knots[j])/dk[j]
	}
	if (min(x) <= knots[1]){
		tmp <- (x <= knots[1])
		Dx[tmp,1] <- (knots[2] - x[tmp])/dk[1]
	}
	if (max(x) > knots[m]){
		tmp <- (x > knots[m])
		Dx[tmp,m] <- (x[tmp] - knots[m-1])/dk[m-1]
	}
	return(Dx)
}

LocalLinearSplines2.1A <- function(knots,x=knots,w=rep(1,length(x)))
# For linear splines with given knots
#   knots[1] < knots[2] < ... < knots[m]
# and arbitrary vectors x and w of equal length with w >= 0,
# this procedure returns a weight vector wk, such that for
# the B-spline basis functions f_1, f_2, ..., f_m,
#    wk[j] = sum_{i=1}^n w[i]*f_j(x[i]) .
# Here we assume that f_1 and f_m are extrapolated linearly
# to the intervals (-Inf,knots[2]] and [knots[m-1],Inf),
# respectively.
{
	Dx <- LocalLinearSplines1.1A(knots,x)
	wk <- colSums(Dx * w)
	return(wk)
}

J00.1A <- function(x,y,delta=0.01)
# J00.1A(x,y) = integral_0^1 exp((1-u)x + uy) du.
{
	J00xy <- exp((x + y)/2)
	tmp1 <- (abs(x - y) > delta)
	z1 <- (y[tmp1] - x[tmp1])/2
	J00xy[tmp1] <- J00xy[tmp1]*(sinh(z1)/z1)
	tmp0 <- !tmp1
	J00xy[tmp0] <- J00xy[tmp0]*
		exp((x[tmp0] - y[tmp0])^2/24)
	return(J00xy)
}

J10.1A <- function(x,y,delta=0.01)
# J10.1A(x,y) = integral_0^1 (1-u)*exp((1-u)x + uy) du.
{
	J10xy <- rep(0,length(x))
	tmp1 <- (abs(x - y) > delta)
	x1 <- x[tmp1]
	y1 <- y[tmp1]
	z1 <- (y1 - x1)/2
	J10xy[tmp1] <- exp((x1 + y1)/2)*
		(sinh(z1) - z1*exp(-z1))/z1^2/2
	tmp0 <- !tmp1
	x0 <- x[tmp0]
	y0 <- y[tmp0]
	J10xy[tmp0] <- exp((2*x0 + y0)/3 +
		(x0 - y0)^2/36 - (x0 - y0)^3/810)/2
	return(J10xy)
}

J20.1A <- function(x,y,delta=0.01)
# J20.1A(x,y) = integral_0^1 (1-u)^2*exp((1-u)x + uy) du.
{
	J20xy <- rep(0,length(x))
	tmp1 <- (abs(x - y) > delta)
	x1 <- x[tmp1]
	y1 <- y[tmp1]
	z1 <- (y1 - x1)/2
	J20xy[tmp1] <- exp((x1 + y1)/2)*
		(sinh(z1)/z1 - (1+z1)*exp(-z1))/z1^2/2
	tmp0 <- !tmp1
	x0 <- x[tmp0]
	y0 <- y[tmp0]
	J20xy[tmp0] <- exp((3*x0 + y0)/4 +
		3*(x0 - y0)^2/160 - (x0 - y0)^3/960)/3
	return(J20xy)
}

J11.1A <- function(x,y,delta=0.01)
# J11.1A(x,y) = integral_0^1 (1-u)*u*exp((1-u)x + uy) du.
{
	J11xy <- exp((x + y)/2)/2
	tmp1 <- (abs(x - y) > delta)
	z1 <- (y[tmp1] - x[tmp1])/2
	J11xy[tmp1] <- J11xy[tmp1]*
		(cosh(z1) - sinh(z1)/z1)/z1^2
	tmp0 <- !tmp1
	J11xy[tmp0] <- J11xy[tmp0]*
		exp((x[tmp0] - y[tmp0])^2/40)/3
	return(J11xy)
}

LocalInterpolate.1A <- function(x,d)
# For a vector x with n = length(x) >= 2, this function returns
# a vector xx with n*d - d + 1 components, such that
# xx[(j*d-d+1):(j*d)] == x[j] + (x[j+1] - x[j])*(1:d)/d
# for 1 <= j < n.
{
	n <- length(x)
	xx <- rep(NA,n*d-d+1)
	xx[1] <- x[1]
	for (k in 1:d){
		lambda <- k/d
		xx[1 + k + (0:(n-2))*d] <-
			(1 - lambda)*x[1:(n-1)] + lambda*x[2:n]
	}
	return(xx)
}


### Auxiliary program 0:

LocalPrepareData <- function(X,W)
# Replace a data vector X with weights W
# by a vector x with strictly increasing
# components such that {x[.]} = {X[.]}
# and corresponding weights w.
{
	n <- length(X)
	tmp <- order(X)
	X <- X[tmp]
	W <- W[tmp]
	XX <- c(X,Inf)
	WW <- cumsum(W)
	tmp <- which(XX[1:n] < XX[2:(n+1)])
	x <- X[tmp]
	n <- length(x)
	WW <- WW[tmp]
	w <- WW - c(0,WW[1:(n-1)])
	w <- w/WW[n]
	return(list(x=x,w=w,n=n))
}