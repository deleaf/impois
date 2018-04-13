# plausibility functions implied by the constrained Poisson confidence intervals in the physics literature:
# Roe, B. P. and Woodroofe, M. (1999), "Improved probability method for estimating signal in the presence of background," Phys. Rev. D, 60, 053009.
# Mandelkern, M. and Schultz, J. (2000), "Coverage of confidence intervals based on conditional probability," J. High Energy Phys., 11, 036.
# Martin, R. (2017), "A mathematical characterization of confidence as valid belief," arXiv:1707.00486

library(plyr)

source("gen_belt.R")
source("plfill.R")

# conditional probability defined in RW`99 (eqn 4)
# NUMERATOR ONLY -- excludes ppois(n,b) denominator which cancels in LR calculation
rw.cond.prob.numer <- function(k,n,b,lambda)
{
	if(k <= n) {
		return( dpois(k,b+lambda) )
	} else {
		j = 0:n
		return( sum(dpois(j,b)*dpois(k-j,lambda)) )
	}
	
}

# conditional probability defined in RW`99 (eqn 4)
rw.cond.prob <- function(k,n,b,lambda)
{
	return(exp(log(rw.cond.prob.numer(k,n,b,lambda)) - log(ppois(n,b))))
}


# conditional likelihood ratio defined in RW`99
rw.cond.LR <- function(k,n,b,lambda)
{
	numer = rw.cond.prob.numer(k,n,b,lambda)
	
	delta = 2
	denom = optimize(function(l,k,n,b) {rw.cond.prob.numer(k,n,b,l)}, interval=c(max(0,k-b-delta),k+delta), maximum=T, k=k, n=n, b=b )$objective
	
	return(exp(log(numer) - log(denom)))
}

# returns minimum confidence level at which the RW99 acceptance region would include n
# this function's name is misleading
rw.cond.rank <- function(n, lambda, b) {
	K = 0:170
	R = sapply(K,rw.cond.LR, n, b, lambda)
		
	if(R[length(K)] >= R[which(K==n)]) { stop("RW99 ranking function goes to the end of K values before reaching n, might need to run with larger K range") }
		
	R.sort.idx = sort(R, decreasing=T, index.return=T)$ix
	n.sort.idx = which(R.sort.idx==which(K==n))	
	K.sorted = K[R.sort.idx]
	cond.pk.sorted = sapply(K.sorted, rw.cond.prob, n, b, lambda)
	
	if(n.sort.idx==1) {
		return(0)
	} else {
		return(cumsum(cond.pk.sorted)[n.sort.idx-1])
	}	
}

# get RW99 acceptance region for lambda
# gamma = confidence level
# return (left, right) bounds
rw.cond.ar <- function(lambda, b, gamma) {
	N = 0:50
	R = sapply(N, rw.cond.rank, lambda, b)
	
   # RW99 interval is min/max of N[R <= gamma]
   rw99.left=min(N[R <= gamma])
   rw99.right=max(N[R <= gamma])
   
   
	return(data.frame(lambda, rw99.left, rw99.right))
}
# Mandelkern and Schultz adjustment to the right bound of the belt to get at least gamma coverage probability
ms00.adj <- function(oldbelt, b, gamma) {
	newbelt = NULL
	right.old = -Inf
	for(i in 1:nrow(oldbelt)) {
		lambda = oldbelt[i,1]
		theta = lambda + b
		left = oldbelt[i,2]
		right = oldbelt[i,3]
		right.new = right
		cdf.left = ppois(left-1, theta)
		old.covg = ppois(right, theta) - cdf.left
		new.covg = ppois(right.new, theta) - cdf.left
		if(new.covg < gamma) {
			right.new = qpois(cdf.left + gamma, theta)
		} 
			
		# verify
		new.covg = ppois(right.new, theta) - cdf.left
		if(new.covg < gamma) {
			stop("MS00 coverage adjustment error")
		}
		 	
		newbelt = rbind(newbelt, c(lambda, left, right, right.new, old.covg, new.covg))
	}
	newbelt = as.data.frame(newbelt)

	names(newbelt) = c(names(oldbelt), "ms00.right", "old.covg", "new.covg")
	return(newbelt)
}

# find acceptance region for RW99 method with MS00 adjustment for fixed lambda and b
# values with confidence level 1-alpha
ar.rw99ms00 <- function(lambda, b, alpha) {
	conflev = 1-alpha
	rw99.ar = rw.cond.ar(lambda, b, conflev)
	rw99ms00.ar = ms00.adj(rw99.ar, b, conflev)
	return(list(left=rw99ms00.ar$rw99.left,right=rw99ms00.ar$ms00.right))
}

# compute plausibilities for a single lambda (signal rate) value with fixed b 
# (background rate) and oberved x using Roe and Woodroofe 99 method with Mandelkern 
# and Schultz 00 adjustment
# 	tol is a tolerance; the returned plausibilties will be within tol of the true plausibility
pl.rw99ms00.raw.lambda <- function(lambda, x, b, tol) {
	alpha.left = 0.0
	alpha.right = 1.0
	alpha.mid = (alpha.left + alpha.right)/2

	# binary search for largest alpha that does not include x in acceptance region
	while(alpha.right - alpha.left > tol) {
		ar = ar.rw99ms00(lambda, b, alpha.mid)
		found.mid = ar$left <= x & x <= ar$right
		
		if(found.mid) {
			alpha.left = alpha.mid
			alpha.mid = (alpha.left + alpha.right)/2
		} else {
			alpha.right = alpha.mid
			alpha.mid = (alpha.left + alpha.right)/2
		}
	}
	# return the largest alpha that does not include x in the acceptance region
	return(ifelse(found.mid, alpha.right, alpha.mid))
}

# compute plausibilities over a range of lambda (signal rate) values with fixed b 
# (background rate) and oberved x using Roe and Woodroofe 99 method with Mandelkern 
# and Schultz 00 adjustment
# 	tol is a tolerance; the returned plausibilties will be within tol of the true plausibility
#	  NP is the number of parallel processes to use; if NP>0 the doParallel library will used
pl.rw99ms00.raw <- function(lambda, x, b, tol, NP=0) {
	# set up parallel cluster if needed
	cl = NULL
	if(NP > 0) {
		if(require(doParallel)==0) { stop("doParallel package required for RW99MS00 method with NP>0") }
		cl <- makeCluster(NP)
		registerDoParallel(cl)
	}
	
	pls = laply(.data=lambda, .fun=pl.rw99ms00.raw.lambda, .parallel=!is.null(cl), .paropts=list(.export=ls(envir=globalenv())), x=x, b=b, tol=tol)
	
	# shut down parallel cluster if it was created
	if(!is.null(cl)) {stopCluster(cl)}
	
	return(data.frame(lambda,pl=pls))
}

# compute plausibilities over a range of lambda (signal rate) values with fixed b (background rate) and oberved x using Roe and Woodroofe 99 method with Mandelkern 
# and Schultz 00 adjustment
# tol is a tolerance; the returned plausibilties will be within tol of the true plausibility
#	NP is the number of parallel processes to be used by pl.rw99ms00.raw; if NP>0 the doParallel library will used
pl.rw99ms00 <- function(lambda, x, b, tol, NP=0) {
	return(adj.pl.fill(pl.rw99ms00.raw(lambda, x, b, tol, NP)))
}

plot.rw99ms00.rawadj <- function(x,b, tol) {
	lambda = seq(0,max(2*x,5),0.001)
	pl.rw99ms00.lambda = pl.rw99ms00(lambda,x,b,tol)$pl
	pl.rw99ms00.lambda.raw = pl.rw99ms00.raw(lambda,x,b,tol)$pl
	plot(lambda, pl.rw99ms00.lambda, type="l", xlim=c(min(lambda),max(lambda)), ylim=c(0,1), ylab="plausibility", xlab="signal rate (lambda)", main=paste("x=",x," b=",b,sep=""))
	lines(lambda, pl.rw99ms00.lambda.raw,lty=2)	
}

# Roe and Woodroofe (1999) 100*gamma % confidence interval for fixed b
rw99.int <- function(lambda.max, b, gamma) {
  lambda = seq(0,lambda.max,0.001)
  belt = gen.belt(lambda, b, gamma, rw.cond.ar)
  names(belt)[2:3] = c("left","right")
  return(invert.belt(belt))
}
