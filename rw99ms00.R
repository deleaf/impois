# plausibility functions implied by the constrained Poisson confidence intervals in the physics literature:
# Roe, B. P. and Woodroofe, M. (1999), "Improved probability method for estimating signal in the presence of background," Phys. Rev. D, 60, 053009.
# Mandelkern, M. and Schultz, J. (2000), "Coverage of confidence intervals based on conditional probability," J. High Energy Phys., 11, 036.
# Martin, R. (2017), "A mathematical characterization of confidence as valid belief," arXiv:1707.00486

library(plyr)

source("gen_belt.R")
source("plordered.R")
source("plfill.R")
source("rw99core.R")

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
  if(x > 50) {stop("RW99MS00 plausibility not supported for x > 50")}
  
  # find the the position of x in the RW99 ordering
  N = 0:50
  R = sapply(N, rw.cond.rank, lambda, b)
  if(R[length(N)] <= R[which(N==x)]) { stop("RW99MS00 pl function goes to the end of N values before reaching x, might need to run with larger N range") }
  R.sort.idx = sort(R, index.return=T)$ix
  N.sorted = N[R.sort.idx]	
  
  # find the cumulative RW99 conditional probability over the ordering up to x
  cl.x = R[which(N==x)]
  
  # if x <= (left of) the value ranked first by RW99 then the MS00 adjustment doesn't matter
  if(cl.x == 0 | x <= N.sorted[1]) {
    # cl.x is the minimum confidence level at which x would be included in the acceptance region
    return(1-cl.x)
  }
  # if x is greater than (right of) the value ranked first by RW99 then we need to account for the MS00 
  # adjustment. Let S(lambda) be the largest RW99 acceptance region that excludes x.
  # Let P1 be the RW99 conditional probability and P2 be unconditional probability used by MS00 for 
  # computing coverage.  We have:
  #   P2(S(lambda)) + P2(MS00 adjustment tail) >= P1(S(lambda))
  # If P2(S(lambda)) >= P1(S(lambda)), then P2(MS00 adjustment tail)=0 (MS00 adjustment not needed)
  # If P2(S(lambda)) < P1(S(lambda)), then the MS00 adjustment tail will include x and we need to find a
  # smaller RW99 acceptance region such that the corresponding MS00 adjustment tail does not include x.  
  # The new RW99 acceptance region will have conditional probability less than P1(S(lambda)).  Therefore,
  # we only need to search alpha in [1-P1(S(lambda)),1] instead of alpha in [0,1].
  else {
  
  	alpha.left = 1-cl.x
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


r.rw99 <- function(X, lambda, b) {
  R = as.numeric(laply(.data=X, .fun=function(x, lambda, b) {
      # the functions in rw99core.R need to be in this function's scope in order for the laply to work
      source("rw99core.R")
      return(rw.cond.rank(x,lambda,b))
    }, lambda=lambda, b=b))
  return(1-R)
}

# raw plausibility of lambda range for fixed x and b
# determined only by ordering function, could have nonmonotone ("valleys" or "teeth")
pl.rw99.raw <- function(lambda, x, b, NP=0) {
  return(raw.pl.ordered(lambda, x, b, r.rw99, NP))
}

# monotone plausibility of lambda range for fixed x and b after filling any valleys
pl.rw99 <- function(lambda, x, b, NP=0) {
  pl.raw = raw.pl.ordered(lambda, x, b, r.rw99, NP)

  # set up parallel cluster if needed
  cl = NULL
  if(length(x) > 1 && NP > 0) {
    if(require(doParallel)==0) { stop("doParallel package required for NP>0") }
    cl <- makeCluster(min(length(x),NP))
    registerDoParallel(cl)
  }
   
  pl.adj = ddply(.data=pl.raw, .variables="x", .fun=adj.pl.fill, .parallel=!is.null(cl), .paropts=list(.export=ls(envir=globalenv())))

  # shut down parallel cluster if it was created
  if(!is.null(cl)) {stopCluster(cl)}
  
  return(pl.adj)
}

# plots the raw and filled plausibility in the same figure
plot.rw99.rawadj <- function(x,b,NP=0) {
  lambda = seq(0,max(2*x,5),0.001)
  pl.rw99.lambda = pl.rw99(lambda,x,b,NP)$pl
  pl.rw99.lambda.raw = pl.rw99.raw(lambda,x,b,NP)$pl
  plot(lambda, pl.rw99.lambda, type="l", xlim=c(min(lambda),max(lambda)), ylim=c(0,1), ylab="plausibility", xlab="signal rate (lambda)", main=paste("x=",x," b=",b,sep=""))
  lines(lambda, pl.rw99.lambda.raw,lty=2)	
}

