# plausibility functions implied by the constrained Poisson confidence intervals in the physics literature:
# Giunti, C. (1999), "New ordering principle for the classical statistical analysis of Poisson processes with background," Phys. Rev. D, 59, 053001.
# Martin, R. (2017), "A mathematical characterization of confidence as valid belief," arXiv:1707.00486

source("ar_ordered.R")
source("gen_belt.R")
source("invert_belt.R")
source("plordered.R")
source("plfill.R")

# determines ranking of X values for fixed lambda and b
# this function will error with max(X) > 170 because the factorial becomes too large
r.g99 <- function(X, lambda, b) {

  est.new = sapply(X, function(x,lambda,b) {
    K = 0:x
    K.factorials = factorial(K)
	  if(max(K.factorials)==Inf) { 
		  first.overflow = K[min(which(K.factorials==Inf))]
		  stop(paste("G99 rho: infinite factorial for k=",first.overflow," with x=",x,sep='')) 
	  }
	
	  # compute reference mean as ratio of cumulative sums
	  temp = exp(K*log(b) - log(K.factorials))
	  est.new.numer = sum(K*temp)
	  est.new.denom = sum(temp)
	
	  est.new = x + 1 - exp(log(est.new.numer) - log(est.new.denom))
	  return(est.new)
  }, lambda, b)
  
	# the numerator becomes zero for large X and small lambda+b, so ranks are determined by logs instead
	#return(dpois(X,lambda+b) / dpois(X,est.new+b))

	numer.theta = lambda+b
	denom.theta = est.new+b
	res = -numer.theta+denom.theta + X*(log(numer.theta)-log(denom.theta))
	return(res)
}

# raw plausibility of lambda range for fixed x and b
# determined only by ordering function, could have nonmonotone ("valleys" or "teeth")
pl.g99.raw <- function(lambda, x, b, NP=0) {
	return(raw.pl.ordered(lambda, x, b, r.g99, NP))
}

# monotone plausibility of lambda range for fixed x and b after filling any valleys
pl.g99 <- function(lambda, x, b, NP=0) {
  pl.raw = raw.pl.ordered(lambda, x, b, r.g99, NP)
  
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
plot.g99.rawadj <- function(x,b,NP=0) {
	lambda = seq(0,max(2*x,5),0.001)
	pl.g99.lambda = pl.g99(lambda,x,b,NP)$pl
	pl.g99.lambda.raw = pl.g99.raw(lambda,x,b,NP)$pl
	plot(lambda, pl.g99.lambda, type="l", xlim=c(min(lambda),max(lambda)), ylim=c(0,1), ylab="plausibility", xlab="signal rate (lambda)", main=paste("x=",x," b=",b,sep=""))
	lines(lambda, pl.g99.lambda.raw,lty=2)	
}

# Giunti (1999) acceptance region for fixed lambda and b
ar.g99 <- function(lambda, b, gamma) {
  return(ar.ordered(lambda, b, gamma, r.g99))
}

# Giunti (1999) 100*gamma % confidence interval
g99.int <- function(lambda.max, b, gamma) {
  lambda = seq(0,lambda.max,0.001)
  belt = gen.belt(lambda, b, gamma, ar.g99)
  return(invert.belt(belt))
}

