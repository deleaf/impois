# plausibility functions implied by the constrained Poisson confidence intervals in the physics literature:
# Feldman, G. J. and Cousins, R. D. (1998), "Unified approach to the classical statistical analysis of small signals," Phys. Rev. D, 57, 3873–3889.
# Martin, R. (2017), "A mathematical characterization of confidence as valid belief," arXiv:1707.00486

source("ar_ordered.R")
source("gen_belt.R")
source("invert_belt.R")
source("plordered.R")
source("plfill.R")

# determines ranking of X values for fixed lambda and b
r.fc98 <- function(X, lambda, b) {
	# the numerator becomes zero for large X and small lambda+b, so ranks are determined by logs instead
	#return(dpois(X,lambda+b) / dpois(X,pmax(X-b,rep(0,length(X)))+b))

	numer.theta = lambda+b
	denom.theta = pmax(X-b,rep(0,length(X)))+b
	res = -numer.theta+denom.theta + X*(log(numer.theta)-log(denom.theta))
	return(res)
}

# raw plausibility of lambda range for fixed x and b
# determined only by ordering function, could have nonmonotone ("valleys" or "teeth")
pl.fc98.raw <- function(lambda, x, b) {
	return(raw.pl.ordered(lambda, x, b, r.fc98))
}

# monotone plausibility of lambda range for fixed x and b after filling any valleys
pl.fc98 <- function(lambda, x, b) {
	return(adj.pl.fill(raw.pl.ordered(lambda, x, b, r.fc98)))
}

# plots the raw and filled plausibility in the same figure
plot.fc98.rawadj <- function(x,b) {
	lambda = seq(0,max(2*x,5),0.001)
	pl.fc98.lambda = pl.fc98(lambda,x,b)$pl
	pl.fc98.lambda.raw = pl.fc98.raw(lambda,x,b)$pl
	plot(lambda, pl.fc98.lambda, type="l", xlim=c(min(lambda),max(lambda)), ylim=c(0,1), ylab="plausibility", xlab="signal rate (lambda)", main=paste("x=",x," b=",b,sep=""))
	lines(lambda, pl.fc98.lambda.raw,lty=2)	
}

# Feldman and Cousins acceptance region
ar.fc98 <- function(lambda, b, gamma) {
  return(ar.ordered(lambda, b, gamma, r.fc98))
}

# Feldman and Cousins 100*gamma % confidence interval
fc98.int <- function(lambda.max, b, gamma) {
  lambda = seq(0,lambda.max,0.001)
  belt = gen.belt(lambda, b, gamma, ar.fc98)
  return(invert.belt(belt))
}




