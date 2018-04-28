library(plyr)

raw.pl.ordered.lambda <- function(l, xs, b, r) {
  X = 0:(max(qpois(0.9999999,l+b), max(xs))+10)

  R = r(X, l, b)
  # get sort indices
  R.sort.idx = sort(R, decreasing=T, index.return=T)$ix
  # probability masses in ranked order of X
  px.sort = dpois(X[R.sort.idx], l+b)
  # cumulative probability over ranked X
  cumpx.sort = cumsum(px.sort)
  
  pl = NULL
  for(x in xs) {
    # ranking function evaluated at x
    R.x = r(x,l,b)
    
    if(R.x == max(R)) {
      # if x has highest rank
      belc = 0
    } else if(R.x <= min(R)) {
      # if x is ranked lower than every value in X then x the plausibility is at least 0.9999999
      #   by definition of max(X)
      belc = 1
    } else {
      # if x is somewhere between first and last in the ranking 

      # not sure how to handle ties in the ranking below first place
      if(length(unique(R)) < length(X)) { 
        Rfreqs = table(R)
        Rdup.idx = as.numeric(Rfreqs) > 1
        Rdups = as.numeric(names(Rfreqs)[Rdup.idx])
        # ties for first place are allowed
        # if there are ties that come after x in the ranking, we can ignore the problem
        R.max = max(R)
        if(min(Rdups) < R.max  && R.x <= max(Rdups[Rdups < R.max])) {
          stop("ranking function produced ties between first place and x") 
        }
      }

      x.sort.idx = which(R.sort.idx==which(X==x))
      belc = cumpx.sort[x.sort.idx-1]
    }
    if( length(1-belc) == 0 || length(l) == 0 || length(b) == 0 || length(x) == 0) {
      stop(paste("zero length error: lambda=",l," b=",b," x=",x," pl=",1-belc,sep=''))
    }
    
    pl = rbind(pl,data.frame(lambda=l, b, x, pl=1-belc))
  }
  return(pl)
}


# plausibility function for an IM based on an ordering of the Xs at each parameter value
# NOTE: this is the raw plausbility that does not adjust for connectivity/convexety issues
#       in the confidence belt.
# lambda = signal rates at which to compute plausibility
# xs = the x values at which plausibility is calculated
# b = background rate
# r = ranking function for Xs
# NOTE: max(lambda)+b <= 200 and x <= 170 will be enforced to makes sure the G99 method works
raw.pl.ordered <- function(lambda, xs, b, r, NP=0) {
	if(max(lambda) + b > 200) { stop("max(lambda) + b > 200 not allowed") }
	if(max(xs) > 170) { stop("x > 170 not allowed") }
  
  # set up parallel cluster if needed
  cl = NULL
  if(NP > 0) {
    if(require(doParallel)==0) { stop("doParallel package required for NP>0") }
    cl <- makeCluster(NP)
    registerDoParallel(cl)
  }
  
	pls = ldply(.data=lambda, .fun=raw.pl.ordered.lambda, .parallel=!is.null(cl), .paropts=list(.export=ls(envir=globalenv())), xs=xs, b=b, r=r)

	# shut down parallel cluster if it was created
	if(!is.null(cl)) {stopCluster(cl)}
	
	return(pls)
}

