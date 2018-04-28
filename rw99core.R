# functions to implement the RW99 method

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
# this function's name is misleading because lower confidence is higher rank
rw.cond.rank <- function(n, lambda, b) {
  K = 0:(max(qpois(0.9999999,lambda+b), n)+10)
  R = sapply(K,rw.cond.LR, n, b, lambda)
  
  R.sort.idx = sort(R, decreasing=T, index.return=T)$ix
  n.sort.idx = which(R.sort.idx==which(K==n))	
  K.sorted = K[R.sort.idx]
  cond.pk.sorted = sapply(K.sorted, rw.cond.prob, n, b, lambda)
  cumpk.sorted = cumsum(cond.pk.sorted)

  if(R[length(K)] >= R[which(K==n)]) {
    if(cumpk.sorted[which(K.sorted == K[length(K)])] >= 0.9999999) {
      return(1)
    } else {
      stop(paste("RW99 ranking function goes to the end of K values (",max(K),") before reaching n=",n," with (lambda,b)=(",lambda,',',b,"); might need to run with larger K range",sep=''))
    }
  } else if(n.sort.idx==1) {
    return(0)
  } else {
    return(cumpk.sorted[n.sort.idx-1])
  }	
}

