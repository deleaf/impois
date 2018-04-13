# level gamma acceptance region from ordering of X's at fixed lambda and b
# r is the ranking function
ar.ordered <- function(lambda, b, gamma, r) {
  X = 0:170
  R = r(X, lambda, b)
  # get sort indices
  R.sort.idx = sort(R, decreasing=T, index.return=T)$ix
 
  # probability masses
  px.sort = dpois(X[R.sort.idx], lambda+b)
  cumprob = cumsum(px.sort)
  # first index with cumulative probability >= gamma
  last.x.idx = min(which(cumprob >= gamma))
  left = min(X[R.sort.idx[1:last.x.idx]])
  right= max(X[R.sort.idx[1:last.x.idx]])
  return(data.frame(lambda,left,right))
}

