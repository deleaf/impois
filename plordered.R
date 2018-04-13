# plausibility function for an IM based on an ordering of the Xs at each parameter value
# NOTE: this is the raw plausbility that does not adjust for connectivity/convexety issues
#       in the confidence belt.
# lambda = signal rates at which to compute plausibility
# x = the observed x value at which plausibility is calculated
# b = background rate
# r = ranking function for Xs
# NOTE: max(lambda)+b <= 200 and x <= 170 will be enforced to makes sure the G99 method works
raw.pl.ordered <- function(lambda, x, b, r) {
	if(max(lambda) + b > 200) { stop("max(lambda) + b > 200 not allowed") }
	if(x > 170) { stop("x > 170 not allowed") }
	pls = sapply(lambda, function(l, x, b, r) {
		X = 0:170
		R = r(X, l, b)

		# if the ordering puts max(X) before x, then we might need to extend the range of X and do the ordering again
		if(R[length(X)] >= R[which(X==x)]) { stop("ranking function goes to the end of X values before reaching x, might need to run with larger X range") }

		# not sure how to handle ties in the ranking below first place
		if(length(unique(R)) < length(X)) { 
			Rfreqs = table(R)
			Rdup.idx = as.numeric(Rfreqs) > 1
			if(min(as.numeric(names(Rfreqs)[Rdup.idx])) < max(R)) {
				stop("ranking function produced ties after first place") 
			}
		}

		if(R[which(X==x)] == max(R)) {
			belc = 0
		} else {		
			# get sort indices
			R.sort.idx = sort(R, decreasing=T, index.return=T)$ix
			x.sort.idx = which(R.sort.idx==which(X==x))
		
			# probability masses
			px.sort = dpois(X[R.sort.idx], l+b)
			belc = cumsum(px.sort)[x.sort.idx-1]
		}
		return(1-belc)
	}, x, b, r)
	return(data.frame(lambda, pl=pls))
}
