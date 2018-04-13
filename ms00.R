# Mandelkern, M. and Schultz, J. (2000), "The statistical analysis of Gaussian and Poisson signals near physical boundaries," J. Math. Phys., 41, 5701-5709.

# compute plausibilities over a range of lambda (signal rate) values with fixed b (background rate) and oberved x using Mandelkern and Schultz 00 method
# tol is a tolerance; the returned plausibilties will be within tol of the true plausibility
pl.ms00 <- function(lambda, x, b, tol) {
	pls = sapply(lambda, function(l, x, b, tol) {
		alpha.left = 0.0
		alpha.right = 1.0
		alpha.mid = (alpha.left + alpha.right)/2

		while(alpha.right - alpha.left > tol) {
			#print(c(alpha.left,alpha.right))
			ar = ar.ms00(l, b, alpha.mid)
			found.mid = ar$left <= x & x <= ar$right

			if(found.mid) {
				alpha.left = alpha.mid
				alpha.mid = (alpha.left + alpha.right)/2
			} else {
				alpha.right = alpha.mid
				alpha.mid = (alpha.left + alpha.right)/2
			}
		}
		# return the largest alpha that does not include X in the acceptance region
		return(ifelse(found.mid, alpha.right, alpha.mid))
	}, x, b, tol)
	
	return(data.frame(lambda,pl=pls))
}


# find acceptance region in X-space for fixed lambda and b values and confidence level 1-alpha
ar.ms00 <- function(lambda, b, alpha) {
	conflev = 1-alpha
	right = qpois((1+conflev)/2, lambda+b)
	while(ppois(right, lambda+b) < (1+conflev)/2){
		right = right + 1
	}
	if(right < b) { right = floor(b)+1 }
		
	left = qpois(ppois(right, lambda+b)-conflev, lambda+b)+1
	
	#if interval has coverage < conflev, shift left bound
	while((left >= 0) && (ppois(right, lambda+b) - ppois(left-1, lambda+b) < conflev)) {
		left = left - 1
	}

	#if left bound < b, then left bound := 0 and adjust right bound for coverage
	if(left <= floor(b)) { 
		left = 0
			
		right = qpois(conflev, lambda+b)
		#while( ppois(right, lambda+b) < cl ) {
		#	right = right + 1
		#}
	}
	return(list(left=left,right=right))
}


### probability mass function for MS00 MLE
ms00lik <- function(lambda.hat, lambda, b) {
	ifelse(lambda.hat==0, ppois(floor(b), lambda+b), dpois(lambda.hat+b, lambda+b))
}



