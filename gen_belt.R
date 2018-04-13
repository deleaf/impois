# create level gamma confidence belt from over a range of lambda values for fixed b
# ar is the function that returns the acceptance region bounds for each lambda
gen.belt <- function(lambda, b, gamma, ar) {
  belt = ldply(lambda, ar, .progress="text", b=b, gamma=gamma)
  return(belt)
}

