# invert confidence belt to create confidence interval
# belt is a data frame with columns lambda, left, right
# returns a data frame with columns x, lower, upper where 
#   lower is the smallest lambda with left <= x <= right in the input
#   and upper is the largest lambda with left <= x <= right
invert.belt <- function(belt) {
  max.x = max(belt$left)-1
  res = NULL
  for(x in 0:max.x) {
    lambda.x = belt[belt$left <= x & x <= belt$right,"lambda"]
    lower = min(lambda.x)
    upper = max(lambda.x)
    res = rbind(res, data.frame(x,lower,upper))
  }
  return(res)
}