# comparison of constrained Poisson plausibility functions derived from the physics
# literature and the elastic belief score balance (EBSB) IM

# POINT THIS TO YOUR WORKING DIRECTORY THAT HAS THE R PROGRAMS FROM THE REPOSITORY
setwd("~/impois_update")

#####################################################################
################ REPLICATION OF PUBLISHED INTERVALS #################
#####################################################################

source("fc98.R")
# Table IV in Feldman and Cousins (1998), b=3.0 column
# Table 2 of Mandelkern (2002)
round(fc98.int(26, 3, 0.9),2)
# FC98 90% acceptance region at lambda=1.08
ar.fc98(1.08,3,0.9)

source("g99.R")
# Fig. 3 of Giunti (1999)
round(g99.int(15, 2, 0.9),2)
round(g99.int(15, 3, 0.9),2)
round(g99.int(15, 4, 0.9),2)
round(g99.int(15, 5, 0.9),2)

source("rw99ms00.R")
# Table 1 in Roe and Woodroofe (1999)
# Table I in Mandelkern and Schultz (2000)
round(rw99.int(19, 3, 0.9), 2)

# clean up
rm(list=ls())

#####################################################################
###### SHOW WHY MS00 METHOD CAN'T BE CONVERTED TO PL. FUNCTION ######
#####################################################################

source("ms00.R")
library(plyr)

# plot MS00 likelihood at different lambda values with b=3
# overlay 90% acceptance region to demonstrate why MS00 did not encounter this problem
b = 3
lambda = seq(1.9,3.6,length.out=15)
pdf(file="ms00bimodal.pdf",width=9,height=7)
par(mfrow=c(3,5))
for(l in lambda) {
	ms00.ar90 = ar.ms00(l, b, 0.1)
	#print(ar)
	lambda.hat = c(0,((floor(b)+1):12)-b)
	lik = sapply(lambda.hat, ms00lik, lambda=l, b)
	plot(lambda.hat, lik, ylim=c(0,0.35), main=bquote(lambda==.(round(l,2))), xlab=expression(hat(lambda)), ylab="Probability")
	abline(v=max(0,ms00.ar90$left-b), lty=2)
	abline(v=max(0,ms00.ar90$right-b), lty=2)
}
## MS00 likelihood is bimodal, so not clear how to construct low confidence intervals
dev.off()

# clean up
rm(list=ls())

#####################################################################
####### SHOW WHY RW99MS00 DID NOT REPLICATE PUBLISHED RESULTS #######
#####################################################################

# bounds from Table 1 of 
# Mandelkern, M. and Schultz, J. (2000), "Coverage of confidence intervals based on conditional probability," J. High Energy Phys., 11, 036.
rw99ms00.int = data.frame(
	x = 0:12,
	lower = c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.16, 0.90,  1.66,  2.44,  2.98,  3.75,  4.52),
	upper = c(2.44, 2.95, 3.75, 4.80, 6.01, 7.28, 8.42, 9.58, 11.02, 12.23, 13.51, 14.77, 16.01)
)

# invert bounds to get coverage
lambda = seq(0,4.5,0.0001)
rw99ms00.belt = ldply(lambda, function(l, rw99ms00.int) {
	cov.belt = rw99ms00.int[rw99ms00.int$lower <= l & l <= rw99ms00.int$upper,]
	x.left = min(cov.belt$x)
	x.right = max(cov.belt$x)
	return(data.frame(lambda=l, left=x.left, right=x.right))
}, rw99ms00.int=rw99ms00.int)

# compute coverage probability
rw99ms00.covg = ddply(rw99ms00.belt, .variables="lambda", function(df) {
	theta = df$lambda + 3
	return(data.frame(coverage=ppois(df$right, theta) - ppois(df$left-1, theta)))
})

# plot coverage by lambda value
pdf(file="rw99ms00undercoverage.pdf",width=8,height=7)
plot(rw99ms00.belt$lambda, rw99ms00.covg$coverage, type="l", xlab=expression(lambda), ylab="Coverage probability")
abline(h=0.9, col="grey")
dev.off()

# clean up
rm(list=ls())

#####################################################################
################# SHOW HOW RW99MS00 MIGHT BE NESTED #################
#####################################################################

lambda = 4.5
b = 3

# 41% RW99 acceptance region
rw41 = rw.cond.ar(lambda, b, 0.41)
ppois(rw41$rw99.right, lambda+b) - ppois(rw41$rw99.left-1, lambda+b)
# coverage probability exceeds 0.41 => MS00 adjustment is not needed
# 41% RW99MS00 acceptance region
rwms41 = list(left = rw41$rw99.left, right=rw41$rw99.right)

# 40% RW99 acceptance region
rw40 = rw.cond.ar(lambda, b, 0.40)
ppois(rw40$rw99.right, lambda+b) - ppois(rw40$rw99.left-1, lambda+b)
# coverage probability is less than 0.40 => MS00 adjustment is needed
# 40% RW99MS00 acceptance region
rwms40 = ar.rw99ms00(lambda, b, 1-0.4)

# the 40% acceptance region is not nested within the 41% acceptance region
rwms40
rwms41
nested = rwms41$left <= rwms40$left & rwms40$right <= rwms41$right
nested

# clean up
rm(list=ls())

#####################################################################
############## MONOTONICITY ISSUE IN ORDERING METHODS  ##############
#####################################################################

# Use Giunti (1999) ordering at lambda=1.158 and 1.441 with b=3 for example
source("g99.R")
source("fc98.R")
b=3
X = 0:20
outtab = NULL
for(lambda in c(1, 1.158, 1.441)) {
  R = r.g99(X, lambda, b)
  R.sort.idx = sort(R, decreasing=T, index.return=T)$ix

  # ordered X's
  X.sorted = X[R.sort.idx]

  # cumulative probability over ordered X's
  cumprob = cumsum(dpois(X.sorted, lambda+b))

  # plausibility of lambda at each ordered X
  pl = c(1,1-cumprob[1:(length(cumprob)-1)])  
  
  res = round(data.frame(rho=1:length(X.sorted), x=X.sorted, cumprob, pl)[1:7,],2)
  print(paste("lambda =",lambda,"b =",b))
  print(res)
  if(is.null(outtab)) { outtab = res}
  else { outtab = cbind(outtab, res) }
}
outtab


#####################################################################
################### COMPARE EBSB, FC98, G99, RW99 ###################
#####################################################################

library(plyr)
library(ggplot2)

source("fc98.R")
source("g99.R")
source("rw99ms00.R")
# the following file can be downloaded from https://www4.stat.ncsu.edu/~rmartin/Codes/impois.R
source("impois.R")

# wrapper function for EB-SB method: compute raw plausibility then fill
pl.ebsb <- function(lambda, x, b, NP) {
  # set up parallel cluster if needed
  cl = NULL
  if(length(x) > 1 && NP > 0) {
    if(require(doParallel)==0) { stop("doParallel package required for NP>0") }
    cl <- makeCluster(min(length(x),NP))
    registerDoParallel(cl)
  }
  
  pl.raw = ldply(x, function(x,lambda,b) {
      # the pl.pois function must be within scope of this function for ldply to work
      source("impois.R")
      return(cbind(data.frame(x),pl.pois(lambda+b,x,PRS="Score-Balanced",b)))
    }, .parallel=!is.null(cl), .paropts=list(.export=ls(envir=globalenv())), lambda=lambda, b=b) 
  pl.raw$lambda = pl.raw$theta-b
  pl.raw$b = b
  
  pl.adj = ddply(.data=pl.raw, .variables="x", .fun=adj.pl.fill, .parallel=!is.null(cl), .paropts=list(.export=ls(envir=globalenv())))
  
  # shut down parallel cluster if it was created
  if(!is.null(cl)) {stopCluster(cl)}
  
  return(pl.adj)
}

# generate plausibiliies for fixed x and b across methods
#	pltol is the tolerance for methods that approximate the plausibility
#		the approximate plausibility will be within pltol of the true plausilibility
#	NP is the number of parallel processes to use
#	verbose==TRUE prints progress updates
# returns a list containing the input x and b values, the lambda (signal rate) values
#	at which plausibility was computed, and the plausibilities
#	pl.ebsb = EBSB method
#	pl.fc98 = Feldman and Cousins (1998)
#	pl.g99 = Giunti (1999)
#	pl.rw99ms00 = Roe and Woodroofe (1999) with Mandelkern and Schultz (2000) coverage adjustment
gen.pl.phys <- function(x, b, NP=0, verbose=TRUE) {
	lambda = seq(0,max(2.3*max(x),5),0.001)

  pl.rw99.lambda = pl.rw99(lambda,x,b, NP)
  if(verbose) {print("RW99MS00 done")}

  pl.ebsb.lambda = pl.ebsb(lambda, x, b, NP)
  if(verbose) {print("EB-SB done")}
  
  pl.fc98.lambda = pl.fc98(lambda,x,b, NP)
  if(verbose) {print("FC98 done")}
  
  pl.g99.lambda = pl.g99(lambda,x,b, NP)
  if(verbose) {print("G99 done")}
  
  return(list(x=x, b=b, lambda=lambda, 
		pl.ebsb=pl.ebsb.lambda,
		pl.fc98=pl.fc98.lambda,
		pl.g99=pl.g99.lambda,
		pl.rw99 = pl.rw99.lambda
	))
}


# plot plausibility functions from the output of gen.pl.phys()
plot.pl.phys <- function(pl.phys) {
	plot(pl.phys$lambda, pl.phys$pl.ebsb$pl, type="l", xlim=c(min(pl.phys$lambda),max(pl.phys$lambda)), ylim=c(0,1), ylab="Plausibility", xlab=expression(lambda), main=paste("x=",pl.phys$x," b=",pl.phys$b,sep=""))
	lines(pl.phys$lambda, pl.phys$pl.fc98$pl,lty=2)
	lines(pl.phys$lambda, pl.phys$pl.g99$pl,lty=3)
	lines(pl.phys$lambda, pl.phys$pl.rw99$pl,lty=4)
	legend(legend=c("EB-SB","FC98","G99","RW99"),lty=1:4,"topright")
}

# maximum number of parallel processes to use
NP = 8
# print status messages while processing
verbose = TRUE

b=3
x=7
pl.phys.7.3 = gen.pl.phys(x,b,NP,verbose)
pdf(file="pl_phys_7_3.pdf",width=8,height=7)
plot.pl.phys(pl.phys.7.3)

x=2
pl.phys.2.3 = gen.pl.phys(x,b,NP,verbose)
pdf(file="pl_phys_2_3.pdf",width=8,height=7)
plot.pl.phys(pl.phys.2.3)

bpoints = c(0.5,1,3)
xpoints = c(0,1,2,3,4,5,7,10)

pl.phys.all = ldply(.data=bpoints, .fun =
	function(b,x,NP,verbose) {
		print(paste("b =",b))
		pl.phys = gen.pl.phys(x,b,NP,verbose)
		
		#return(data.frame(lambda=pl.phys$lambda, EBSB=pl.phys$pl.ebsb, 
		#	FC98=pl.phys$pl.fc98, G99=pl.phys$pl.g99, RW99MS00=pl.phys$pl.rw99ms00))
		
		# create stacked data frame for plotting with ggplot
		df.ebsb = cbind(data.frame(Method=rep("EB-SB",nrow(pl.phys$pl.ebsb)), pl.phys$pl.ebsb))
		print(names(df.ebsb))
		df.fc98 = cbind(data.frame(Method=rep("FC98",nrow(pl.phys$pl.fc98)), pl.phys$pl.fc98))
		print(names(df.fc98))
		df.g99  = cbind(data.frame(Method=rep("G99",nrow(pl.phys$pl.g99)), pl.phys$pl.g99))
		print(names(df.g99))
	  df.rw99 = cbind(data.frame(Method=rep("RW99",nrow(pl.phys$pl.rw99)), pl.phys$pl.rw99))
	  print(names(df.rw99))
		
		return(rbind(df.fc98, df.g99, df.rw99, df.ebsb[,names(df.ebsb) != "theta"]))	
	}
	, .progress="text", x=xpoints, NP=NP, verbose=verbose)
	
save(pl.phys.all,file="pl.phys.all.Rdata")

# turn method name into factor to control the legend order in the plot
pl.phys.all$Method = factor(pl.phys.all$Method, levels=c("G99","RW99","FC98","EB-SB"))

# plot color for each method
colors = c("G99" = "lightgrey", "RW99" = "grey", "FC98" = "darkgrey", "EB-SB" = "black")

ggplot(pl.phys.all, aes(lambda, pl, color=Method)) + geom_line() + facet_grid(x ~ b) + xlab(expression(lambda)) + ylab("Plausibility") + ylim(0,1) + xlim(0,15) + theme(legend.position="bottom") + theme_classic() + scale_colour_manual(values=colors)
ggsave("plphys.pdf", width=8, height=7, units="in")


