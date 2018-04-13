# adjust the raw plausibilities to fill in valleys
# this has the same effect as connecting the plausibility/confidence interval
# at each x when the ordering function does not produce a convex/connected 
# interval
# pldf = output from raw.pl.ordered
adj.pl.fill <- function(pldf) {
	for(gamma in sort(pldf$pl)) {
		t1 = min(pldf[pldf$pl >= gamma,"lambda"])
		t2 = max(pldf[pldf$pl >= gamma,"lambda"])
		pldf[t1 <= pldf$lambda & pldf$lambda <= t2 & pldf$pl < gamma,"pl"] = gamma
	}
	return(pldf)
}
