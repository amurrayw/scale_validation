library(pcalg)
library(graph)
library(igraph)
library(nFactors)
library(rlist)

##TODO: changing edge weights doesn't produce results that make sense for either FA or FOFC. 
## There must be something going wrong in the igraph->graphNEL part of the project. Continue sims without varying edge coef. for now.

#setwd("~/Dropbox/school/grad. school/gesis/2017/upload/scale_validation/")
source("call_fofc_from_r.R")
source("generate_random_model.R")
source("convert_graph_to_lavaan/convert_graph_to_lavaan.R")

create.dataset.from.file <- function(file="sim_graph_2_lat_pure_measure.r.txt", n=1000, unif.min =1, unif.max=1){

	results <- generate.data.from.dag(read.dag(file=file, unif.min =unif.min, unif.max=unif.max), n=n)

	##Removing latent variables from dataset.
	results[,-grep(names(results), pattern="L")]

}


create.dataset.from.graph <- function(graph, n=1000){

	results <- generate.data.from.dag(graph, n=n)

	##Removing latent variables from dataset.
	results[,-grep(names(results), pattern="L")]

}


build.fofc.model <- function(data, TestType = "TETRAD_WISHART", fofcAlgorithm = "GAP", 
    alpha = .01){

	fofc <- fofc(df=data, TestType, fofcAlgorithm, alpha)	

	nodes <- fofc$nodes


	if(length(fofc$edges)>0){
		edges <- do.call(rbind, lapply(strsplit(fofc$edges, fixed=T, split=" "), function(edgestr){return(edgestr[c(1,3)])}))
	}
	else{edges <- list()}
	adjMat <- matrix(nrow=length(nodes), ncol=length(nodes), data=0)
	adjMat <- data.frame(adjMat)
	names(adjMat) <- nodes

	return(igraph.to.graphNEL(graph.adjacency(as.matrix(assemb.matrix(adjMat, edges)))))
}

##assemb.matrix <- function(empty.mat, edge.list){
##	result <- empty.mat
##	if(length(edge.list)==0){return(result)}
##	for(i in 1:nrow(edge.list)){
##		cause <- which(names(empty.mat)%in%edge.list[i, 1])
##		effect <- which(names(empty.mat)%in%edge.list[i, 2])
##		result[cause, effect] <- 1
##	}
	
##	return(result)
##}


score.fofc <- function(fofc.model, true.graph){
	
	adj.matrix.true <- data.frame(true.graph)
	adj.matrix.fofc <- data.frame(as(fofc.model, "matrix"))
	
	n.latents.truth <- length(grep(pattern="L",
	 names(adj.matrix.true)))

	location.latents.true <- grep(pattern="L", names(adj.matrix.true))


	## Strips out latent-latent edges as FA can't find these.
	true.graph[location.latents.true, location.latents.true] <- 0
	true.graph<-cbind(true.graph, true.graph[,location.latents.true])
	true.graph<-true.graph[,-c(location.latents.true)]
	true.graph<-rbind(true.graph, true.graph[location.latents.true,])
	true.graph<-true.graph[-c(location.latents.true),]

	n.latents.fofc <- length(location.latents.true)

    true.graph<-data.frame(true.graph)


    names(true.graph)[grep(pattern="V", names(true.graph))] <- paste("L", 1:n.latents.truth, sep="")
    rownames(true.graph) <- names(true.graph)   

	true.graph<-igraph.to.graphNEL(graph.adjacency(as(true.graph, "matrix")))

	## Needs to be boolean, not integer for address.mislab function (converting adj.mat to graph fails otherwise).
	fofc.model.backup <- adj.matrix.fofc>0

	#fofc.model <- igraph.to.graphNEL(graph.adjacency(fofc.model))

	# If the mimic.model found has a different numnber of latents than the
	# true graph, then cannot calculate tpr, fpr, tdr. Have therefore treated
	# those as NULL objects. (i.e., as with FA, correct n.latents is assumed)
	if(length(nodes(fofc.model)) == length(nodes(true.graph))){

#		par(mfrow=c(2,2)); plot(fofc.model); plot(true.graph)

		graph.comparison <- compareGraphs(fofc.model, true.graph)
		return(address.mislab(graph.comparison=graph.comparison, fa.model=fofc.model, true.graph=true.graph, fa.model.backup=fofc.model.backup, n.latents=n.latents.fofc, n.calls=n.calls, max.calls=100, best.graph.comparison=best.graph.comparison))
	}
	## Instead of returning null as before, now return a tpr of 0, fpr of 1, and tdr 0. TODO: decide if this makes more sense than returning null.
	return(c(0,1,0))
}

#score.fofc(build.fofc.model(create.dataset.from.file(n=1000, file="graph1.r.txt")), true.graph=as(read.dag("graph1.r.txt"), "matrix"))

score.fa <- function(fa.model, cut.off=.3, true.graph){
	
#	new.order <- permute.latent.col(as.matrix(fa.model), n.latents=ncol(fa.model$loadings))

	fa.mat <- as(fa.model$loadings, "matrix")
	n.latents <- ncol(fa.mat)
	var.names <- row.names(fa.mat)
	n.vars <- nrow(fa.mat)

	fa.model <- prune.fa.paths(fa.model, cut.off=cut.off)

	fa.model.backup <- fa.model	

	location.latents.true <- grep(pattern="^L", x=names(data.frame(true.graph)))

	## Strip out latent-latent edges as FA can't find these.
	true.graph[location.latents.true, location.latents.true] <- 0
	
	true.graph<-cbind(true.graph, true.graph[,location.latents.true])
	true.graph<-true.graph[,-c(location.latents.true)]
	true.graph<-rbind(true.graph, true.graph[location.latents.true,])
	true.graph<-true.graph[-c(location.latents.true),]

    ## Adds in the names of the moved latent variables, which can sometimes be lost.
    true.graph<-data.frame(true.graph)
    names(true.graph)[grep(pattern="V", names(true.graph))] <- paste("L", 1:length(location.latents.true), sep="")
    rownames(true.graph) <- names(true.graph)

	fa.model <- igraph.to.graphNEL(graph.adjacency(fa.model))
	true.graph<-igraph.to.graphNEL(graph.adjacency(as(true.graph, "matrix")))

	if(length(nodes(fa.model))==length(nodes(true.graph))){

		graph.comparison <- (compareGraphs(fa.model, true.graph))

		##Handles the case where, due to a mislabeling (i.e., fa assigns all cluster 2 elements to cluster 1 and vice versa), FA gets more wrong than it should.
		return(address.mislab(graph.comparison, fa.model, true.graph, fa.model.backup, n.latents, n.calls=1, max.calls=permutation(n.latents, n.latents)^2, best.graph.comparison=graph.comparison))
	}
	return(NULL)
}

## Checks every permutation of labels for latent variables in case mislabeling is harming FA unjustly.
address.mislab <- function(graph.comparison, fa.model, true.graph, fa.model.backup, n.latents, n.calls, max.calls=100, best.graph.comparison){

	n.var <- ncol(fa.model.backup)-n.latents

	permutations <- rbind(1:n.latents, permute::allPerms(1:n.latents))+n.var

	comparisons <- t(apply(permutations, 1, function(new.order){compareGraphs(igraph.to.graphNEL(graph.adjacency(fa.model.backup[c(1:n.var, new.order),])), true.graph)}))

	## The third column in the matrix comparisons are true discovery rates. Then select the first max row (in case of a tie).
	best.graph.comparison <- comparisons[which(comparisons[,3] == max(comparisons[,3]))[1], ]
##plot(fa.model)
	return(best.graph.comparison)
}


#score.fa(factanal(x=create.dataset(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.5)


## Keeps trying until the true discovery rate is as high as possible (or if the max number of attempts at maximization has been reached).
##address.mislab <- function(graph.comparison,  fa.model=fa.model, true.graph=true.graph, fa.model.backup=fa.model.backup, n.latents=n.latents, n.calls=1, max.calls=100, best.graph.comparison){
##		## The third element in the vector graph.comparison is the true discovery rate. 
##		if(graph.comparison[3]>best.graph.comparison[3]){best.graph.comparison <- graph.comparison}


##		if(n.calls<max.calls && graph.comparison[3]<=best.graph.comparison[3]){
##				new.order <- permute.latent.col(fa.model.backup, n.latents=n.latents)
##				fa.model <- igraph.to.graphNEL(graph.adjacency(fa.model.backup[,new.order]))
##			return(address.mislab(compareGraphs(fa.model, true.graph), fa.model, true.graph, fa.model.backup, n.latents, n.calls+1, best.graph.comparison=best.graph.comparison))
##		}
##		else{
##			return(best.graph.comparison)
##		}
##}


## Need to do this so that a misslabeling won't harm results of FA.
##permute.latent.col <- function(fa.model, n.latents){

##	n.var <- ncol(fa.model)-n.latents
##	permute <- sample((n.var+1):ncol(fa.model), replace=FALSE, size=n.latents)

##	new.order <- c(1:n.var, permute)

##	return(new.order)
##}

##permutation = function(n, x) {
##  return(factorial(n)/factorial(n-x))
##}


prune.fa.paths <- function(fa.model, cut.off=.3){
	
	fa.loadings.matrix <- as(fa.model$loadings, "matrix")
	
	n.latents <- ncol(fa.loadings.matrix)
	var.names <- row.names(fa.loadings.matrix)
	n.vars <- length(var.names)
	
	adj.mat <- matrix(FALSE, nrow=(n.vars+n.latents),
	 ncol=(n.vars+n.latents), dimnames=list("row"=c(var.names,
		 1:n.latents), "col"=c(var.names, 1:n.latents)))
	
	for(i in 1:n.latents){

		adj.mat[abs(fa.model$loadings[,i])>cut.off, 
		i+length(var.names)] <- TRUE
	}
	
	adj.mat[(length(var.names)+1):(length(var.names)+n.latents),
	 (length(var.names)+1):(length(var.names)+n.latents)] <- FALSE
	
##Return transpose of matrix so row/column cause-effect agrees with format used by igraph/graphNEL.
	return(t(adj.mat))
	
}


# converts the Tetrad .r.txt representation to a graphNEL object.
read.dag <- function(file, unif.min =1, unif.max=1){
	orig.mat <- read.table(file=file)
	
	orig.mat[orig.mat==1] <-0
	orig.mat[orig.mat==-1] <-1


	ind <- which(orig.mat==1, arr.ind=TRUE)
	## Allows edge weights to vary randomly (drawn from uniform distribution).
	for(i in 1:nrow(ind)){
		orig.mat[ind[i,1], ind[i,2]] <- runif(n=1, unif.min, unif.max)
	}

#	print(E(graph.adjacency(as.matrix(orig.mat), weighted=TRUE))$weight)

	

	final.graph <- igraph.to.graphNEL(graph.adjacency(as.matrix(orig.mat), weighted=TRUE))
	return(final.graph)
}

#read.dag("graph1.r.txt")


# Takes in a graphNEL object, and generates normally distributed data from it.
generate.data.from.dag <- function(graph, n=100, errDist="normal"){	

	top.sort <- topological.sort(igraph.from.graphNEL(graph))
	var.names <- nodes(graph)[top.sort]

	graph<-igraph.to.graphNEL(graph.adjacency(as(graph, "matrix")[var.names,
	 var.names], weighted=TRUE))
##	plot(graph)
	generated.data<-data.frame(rmvDAG(dag=graph, n=n, errDist=errDist))
	names(generated.data) <- var.names
	
	return(generated.data)
}

##generate.data.from.dag(read.dag(file="sim_graph_2_lat_pure_measure.r.txt", unif.min =.5, unif.max=1))










### Test portion for CFA. lavaan converter seems to work.

set.seed(123)
tmp.graph <- generate.measurement.model(n.latents=2)
sem(convert.igraph.to.lavaan(igraph.from.graphNEL(tmp.graph)), generate.data.from.dag(tmp.graph))

## Based, on plot, model seems to be read in correctly.
semPlot::semPaths(sem(convert.igraph.to.lavaan(igraph.from.graphNEL(tmp.graph)), generate.data.from.dag(tmp.graph)))

###Test to see if chi-square test in lavaan can distinguish models with l1->l2 edge. 
## Note, sample size is only 100.
##false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph)))@test[[1]]$pvalue<=.05, n=1000))

##true.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+L2", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph)))@test[[1]]$pvalue<=.05, n=1000))

##> false.model
##[1] 71
##> true.model
##[1] 67

## Evidently, despite the false model implying that e.g., X4 is indep. of X1, chi-square rejects model about as often as true model. I hope I'm interpreting something wrong, as this is not a good result for CFA. As sample size was very small (100 obs), will try again with 1000 obs.
##false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph, n=1000)))@test[[1]]$pvalue<=.05, n=1000))

##true.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+L2", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph, n=1000)))@test[[1]]$pvalue<=.05, n=1000))


##> false.model
##[1] 57
##> true.model
##[1] 51

##false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph, n=10000)))@test[[1]]$pvalue<=.05, n=1000))

##true.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+L2", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph, n=10000)))@test[[1]]$pvalue<=.05, n=1000))

##> false.model
##[1] 62
##> true.model
##[1] 53

##Turns out these results are due to lavaan treating latent variables as correlated by default. Need to set: orthogonal=TRUE in order to turn this behavior off.

false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph,n=10000),orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

true.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+L2", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

##At last, proper behavior! So, it looks like independence facts still matter.

##> false.model
##[1] 1000
##> true.model
##[1] 48


##Lets try one where the false model wrongly swaps two variables between clusters (X0 and X4), but no indep. is violated.

false.model <- sum(replicate((sem(c( "L1 =~ X4+X1+X2+X3+L2", "L2 =~ X0+X5+X6+X7"), generate.data.from.dag(tmp.graph),n=10000,orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

true.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+L2", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

##> false.model
##[1] 990
##> true.model
##[1] 44

## Interesting, the rejection doesn't seem to be due to independence.

## maybe result was due to first measure of each latent being set to 1 (so as to scale)?
false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X7+L2", "L2 =~ X4+X5+X6+X3"), generate.data.from.dag(tmp.graph,n=10000),orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

##> false.model
##[1] 1000

#Nope! Seems to be able to pick up some cluster errors.

## Test edge direction error: L2-> L1, not L1->L2

false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3", "L2 =~ X4+X5+X6+X7+L1"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

## Doesn't seem to be able to figure out (some) errors in latent-latent direction.
##> false.model
##[1] 52


## Clustering error where messes up 2 variables:

false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X6+X7+L2", "L2 =~ X4+X5+X2+X3"), generate.data.from.dag(tmp.graph,n=10000),orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))
##> false.model
##[1] 1000


## False impurity model:

false.model <- try(sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+L2", "L2 =~ X4+X5+X6+X7", "X4~X5"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000)))

## Doesn't seem able to detect false measure measure edge within the same cluster.
##> false.model
##[1] 61

## False impurity, different clusters.
false.model <- try(sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+L2", "L2 =~ X4+X5+X6+X7", "X4~X1"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000)))

## Doesn't seem able to detect false measure measure edge between clusters.
##> false.model
##[1] 54


## Clustering error where one variable is assigned to two clusters.
false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+L2", "L2 =~ X4+X5+X6+X7+X2"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

## Doesn't seem to be able to pick up this kind of error.
##> false.model
##[1] 54

## Clustering error where one variable is assigned to two clusters (other cluster has problem)
false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+X6+L2", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

## Which cluster the error is made with doesn't seem to matter.
##> false.model
##[1] 47


## Clustering error where two variables are assigned to two clusters.
false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+X6+X7+L2", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

## Doesn't seem able to detect the false double clustering.
##> false.model
##[1] 42


## Clustering error where all variables from one clsuter are assigned to two clusters.
false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+X4+X5+X6+X7+L2", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

## doesn't seem able to detect this, but note that produced 50+ warnings stating:
##50: In lav_model_vcov(lavmodel = lavmodel, lavsamplestats = lavsamplestats,  ... :
##  lavaan WARNING: could not compute standard errors!
##  lavaan NOTE: this may be a symptom that the model is not identified.
##> false.model
##[1] 65

## Since the model is evidently not identifiable (I assume the edge parameters can't be narrowed down to a finite set), this error isn't too serious, as the fitting program will warn if it happens.

## Trying case where latent-latent edge is mistakenly omitted but all latents of L2 are also assigned to L1.
false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+X4+X4+X5+X6+X7", "L2 =~ X4+X5+X6+X7"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

## Looks like chi-square also can't detect this situation.
##> false.model
##[1] 50



## False two factor model:
false.model <- sum(replicate((sem(c( "L1 =~ X0+X1+X2+X3+X4+X5+X6+X7"), generate.data.from.dag(tmp.graph,n=10000), orthogonal=TRUE))@test[[1]]$pvalue<=.05, n=1000))

## Seems able to pick up this sort of error.
##> false.model
##[1] 1000





### So argument would be: EFA doesn't work well for det. number of latents, clustering with promax is somewhat worse than FOFC (though still need to try adding impurities). CFA is just SEM. Assuming number of latents is correct, preliminary results indicate can pick up false indep. claims and some misclustering, but cannot pick up (non-dep/indep adding) errors in latent-latent direction. This is all based on using chi-square test statistic (the p-value) to decide accept/reject model. Should then argue that such a method suffers from computational explosion (number of alternative models capable of being compared is large).  Note: This objection requires the misclustering detection (or impurity detection) to be weak, which it seems to be.

##Since  chi-square test statistics compares the null covariance matrix to the observed cov.mat, these results make sense. Indep is strongest, followed by by changes in the graph that add latent-latent edge coefficients to the path (i.e., tetrad const. violation stuff). 












