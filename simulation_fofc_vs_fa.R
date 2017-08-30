library(pcalg)
library(graph)
library(igraph)
library(nFactors)
library(rlist)

##TODO: changing edge weights doesn't produce results that make sense for either FA or FOFC. 
## There must be something going wrong in the igraph->graphNEL part of the project. Continue sims without varying edge coef. for now.

#setwd("~/Dropbox/school/grad. school/gesis/2017/upload/scale_validation/")
source("call_fofc_from_r.R")



create.dataset <- function(file="sim_graph_2_lat_pure_measure.r.txt", n=1000, unif.min =1, unif.max=1){

	results <- generate.data.from.dag(read.dag(file=file, unif.min =unif.min, unif.max=unif.max), n=n)

	##Removing latent variables from dataset.
	results[,-grep(names(results), pattern="L")]

}


build.fofc.model <- function(data, TestType = "TETRAD_WISHART", fofcAlgorithm = "GAP", 
    alpha = .01){

	fofc <- fofc(loadContinuousData(df=data), TestType, fofcAlgorithm, alpha)	

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

assemb.matrix <- function(empty.mat, edge.list){
	result <- empty.mat
	if(length(edge.list)==0){return(result)}
	for(i in 1:nrow(edge.list)){
		cause <- which(names(empty.mat)%in%edge.list[i, 1])
		effect <- which(names(empty.mat)%in%edge.list[i, 2])
		result[cause, effect] <- 1
	}
	
	return(result)
}


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

	true.graph<-igraph.to.graphNEL(graph.adjacency(true.graph))

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


#score.fofc(build.fofc.model(create.dataset(n=1000, file="graph1.r.txt")), true.graph=as(read.dag("graph1.r.txt"), "matrix"))



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


	fa.model <- igraph.to.graphNEL(graph.adjacency(fa.model))
	true.graph <- igraph.to.graphNEL(graph.adjacency(true.graph))

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
	 var.names]))
##	plot(graph)
	generated.data<-data.frame(rmvDAG(dag=graph, n=n, errDist=errDist))
	names(generated.data) <- var.names
	
	return(generated.data)
}




