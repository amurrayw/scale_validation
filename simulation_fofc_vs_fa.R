library(pcalg)
library(graph)
library(igraph)
library(nFactors)
library(rlist)

#setwd("~/Dropbox/school/grad. school/gesis/2017/upload/scale_validation/")
source("call_fofc_from_r.R")

#score.fa(factanal(x=create.dataset(), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"))


##table(replicate(unlist(nScree(x=create.dataset(n=1000), model="factors")$Components), n=1000))


table(replicate(as.numeric(names(which.max(table(unlist(nScree(x=create.dataset(n=1000), model="factors")$Components))))), n=1000))


table(replicate(as.numeric(names(which.max(table(unlist(nScree(x=create.dataset(n=1000, file="graph1.r.txt"), model="factors")$Components))))), n=1000))

par(mfrow=c(2,2))
hist(replicate(score.fa(factanal(x=create.dataset(n=1000, file="graph1.r.txt"), factors=3), true.graph=as(read.dag("graph1.r.txt"), "matrix"), cut.off=.3), n=1000)[3,])


hist(replicate(score.fofc(build.fofc.model(create.dataset(n=1000, file="graph1.r.txt")), true.graph=as(read.dag("graph1.r.txt"), "matrix")), n=1000)[3,])



##TODO: Need to also account for getting rotation correct (promax is supposed to be used when latents are "correlated" with one another). Should test with impure measures (i.e., correlation is through observed variables, not latents).
 plot(igraph.to.graphNEL(graph.adjacency(prune.fa.paths(factanal(x=create.dataset(n=1000, file="graph1.r.txt"), factors=3, rotation="promax"), .3))))

score.fa(factanal(x=create.dataset(n=1000, file="graph1.r.txt"), factors=3, rotation="promax"), true.graph=as(read.dag("graph1.r.txt"), "matrix"), cut.off=.3) 



hist(replicate(score.fa(factanal(x=create.dataset(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.3), n=1000)[3,])
hist(replicate(score.fa(factanal(x=create.dataset(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.5), n=1000)[3,])
hist(replicate(score.fa(factanal(x=create.dataset(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.7), n=1000)[3,])





handle.null.fofc <- function(fofc.sim.object){

	if(class=="matrix"){return(fofc.sim.object)}
	else{
		fofc.sim.object <- list.clean(fofc.sim.object, is.null)
				
	}
}







#hist(replicate(score.fa(factanal(x=create.dataset(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.3), n=1000)[3,])


#set.seed(100);prune.fa.paths(factanal(x=create.dataset(n=1000), factors=2), cut.off=.3)[, c(10,9)]





create.dataset <- function(file="sim_graph_2_lat_pure_measure.r.txt", n=1000){

	results <- generate.data.from.dag(read.dag(file=file), n=n)

	##Removing latent variables from dataset.
	results[,-grep(names(results), pattern="L")]

}


build.fofc.model <- function(data, TestType = "TETRAD_WISHART", fofcAlgorithm = "GAP", 
    alpha = .01){

	fofc <- fofc(loadContinuousData(df=data), TestType, fofcAlgorithm, alpha)	

	nodes <- fofc$nodes

	#print(nodes)
	edges <- do.call(rbind, lapply(strsplit(fofc$edges, fixed=T, split=" "), function(edgestr){return(edgestr[c(1,3)])}))
##
#	print(edges)


	adjMat <- matrix(nrow=length(nodes), ncol=length(nodes), data=0)
	adjMat <- data.frame(adjMat)
	names(adjMat) <- nodes

	return(igraph.to.graphNEL(graph.adjacency(as.matrix(assemb.matrix(adjMat, edges)))))
}

assemb.matrix <- function(empty.mat, edge.list){
	result <- empty.mat

	for(i in 1:nrow(edge.list)){
		cause <- which(names(empty.mat)%in%edge.list[i, 1])
		effect <- which(names(empty.mat)%in%edge.list[i, 2])
		result[cause, effect] <- 1
	}
	
	return(result)
}


## TODO: Need to fix scoring. Currently there seems to be an issue with the way latents are named. Both search model and true graph need to use same names for all nodes. Should also change/check this in FA scoring function.
score.fofc <- function(fofc.model, true.graph){
	
	adj.matrix.true <- data.frame(true.graph)
	
	adj.matrix.fofc <- data.frame(as(fofc.model, "matrix"))
	
	n.latents.truth <- length(grep(pattern="L",
	 names(adj.matrix.true)))

	location.latents.true <- grep(pattern="L", names(adj.matrix.true))

	## Strip out latent-latent edges as FA can't find these.
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

##print(true.graph)




	location.latents.true <- grep(pattern="^L", x=names(data.frame(true.graph)))
#print(location.latents.true)
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
##print("fa.model.backup")
##print(fa.model.backup)

##print("true.graph")
##print(as(true.graph, "matrix"))

	## The third column in the matrix comparisons are true discovery rates. Then select the first max row (in case of a tie).
	best.graph.comparison <- comparisons[which(comparisons[,3] == max(comparisons[,3]))[1], ]
plot(fa.model)
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
read.dag <- function(file){
	orig.mat <- read.table(file=file)
	
	orig.mat[orig.mat==1] <-0
	orig.mat[orig.mat==-1] <-1
	

	final.graph <- igraph.to.graphNEL(graph.adjacency(as.matrix(orig.mat)))

	return(final.graph)
}


# Takes in a graphNEL object, and generates normally distributed data from it.
generate.data.from.dag <- function(graph, n=100, errDist="normal"){	

	top.sort <- topological.sort(igraph.from.graphNEL(graph))
	var.names <- nodes(graph)[top.sort]
	
	##TODO: Let edge weights vary (currently all are 1). Error distribution is always std. normal.
	## Idea: if adj.mat has a 1, replace with runif(0.001,1,1).
	graph<-igraph.to.graphNEL(graph.adjacency(as(graph, "matrix")[var.names,
	 var.names]))
	
	generated.data<-data.frame(rmvDAG(dag=graph, n=n, errDist=errDist))
	names(generated.data) <- var.names
	
	return(generated.data)
}




