library(pcalg)
library(graph)
library(igraph)

setwd("~/Dropbox/school/grad. school/gesis/2017/upload/scale_validation/")


#score.fa(factanal(x=create.dataset(), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"))


hist(replicate(score.fa(factanal(x=create.dataset(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.3), n=1000)[3,])


#set.seed(100);prune.fa.paths(factanal(x=create.dataset(n=1000), factors=2), cut.off=.3)[, c(10,9)]


create.dataset <- function(file="sim_graph_2_lat_pure_measure.r.txt", n=1000){

	results <- generate.data.from.dag(read.dag(file=file), n=n)

	##Removing latent variables from dataset.
	results[,-grep(names(results), pattern="L")]

}

build.fofc model <- function(){


}

## TODO: Implement this (is still the old MIMIC scoring function).
score.fofc <- function(fofc.model, true.graph){
	
	adj.matrix.true <- data.frame(true.graph)
	
	adj.matrix.mimic <- data.frame(as(mimic.model, "matrix"))
	
	n.latents.truth <- grep(pattern="[L:digit:]",
	 names(adj.matrix.true))

	n.latents.mimic <- grep(pattern="[L:digit:]",
	 names(adj.matrix.mimic))
	
	true.graph<-igraph.to.graphNEL(graph.adjacency(true.graph))


	# If the mimic.model found has a different numnber of latents than the
	# true graph, then cannot calculate tpr, fpr, tdr. Have therefore treated
	# those as NULL objects. (i.e., as with FA, correct n.latents is assumed)
	if(length(nodes(mimic.model)) == length(nodes(true.graph))){
		return(compareGraphs(mimic.model, true.graph))
	}
	return(NULL)
}

score.fa <- function(fa.model, cut.off=.3, true.graph){
	
#	new.order <- permute.latent.col(as.matrix(fa.model), n.latents=ncol(fa.model$loadings))

	fa.mat <- as(fa.model$loadings, "matrix")
	n.latents <- ncol(fa.mat)
	var.names <- row.names(fa.mat)
	n.vars <- nrow(fa.mat)

	fa.model <- prune.fa.paths(fa.model, cut.off=cut.off)

	fa.model.backup <- fa.model	

	fa.model <- igraph.to.graphNEL(graph.adjacency(fa.model))
	true.graph <- igraph.to.graphNEL(graph.adjacency(true.graph))

	if(length(nodes(fa.model))==length(nodes(true.graph))){

		graph.comparison <- (compareGraphs(fa.model, true.graph))

		##Handles the case where, due to a mislabeling (i.e., fa assigns all cluster 2 elements to cluster 1 and vice versa), FA gets everthing wrong instead of everything right.
		## The third element in the vector graph.comparison is the true discovery rate. Note that this fix only solves the problem for a pair of latents. A triple will require something more elaborate. 

		return(address.mislab(graph.comparison, fa.model, true.graph, fa.model.backup, n.latents))
	}
	return(NULL)
}

## Keeps trying until the true discovery rate is greater than 0.
address.mislab <- function(graph.comparison,  fa.model=fa.model, true.graph=true.graph, fa.model.backup=fa.model.backup, n.latents=n.latents){
		if(graph.comparison[3]==0){
				new.order <- permute.latent.col(fa.model.backup, n.latents=n.latents)
				fa.model <- igraph.to.graphNEL(graph.adjacency(fa.model.backup[,new.order]))
			return(address.mislab(compareGraphs(fa.model, true.graph), fa.model, true.graph, fa.model.backup, n.latents))
		}
		else{
			return(graph.comparison)
		}
}

#score.fa(factanal(x=create.dataset(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.5)


## Need to do this so that a misslabeling won't harm results of FA.
permute.latent.col <- function(fa.model, n.latents){

	n.var <- ncol(fa.model)-n.latents
	permute <- sample((n.var+1):ncol(fa.model), replace=FALSE, size=n.latents)

	new.order <- c(1:n.var, permute)

	return(new.order)
}


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
	
	return(adj.mat)
	
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
	
	graph<-igraph.to.graphNEL(graph.adjacency(as(graph, "matrix")[var.names,
	 var.names]))
	
	generated.data<-data.frame(rmvDAG(dag=graph, n=n, errDist=errDist))
	names(generated.data) <- var.names
	
	return(generated.data)
}
