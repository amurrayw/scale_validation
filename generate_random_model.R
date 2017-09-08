

generate.measurement.model <- function(n.latents=1, n.measures.per.latent=4, n.impurities=list(n.sharing.latent=0, n.output.output.edges=0), n.latent.latent.edges=0){

	created.latents <- list()
	for(i in 1:n.latents){
		created.latents[[i]] <- create.latent(latent.name=paste("L", i, sep=""),n.measures=n.measures.per.latent, naming.index=n.measures.per.latent*(i-1))
	}

	possible.permutations <- produce.permutations(1:n.latents, rand.remove.mirrors=TRUE)

	chosen.permutations <- possible.permutations[sample(1:nrow(possible.permutations), replace=FALSE, size=n.latent.latent.edges), ]

	if(is.vector(chosen.permutations)){

		created.latents[[chosen.permutations[1]]] <- connect.latents(created.latents[[chosen.permutations[1]]], created.latents[[chosen.permutations[2]]])
	}
	else{
		for(i in 1:nrow(chosen.permutations)){
			if(nrow(chosen.permutations)==0){break}
			cause <- chosen.permutations[i, 1]
			effect <- chosen.permutations[i, 2]

			created.latents[[cause]] <- connect.latents(created.latents[[cause]], created.latents[[effect]])

		}
	}

	master.nodes <- unique(unlist(lapply(created.latents, function(cluster){return(cluster$nodes)})))
	master.edges <- unlist(lapply(created.latents, function(cluster){return(cluster$edges)}))

	result.graph <- convert.list.to.graph(list(nodes=master.nodes, edges=master.edges))

	## Graph has a cycle if can't tsort, tsort won't return all nodes, so length diff. implies regen. graph. Will have to see if performance is good enough as complexity of model increases.
	if(length(topological.sort(igraph.from.graphNEL(result.graph)))!=length(nodes(result.graph))){
		return(generate.measurement.model(n.latents=n.latents, n.measures.per.latent=n.measures.per.latent, n.impurities=n.impurities, n.latent.latent.edges=n.latent.latent.edges))
	}
	else{return(result.graph)}



}

#generate.measurement.model(2, n.measures.per.latent=4)
#generate.measurement.model(2, n.measures.per.latent=4, n.latent.latent.edges=1)
#generate.measurement.model(3, n.measures.per.latent=4, n.latent.latent.edges=1)
#plot(generate.measurement.model(3, n.measures.per.latent=4, n.latent.latent.edges=2))


create.latent <- function(latent.name="L1", n.measures=4, naming.index=1){

	measure.var <- c()
	for(i in naming.index:(naming.index+n.measures-1)){
		measure.var <- c(measure.var, paste("X", i, sep=""))
	} 
	
	edges.var <- c()
	for(i in measure.var){edges.var <- c(edges.var, create.edge(latent.name, i))}
	
	return(list(nodes=c(latent.name, measure.var), edges=edges.var))
	
}


connect.latents <- function(latent.cause=create.latent("L1"), latent.effect=create.latent("L2", naming.index=5)){

	## As the first variable created in a latent is it's own name, the first match for ^L is the name of the latent.
	name.cause <- latent.cause$nodes[grep(latent.cause$nodes, pattern="^L")[1]]
	name.effect <- latent.effect$nodes[grep(latent.effect$nodes, pattern="^L")[1]]

	latent.cause$edges <- c(latent.cause$edges, create.edge(name.cause, name.effect))
	
	latent.cause$nodes <- c(latent.cause$nodes, name.effect)

	return(latent.cause)
}


add.impurities <- function(latent.1, latent.2, shared.latent, output.output.edge){




}




create.edge <- function(cause, effect){return(paste(c(cause, " --> ", effect), sep="", collapse=""))}

convert.list.to.graph <- function(latent.structure=create.latent()){
	

	nodes <- latent.structure$nodes


	if(length(latent.structure$edges)>0){
		edges <- do.call(rbind, lapply(strsplit(latent.structure$edges, fixed=T, split=" "), function(edgestr){return(edgestr[c(1,3)])}))
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


produce.permutations <- function(vector, choose=2, rand.remove.mirrors=FALSE){
	perm <- c()
	for(i in vector){
		for(j in vector){
			if(i==j){next}
			else{perm <- rbind(perm, c(i, j))}
		}
	}

	## Removes cyclic inducing edges from being present at the same time (e.g., 1->3 and 3->1).
	might.remove <- c()
	if(rand.remove.mirrors){
		for(i in 1:nrow(perm)){
			for(j in 1:nrow(perm)){
				if(i>=j){next}
				else if(sum(perm[i,] %in% perm[j,])==choose){might.remove <- c(might.remove, sample(c(i, j), size=1))}
			}			
		}
		return(perm[-unique(might.remove),])
	}
	return(perm)
}




