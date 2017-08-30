

generate.measurement.model <- function(n.latents=1, n.measures.per.latent=4, n.impurities=list(n.sharing.latent=0, n.output.output.edges=0), n.latent.latent.edges=0){

	created.latents <- list()
	for(i in 1:n.latents){
		created.latents[[i]] <- create.latent(latent.name=paste("L", i, sep=""),n.measures=n.measures.per.latent, naming.index=n.measures.per.latent*(i-1))
	}


	
	
	latent.edges.to.add <- rbind(latent.edges.to.add, )


	return(created.latents)



}

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


add.impurities <- function(latent.1, latent.2, shared.latent, output.output.edge){}




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



