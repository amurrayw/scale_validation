
library(pcalg)
library(igraph)
library(lavaan)

## TODO: Given a graph, convert to lavaan object. (guide to lavaan at: https://m-clark.github.io/docs/sem/graphical-models-1.html)

## Can use topo_sort function to find order that should be used to construct model.
## Need to write: get parent/descendent functions.

## for a given node, parse node names (variable and its parents) & structure into lavaan sem language. 
## E.g., x->y<-z, define "x, y, y~x+z". 
## Note for latents, use "y=~x+z"



## Need to process nodes in topological order.

convert.igraph.to.lavaan <- function(graph.obj){


    adj.matrix <- as.matrix(as_adjacency_matrix(graph.obj))

    node.names <- get.node.names(graph.obj)

    graph.order <- topo_sort(graph.obj)

    lavaan.model <- c()

    for(node.index in graph.order){
        current.node <- node.names[node.index]
        node.parents <- node.names[get.parents(node.index, adj.matrix)]
        
 
        if(length(node.parents)==0){
            next
            #lavaan.model <- c(lavaan.model, paste(current.node, "\n"))
        }
        else{
            lavaan.model <- c(lavaan.model, paste(current.node, "=~", paste(node.parents, sep="", collapse="+"), collapse=""))
    
        }
    }

    return(lavaan.model)

}

#convert.igraph.to.lavaan(g)

#convert.igraph.to.lavaan(igraph.from.graphNEL(randomDAG(n=10, prob=.5)))



#set.seed(123)
#lavaanify(convert.igraph.to.lavaan(igraph.from.graphNEL(randomDAG(n=3, prob=.5, V=c("X", "Y", "Z")))))


#####Helper Methods#####

get.parents <- function(current.node, adj.matrix){
#print(adj.matrix)
    return(which(adj.matrix[current.node,]>0))
}

get.node.names <- function(graph.obj){
    return(as.character(V(graph.obj)$name))
}



