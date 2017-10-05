library(rJava)
library(rcausal)


fofc <- function(df, TestType = "TETRAD_WISHART", fofcAlgorithm = "GAP", 
    alpha = .01, java.parameters = NULL){
     
    df <- loadContinuousData(df)

    if(!is.null(java.parameters)){
        options(java.parameters = java.parameters)
        params <- c(java.parameters = java.parameters)
    }

    fofcAlgorithmPath <- "edu/cmu/tetrad/search/FindOneFactorClusters"
    TestTypePath <- "edu/cmu/tetrad/search/TestType"

    ## Instantiate fofc object.
#    fofc_instance <- .jnew("edu/cmu/tetrad/search/FindOneFactorClusters", df, J(TestTypePath)$TETRAD_WISHART, J(fofcAlgorithmPath)$Algorithm$GAP, as.double(alpha))

    fofc_instance <- .jnew("edu/cmu/tetrad/search/FindOneFactorClusters", df, .jfield(TestTypePath,name=TestType), .jfield(paste(fofcAlgorithmPath, "$Algorithm", sep=""), name=fofcAlgorithm), as.double(alpha))

    ## Search
    fofc_graph <- .jcall(fofc_instance, "Ledu/cmu/tetrad/graph/Graph;", "search")

    ##List to contain results of fofc search.
    fofc <- list()

    if(!is.null(e <- .jgetEx())){
        .jclear()
        fofc$nodes <- colnames(df)
        fofc$edges <- NULL
        print("Java exception was raised")
        print(e)
    }else{
        V <- extractTetradNodes(fofc_graph)

        ## Get nodes.
        fofc$nodes <- V
        
        ## extract edges
        fofc_edges <- extractTetradEdges(fofc_graph)
        
        fofc$edges <- fofc_edges
    }
    
    return(fofc)
}

## Data Frame to Tetrad Dataset
#data<-read.csv(file="~/Dropbox/school/grad. school/gesis/2017/upload/scale_validation/depression_scale.csv")[,-1]
#data<-data[complete.cases(data),]

#fofc(loadMixedData(df=data, 0))




#####Helper methods

############## From R-causal https://github.com/bd2kccd/r-causal ###############
loadContinuousData <- function(df){
	node_names <- colnames(df)
	node_list <- .jnew("java/util/ArrayList")
	for (i in 1:length(node_names)){
		nodname <- .jnew("java/lang/String", node_names[i])
		nodi <- .jnew("edu/cmu/tetrad/data/ContinuousVariable", nodname)
		node_list$add(nodi)
	}
	node_list <- .jcast(node_list, "java/util/List")
	mt <- as.matrix(df)
	mat <- .jarray(mt, dispatch=TRUE)
    
    data <- .jnew("edu/cmu/tetrad/data/DoubleDataBox", mat)
    data <- .jcast(data, "edu/cmu/tetrad/data/DataBox")
    boxData <- .jnew("edu/cmu/tetrad/data/BoxDataSet", 
                            data, node_list)
    boxData <- .jcast(boxData, "edu/cmu/tetrad/data/DataSet")
	return(boxData)
}

loadMixedData <- function(df, numCategoriesToDiscretize = 4){
    node_names <- colnames(df)
    cont_list <- c()
    disc_list <- c()
    node_list <- .jnew("java/util/ArrayList")
    for (i in 1:length(node_names)){
        nodname <- .jnew("java/lang/String", node_names[i])
        cate <- unique(df[[node_names[i]]])
        if(length(cate) > numCategoriesToDiscretize){
            # Continuous variable
            nodi <- .jnew("edu/cmu/tetrad/data/ContinuousVariable", nodname)
            node_list$add(nodi)
            
            cont_list <- c(cont_list, i)
        }else{
            # Discrete variable
            cate <- sort(cate)
            cate_list <- .jnew("java/util/ArrayList")
            for(j in 1:length(cate)){
                cate_list$add(as.character(cate[j]))
            }
            cate_list <- .jcast(cate_list, "java/util/List")
            nodi <- .jnew("edu/cmu/tetrad/data/DiscreteVariable",
            nodname, cate_list)
            node_list$add(nodi)
            
            # Substitute a new categorial value
            cate <- data.frame(cate)
            new_col <- sapply(df[,i],function(x,cate)
            as.integer(which(cate[,1] == x)),cate=cate)
            new_col = as.integer(new_col - 1)
            df[,i] <- (data.frame(new_col))[,1]

            disc_list <- c(disc_list, i)
        }
    }
    
    node_list <- .jcast(node_list, "java/util/List")
    mixedDataBox <- .jnew("edu/cmu/tetrad/data/MixedDataBox", node_list,as.integer(nrow(df)))
    
    for(row in 1:nrow(df)){
        # print(paste("row:",row,sep=" "))
        if(length(cont_list) > 0){
            for(j in 1:length(cont_list)){
                col <- cont_list[j]
                # print(paste("col:",col,sep=" "))
                value <- as.character(df[row,col])
                #print(value)
                value <- .jnew("java/lang/Double", value)
                value <- .jcast(value, "java/lang/Number")
                mixedDataBox$set(as.integer(row-1),as.integer(col-1),value)
            }
        }
        if(length(disc_list) > 0){
            for(j in 1:length(disc_list)){
                col <- disc_list[j]
                # print(paste("col:",col,sep=" "))
                value <- as.character(df[row,col])
                # print(value)
                value <- .jnew("java/lang/Integer", value)
                value <- .jcast(value, "java/lang/Number")
                mixedDataBox$set(as.integer(row-1),as.integer(col-1),value)
            }
        }
    }
    
    data <- .jcast(mixedDataBox, "edu/cmu/tetrad/data/DataBox")
    boxData <- .jnew("edu/cmu/tetrad/data/BoxDataSet",
                data, node_list)
    boxData <- .jcast(boxData, "edu/cmu/tetrad/data/DataSet")
    return(boxData)
}
extractTetradNodes <- function(resultGraph){
    nods <- resultGraph$getNodes()
	V <- sapply(as.list(nods), with, toString())
    return(V)
}

extractTetradEdges <- function(resultGraph){
    eds <- resultGraph$getEdges()
    fgs_edges <- c()
    if(!is.null(eds)){
	   fgs_edges <- sapply(as.list(eds), .jrcall, "toString")
    }
    return(fgs_edges)
}

rCovMatrix2TetradCovMatrix <- function(covmat, node_list, sample_size){
  mat <- .jarray(covmat, dispatch=TRUE)
  tetmat <- .jnew("edu/cmu/tetrad/util/TetradMatrix", mat)
  tetcovmat <- .jnew("edu/cmu/tetrad/data/CovarianceMatrix", node_list, 
                        tetmat, as.integer(sample_size))
  tetcovmat <- .jcast(tetcovmat, "edu/cmu/tetrad/data/ICovarianceMatrix", 
                      check=TRUE)
  return(tetcovmat)
}

dataFrame2TetradConditionalGaussianScore <- function(df,
    numCategoriesToDiscretize = 4, penaltydiscount = 4, structurePrior = 1.0){
    boxData <- loadMixedData(df, numCategoriesToDiscretize)
    score <- .jnew("edu/cmu/tetrad/search/ConditionalGaussianScore",
                boxData, structurePrior, TRUE)
    score$setPenaltyDiscount(penaltydiscount)
    score <- .jcast(score, "edu/cmu/tetrad/search/Score")
    return(score)
}

dataFrame2TetradSemBicScore <- function(df,penaltydiscount = 4.0){
    boxData <- loadContinuousData(df)
    covMat <- .jnew("edu/cmu/tetrad/data/CovarianceMatrixOnTheFly", boxData)
    covMat <- .jcast(covMat, "edu/cmu/tetrad/data/ICovarianceMatrix")
    score <- .jnew("edu/cmu/tetrad/search/SemBicScore", covMat)
	score$setPenaltyDiscount(penaltydiscount)
    score <- .jcast(score, "edu/cmu/tetrad/search/Score")
    return(score)
}



dataFrames2TetradSemBicScoreImages <- function(dfs,penaltydiscount = 4.0){
    datasets <- .jnew("java/util/ArrayList")
    for(i in 1:length(dfs)){
        df <- dfs[[i]]
        boxData <- loadContinuousData(df)
        dataModel <- .jcast(boxData, "edu/cmu/tetrad/data/DataModel")
        datasets$add(dataModel)
        
    }
    datasets <- .jcast(datasets, "java/util/List")
    score <- .jnew("edu/cmu/tetrad/search/SemBicScoreImages", datasets)
    score$setPenaltyDiscount(penaltydiscount)
    score <- .jcast(score, "edu/cmu/tetrad/search/Score")
    return(score)
}

dataFrame2TetradBDeuScore <- function(df,structurePrior = 1.0, 
    samplePrior = 1.0){
    boxData <- loadDiscreteData(df)
    score <- .jnew("edu/cmu/tetrad/search/BDeuScore", boxData)
    score$setStructurePrior(as.double(structurePrior))
    score$setSamplePrior(as.double(samplePrior))
    score <- .jcast(score, "edu/cmu/tetrad/search/Score")
    return(score)
}


ugraphToTetradGraph <- function(ugmat, node_list){
  numNodes <- ncol(ugmat)
  varnames <- strsplit(gsub("\\[|\\]", "", 
                node_list$toString()), 
                split=", ")[[1]]
  edgelist <- c()
  for (i in 2:numNodes){
    for (j in 1:(i-1)){
      if (ugmat[j,i]==1) edgelist <- c(edgelist, 
                                        paste(varnames[j], 
                                        "---", 
                                        varnames[i]))
    }
  }
  
  varstring <- paste(varnames, collapse=" ")
  edgestring <- paste(1:length(edgelist),". ", 
                edgelist, "\n",sep="", collapse="")
  graphstring <- paste("\nGraph Nodes:\n", varstring, 
                        " \n\nGraph Edges: \n", 
                        edgestring, "\n", sep="")
  
  graphfilename <- "impossibly_long_graph_file_name_temporary.txt"
  if ("!"(file.exists(graphfilename))){
    write(graphstring, graphfilename)
    graphfile <- .jnew("java/io/File", graphfilename)
    newug_tetrad <- .jcall("edu/cmu/tetrad/graph/GraphUtils", 
                           "Ledu/cmu/tetrad/graph/Graph;", 
                           "loadGraphTxt", graphfile)
    newug_tetrad <- .jcast(newug_tetrad, "edu/cmu/tetrad/graph/Graph", 
                            check=TRUE)
    rm(graphfile)
    file.remove(graphfilename)
    return(newug_tetrad)
  } else {
    print("Whoops, don't want to overwrite existing file!")
    stop()
  }
}
#######################################3

