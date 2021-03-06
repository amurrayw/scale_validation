library(ggplot2)
library(reshape)
library(plyr)
library(dplyr)

source("simulation_fofc_vs_fa.R")


#nScree(create.dataset.from.graph(generate.measurement.model(n.latent=3, n.measures.per.latent=9, n.latent.latent.edges=3),n=1000))


#(do.call(replicate(data.frame(unlist(nScree(create.dataset.from.graph(generate.measurement.model(n.latent=3, n.measures.per.latent=9, n.latent.latent.edges=3),n=1000))$Components)), n=3), what=rbind))









test.efa <- function(n.latent=3, n.measures.per.latent=4, n.latent.latent.edges=3, sample.size = 1000,n.rep=500){

    test.types <- c("noc", "naf", "nparallel", "nkaiser")

    simulations <- replicate(data.frame(unlist(nScree(create.dataset.from.graph(generate.measurement.model(n.latent=n.latent, n.measures.per.latent=n.measures.per.latent, n.latent.latent.edges=n.latent.latent.edges),n=sample.size))$Components)), n=n.rep)

   # print(simulations)

    simulations <- do.call(simulations, what=rbind)

    row.names(simulations) <- c()

    simulations <- data.frame(simulations)

    names(simulations) <- test.types

#    print(simulations)

    return(simulations)
}

##test.efa(n.rep=2)

##Produces barplot of mean n.latents each EFA method predicts under circumstances (and increasing multiples of n.measures.per.latent up to a maximum of 4 times).
efa.boxplots <- function(n.latent=3, n.measures.per.latent=4, n.latent.latent.edges=3, sample.size = 1000,n.rep=500){

    simulation.result <- list()
    n.measures <- c()
    for(i in 1:4){
        n.measures <- c(n.measures, n.measures.per.latent*i)
        simulation.result[[i]] <- test.efa(n.latent=n.latent, n.measures.per.latent=n.measures.per.latent*i, n.latent.latent.edges=n.latent.latent.edges, sample.size = sample.size, n.rep=n.rep)
    }

    df <- bind_rows(simulation.result, .id=c("ID"))
    
    for(i in sort(1:length(n.measures), decreasing=TRUE)){
            df$ID[df$ID==i] <- n.measures[i]
    }

   
    df.long<-melt(df)

   # df.long <- df.long[base::order(as.numeric(df.long$ID), decreasing=FALSE),]

print(ggplot(df.long) +
  geom_boxplot(aes(ID, value, colour=as.factor(variable)), show.legend=FALSE) +
  facet_grid(variable ~ .) + labs(title = paste("EFA: Predicted number of latents given", n.latent.latent.edges,  "latent edges and", n.latent, "latents" , sep=" ", collapse=""), x="Number of measures per latent", y="Number of predicted latents") + scale_x_discrete(limits = as.character(unique(as.numeric(df.long$ID)))))


}

##efa.boxplots(n.rep=2)


pdf("EFA_n_latents_boxplots.pdf")

for(i in c(3, 2, 1, 0)){
    efa.boxplots(n.latent.latent.edges=i, n.rep=500)
}
dev.off()





## Based on these plots, increasing the number of measures improves EFA significantly. Sample size doesn't really seem to change results (just increases variablity a bit) -- Based on simulations with sample size set to 250 instead of 1000, which are excluded as are otherwise not very informative given number of plots produced. At 12 measures per latent, EFA works very well (almost always correct), at 8 still makes mostly errors when latents are correlated (3 latent-latent edges and 3 latent variable). As number of latent-latent edges decreases, corelation between latents will decrease (more likely to end up with latents indep. of one another). Fewer measures are needed per latent when there are fewer latent-latent edges.



test.fofc.n.latents <- function(n.latent=3, n.measures.per.latent=4, n.latent.latent.edges=3, sample.size = 1000,n.rep=500, alpha=.01){

    fofc.model <- replicate(build.fofc.model(create.dataset.from.graph(generate.measurement.model(n.latent=n.latent, n.measures.per.latent=n.measures.per.latent, n.latent.latent.edges=n.latent.latent.edges),n=sample.size), alpha=alpha), n=n.rep)
    
    return(unlist(lapply(fofc.model, function(graph){return(length(grep(nodes(graph), pattern="L")))})))
}

fofc.boxplots <- function(n.latent=3, n.measures.per.latent=4, n.latent.latent.edges=3, sample.size = 1000,n.rep=500, alpha=1/sample.size){

    simulation.result <- list()
    n.measures <- c()
    for(i in 1:4){
        n.measures <- c(n.measures, n.measures.per.latent*i)
         results.df <- data.frame(test.fofc.n.latents(n.latent=n.latent, n.measures.per.latent=n.measures.per.latent, n.latent.latent.edges=n.latent.latent.edges ,sample.size=sample.size, n.rep=n.rep, alpha=alpha))
        names(results.df) <- (i)

        simulation.result[[i]] <- results.df
    }

  ##  print(simulation.result)

    df <- bind_rows(simulation.result, .id=c("ID"))
  ##  print(df)
        for(i in sort(1:length(n.measures), decreasing=TRUE)){
            df$ID[df$ID==i] <- n.measures[i]
    }

    df.long<-melt(df)

    return(ggplot(df.long, aes(x=ID, y=value),environment = localenv) + geom_boxplot() +labs(title = paste("FOFC: Predicted number of latents given", n.latent.latent.edges,  "latent edges and", n.latent, "latents" , sep=" ", collapse=""), x="Number of measures per latent", y="Number of predicted latents") + scale_x_discrete(limits = as.character(unique(as.numeric(df.long$ID)))))

}

##fofc.boxplots(n.rep=2)




fofc.box.plots.results <- list()
for(i in c(3, 2, 1, 0)){
    fofc.box.plots.results[[i+1]]<- fofc.boxplots(n.latent.latent.edges=i, n.rep=500)
}

localenv <- environment() 

pdf("FOFC_n_latents_boxplots_samp_1000.pdf")

do.call(grid.arrange, c(fofc.box.plots.results[1:2], fofc.box.plots.results[3:4], list(ncol=2, nrow=2)))

dev.off()



fofc.box.plots.results <- list()
for(i in c(3, 2, 1, 0)){
    fofc.box.plots.results[[i+1]]<- fofc.boxplots(n.latent.latent.edges=i, n.rep=500, sample.size=500)
}

localenv <- environment() 

pdf("FOFC_n_latents_boxplots_samp_500.pdf")

do.call(grid.arrange, c(fofc.box.plots.results[1:2], fofc.box.plots.results[3:4], list(ncol=2, nrow=2)))

dev.off()



fofc.box.plots.results <- list()
for(i in c(3, 2, 1, 0)){
    fofc.box.plots.results[[i+1]]<- fofc.boxplots(n.latent.latent.edges=i, n.rep=500, sample.size=250)
}

localenv <- environment() 

pdf("FOFC_n_latents_boxplots_samp_250.pdf")

do.call(grid.arrange, c(fofc.box.plots.results[1:2], fofc.box.plots.results[3:4], list(ncol=2, nrow=2)))

dev.off()

## Performance of FOFC goes down as sample size decreases. Also seems to get worse as number of latent-latent edges increases. Overall performance seems pretty similar to EFA in these cases with ~8-12 measures (though it does much better in the 4-8 measure, 2-3 latent-latent edge case). 

## EFA methods are: noc naf nparallel nkaiser.










#options(java.parameters = "-Xmx8000m")



## Segfault from generate.measurement.model when n.latents=10.
## Segfault from factanal when n.latents is 5 or larger.
compare.fofc.fa <- function(n.latents=1, n.measures.per.latent=4, n.impurities=list(n.sharing.latent=0, n.output.output.edges=0), n.latent.latent.edges=sample(1:choose(n.latents, 2), size=1), unif.min=1, unif.max=1, cut.off=.3, sample.size=1000, rotation="promax"){
#print(n.latents)
#print(n.latent.latent.edges)
#print(0)
    orig.graph <- generate.measurement.model(n.latents=n.latents, n.measures.per.latent=n.measures.per.latent, n.impurities=n.impurities, n.latent.latent.edges=n.latent.latent.edges, unif.min=unif.min, unif.max=unif.max)
#print(1)
    dataset <- create.dataset.from.graph(orig.graph, n=sample.size)
#print(2)
    fofc.graph <- build.fofc.model(data=dataset)
#print(3)
    fofc.score <-  score.fofc(fofc.graph, true.graph=as(orig.graph, "matrix"))
#print(4)
#write.csv(dataset, file="debug/crashing_dataset.csv")
    fa.graph <- factanal(x=dataset, factors=n.latents, rotation=rotation)
#print(5)
    fa.score <-  score.fa(fa.graph, true.graph=as(orig.graph, "matrix"), cut.off=cut.off)
#print(6)
    return(list(fofc.score = fofc.score, fa.score=fa.score))
}



##compare.fofc.fa(n.latents=sample(1:4, size=1))


pdf("edge_adj_results.pdf")

##TODO: Turn this bit into a function. 

promax.run <- replicate(compare.fofc.fa(n.latents=sample(1:4, size=1), unif.min=.01, rotation="promax", sample.size=500), n=100000)


save.image()

toplot<-(data.frame(do.call(promax.run[2,], what=rbind)))
##reshape2::melt(toplot, id.vars = NULL)
toplot<-reshape2::melt(toplot, id.vars = NULL)
print(ggplot(toplot, aes(x=variable, y=value)) + geom_violin()+labs(title = "FA: 500 obs, promax rotation, 100000 replicates.", x="", y="Rate"))

##+geom_jitter(height = 0, width = 0.1)

toplot<-(data.frame(do.call(promax.run[1,], what=rbind)))
##reshape2::melt(toplot, id.vars = NULL)

##Removing cases where fofc got number of latents wrong and as a result, interp becomes difficult (and graph difficult.

toplot <- toplot[-which(toplot[,c(1,3)]==0),] 
##toplot <- toplot[-which(toplot[,2]==1),] 
toplot<-reshape2::melt(toplot, id.vars = NULL)


print(ggplot(toplot, aes(x=variable, y=value)) + geom_violin()+labs(title = "FOFC: 500 obs, alpha=.01, 100000 replicates.", x="", y="Rate"))

##+geom_jitter(height = 0, width = 0.1)

##TODO: Also try version with varimax rotation instead of promax

##varimax

varimax.run <- replicate(compare.fofc.fa(n.latents=sample(1:4, size=1), unif.min=.01, rotation="varimax", sample.size=500), n=1000)


save.image()

toplot<-(data.frame(do.call(varimax.run[2,], what=rbind)))
##reshape2::melt(toplot, id.vars = NULL)
toplot<-reshape2::melt(toplot, id.vars = NULL)
print(ggplot(toplot, aes(x=variable, y=value)) + geom_violin()+ geom_violin()+labs(title = "FA: 500 obs, varimax rotation, 1000 replicates.", x="", y="Rate"))

dev.off()
##INTERP:

   ##  tpr: True Positive Rate: Number of correctly found edges (in
  ##        estimated graph) divided by number of true edges (in true
 ##         graph)
##
  ##   fpr: False Positive Rate: Number of incorrectly found edges
 ##         divided by number of true gaps (in true graph)
##
  ##   tdr: True Discovery Rate: Number of correctly found edges divided
 ##         by number of found edges (both in estimated graph)
##



table(replicate(as.numeric(names(which.max(table(unlist(nScree(x=create.dataset.from.graph(generate.measurement.model(n.latent=3, n.measures.per.latent=9, n.latent.latent.edges=3), n=1000), model="factors")$Components))))), n=1000))


#score.fa(factanal(x=create.dataset.from.file(), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"))


##table(replicate(unlist(nScree(x=create.dataset.from.file(n=1000), model="factors")$Components), n=1000))


table(replicate(as.numeric(names(which.max(table(unlist(nScree(x=create.dataset.from.file(n=1000), model="factors")$Components))))), n=1000))


table(replicate(as.numeric(names(which.max(table(unlist(nScree(x=create.dataset.from.file(n=1000, file="graph1.r.txt"), model="factors")$Components))))), n=1000))

par(mfrow=c(2,2))
hist(replicate(score.fa(factanal(x=create.dataset.from.file(n=1000, file="graph1.r.txt"), factors=3), true.graph=as(read.dag("graph1.r.txt"), "matrix"), cut.off=.3), n=1000)[3,])


hist(replicate(score.fofc(build.fofc.model(create.dataset.from.file(n=1000, file="graph1.r.txt")), true.graph=as(read.dag("graph1.r.txt"), "matrix")), n=1000)[3,])



##TODO: Need to also account for getting rotation correct (promax is supposed to be used when latents are "correlated" with one another). Should test with impure measures (i.e., correlation is through observed variables, not latents).
 plot(igraph.to.graphNEL(graph.adjacency(prune.fa.paths(factanal(x=create.dataset.from.file(n=1000, file="graph1.r.txt"), factors=3, rotation="promax"), .3))))

score.fa(factanal(x=create.dataset.from.file(n=1000, file="graph1.r.txt"), factors=3, rotation="promax"), true.graph=as(read.dag("graph1.r.txt"), "matrix"), cut.off=.3) 



hist(replicate(score.fa(factanal(x=create.dataset.from.file(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.3), n=1000)[3,])
hist(replicate(score.fa(factanal(x=create.dataset.from.file(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.5), n=1000)[3,])
hist(replicate(score.fa(factanal(x=create.dataset.from.file(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.7), n=1000)[3,])





##handle.null.fofc <- function(fofc.sim.object){
##	if(class=="matrix"){return(fofc.sim.object)}
##	else{
##		fofc.sim.object <- list.clean(fofc.sim.object, is.null)	
##	}
##}


#hist(replicate(score.fa(factanal(x=create.dataset.from.file(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.3), n=1000)[3,])


#set.seed(100);prune.fa.paths(factanal(x=create.dataset.from.file(n=1000), factors=2), cut.off=.3)[, c(10,9)]




















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





