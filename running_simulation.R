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

fofc.boxplots(n.rep=2)




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















compare.fofc.fa <- function(){


    fofc.score <-  score.fofc(build.fofc.model(create.dataset.from.file(n=1000, file="graph1.r.txt")), true.graph=as(read.dag("graph1.r.txt"), "matrix"))


    fa.score <- score.fa(factanal(x=create.dataset.from.file(n=1000, file="graph1.r.txt"), factors=3, rotation="promax"), true.graph=as(read.dag("graph1.r.txt"), "matrix"), cut.off=.3)

    return(list(fofc.score = fofc.score, fa.score=fa.score))



}



compare.fofc.fa()







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

