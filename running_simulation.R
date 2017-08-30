
source("simulation_fofc_vs_fa.R")






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





##handle.null.fofc <- function(fofc.sim.object){
##	if(class=="matrix"){return(fofc.sim.object)}
##	else{
##		fofc.sim.object <- list.clean(fofc.sim.object, is.null)	
##	}
##}


#hist(replicate(score.fa(factanal(x=create.dataset(n=1000), factors=2), true.graph=as(read.dag("sim_graph_2_lat_pure_measure.r.txt"), "matrix"), cut.off=.3), n=1000)[3,])


#set.seed(100);prune.fa.paths(factanal(x=create.dataset(n=1000), factors=2), cut.off=.3)[, c(10,9)]

