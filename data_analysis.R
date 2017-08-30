library(foreign)

source("simulation_fofc_vs_fa.R")

####TEMP DATA GSS 2012
gss.df <- read.spss("GSS2012.sav", to.data.frame=TRUE)

temp <- gss.df[, grep(names(gss.df), pattern="^ab")]

abort.names <-names(temp)
#temp <- ifelse(temp=="YES", 1 ,0)
#temp <- temp[complete.cases(temp), ]

gss_all.df<-read.spss("~/Desktop/GSS7216_R1.sav", to.data.frame=TRUE)
temp <- gss_all.df[, abort.names]

temp <- ifelse(temp=="YES", 1 ,0)

temp <-data.frame(temp, gss_all.df$year)
names(temp)[8] <- "year"
temp <- temp[complete.cases(temp), ]

par(mfrow=c(5,5))
dlply(.data=temp, .variables="year", function(data){plot(build.fofc.model(data[,-8], alpha=1/nrow(data)))})


fofc(loadContinuousData(temp[,-8]), alpha=1/nrow(temp))


build.fofc.model(gss_all.df[,grep(names(gss_all.df), pattern="^grnte")][complete.cases(gss_all.df[,grep(names(gss_all.df), pattern="^grnte")]),], alpha=1/2067)





##ALSO look at Bollen dataset in piecwiseSEM
library(sem)
data(Bollen)
build.fofc.model(Bollen, alpha=1/nrow(Bollen))


data(CNES)
fofc(df=loadContinuousData(sapply(CNES, as.numeric)), alpha=.0000001) ## All clustered around one latent.
fofc(df=loadContinuousData(sapply(CNES, as.numeric)), alpha=.01) ## No clustering.
fofc(df=loadContinuousData(sapply(CNES, as.numeric)), alpha=.00001) ## ONly 3 var. clustered. 

data(Tests)

fofc(df=loadContinuousData(sapply(Tests, as.numeric)[complete.cases(Tests),]), alpha=.001) ##  Assigns all to one latent (but only has 16 observations!).


#####
