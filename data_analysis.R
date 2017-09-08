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
library(EFAutilities)
data(BFI228)
BFI228 <- apply(BFI228, 2, as.numeric)

plot(build.fofc.model(BFI228, alpha=.01))

##Associates all 8 with the latent they're supposed to go with. Idea: Only run FOFC on set of variables hypothesized to be related already. May also need to discuss permuting columns (due to a path dependency)?
plot(build.fofc.model(BFI228[,1:8], alpha=.000001)) 

## Following the rule Joe and Erich recommend does not do so (4 remain unclustered).
plot(build.fofc.model(BFI228[,1:8], alpha=1/nrow(BFI228)))

data(CPAI537)
CPAI537 <- apply(CPAI537, 2, as.numeric)
plot(build.fofc.model(CPAI537, alpha=.01))

###
library(OpenMx)
data("HS.ability.data")
#sapply(HS.ability.data[complete.cases(HS.ability.data),-c(1:6)], as.numeric)

plot(build.fofc.model(sapply(HS.ability.data[complete.cases(HS.ability.data),-c(1:6)], as.numeric), alpha=.01)) ###This one looks like it should be a two factor (bifactor model). FTFC would be more appropriate. Can't tell if the orig. paper came to the same conclusion, but the title of it suggests that it might have: "Holzinger, K., and Swineford, F. (1939). A study in factor analysis: The stability of a bifactor solution"







