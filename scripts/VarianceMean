### Calculate variance of the mean for chameleons and synapsids from intraspecific data for use in slouch models
library(tidyverse)



setwd("D:/Work backup/Research Projects/Chameleons/Measurements")



####### CHAMELEONS


DataAll <- read.csv("masterAll.csv")

sps <- tibble()  #make a list of all species
for (i in 1:nrow(DataAll)) {
  if ((DataAll$Species[i] %in% sps) == FALSE) {
    sps <- append(sps, DataAll$Species[i])
  }
  
}

ns <- tibble()  #list of sample sizes for each species
for (k in 1:length(sps)) {
    ns <- append(ns, nrow(DataAll[DataAll$Species == sps[k],]))
  
}

###intraspecific variation pterygoid plate area

varspter <- tibble() #calculate variance of the mean for chameleon pter
for (j in 1:length(sps)) {
  jsp <- tibble()
  for (l in 1:nrow(DataAll)) { #make subset of dataset for each species
    if ((DataAll$Species[l]) == as.character(sps[j])) {
      jsp <- append(jsp, DataAll$logpter[l])
    }
    
  }
  av <- mean(as.numeric(jsp))
  diffs <- tibble()
  for (m in 1:length(jsp)) { #calculate squared differences between measurement and mean
    mdiff <- (as.numeric(unlist(jsp[m]))-av)^2
    diffs <- append(diffs, mdiff)
    
  }
  varsppter <- (sum(as.numeric(diffs))/(as.numeric(ns[j])-1))/as.numeric(ns[j]) #sigma squared/n
  varsppter <- (varsppter/av)*100 #leave this in if we need to plug in variance of the mean as a percentage of the mean
  varspter <- append(varspter, varsppter)
}

avevar <- mean(as.numeric(varspter)) #about 0.35% intraspecific variation for pterygoid plate

###measurement error pterygoid plate area

DataError <- read.csv("error measurement.csv")

vare <- tibble() #calculate variance of the mean for replicates
for (p in 1:nrow(DataError)) {
    averr <- mean(DataError$logpter1[p],DataError$logpter2[p], DataError$logpter3[p])
    sumsqdiff <- sum((DataError$logpter1[p]-averr)^2,(DataError$logpter2[p]-averr)^2,(DataError$logpter3[p]-averr)^2)
    varspe <- (sumsqdiff/2)/3 #sigma squared/n, note that n is 3 for measurement replicates
    
    varspe <- (varspe/averr)*100 #leave this in if we need to plug in variance of the mean as a percentage of the mean
    vare<- append(vare, varspe)
}

avevarerr <- mean(as.numeric(vare)) #about 0.025% measurement error




####### SYNAPSIDS

setwd("D:/Work backup/Research Projects/Chameleons/Allometry/Data/Analyses")

###intraspecific variation reflected lamina area
DataD <- read.csv("DiictodonError.csv")

avD <- mean(as.numeric(DataD$logRL))
Ddiffs <- tibble()
for (d in 1:nrow(DataD)) { #calculate squared differences between measurement and mean
    Ddiff <- (as.numeric(unlist(DataD$logRL[d]))-avD)^2
    Ddiffs <- append(Ddiffs, Ddiff)
    
}
varsD <- (sum(as.numeric(Ddiffs))/(nrow(DataD)-1))/nrow(DataD) #sigma squared/n
varsD <- (varsD/avD)*100 #leave this in if we need to plug in variance of the mean as a percentage of the mean
#about 14% intraspecific variation based on Diictodon


###measurement error reflected lamina area

DataErrorS <- read.csv("SynapsidError.csv")

vares <- tibble() #calculate variance of the mean for replicates
for (s in 1:nrow(DataErrorS)) {
  averrs <- (DataErrorS$logRL1[s]+DataErrorS$logRL2[s])/2
  sumsqdiffs <- sum((DataErrorS$logRL1[s]-averrs)^2,(DataErrorS$logRL2[s]-averrs)^2)
  varspes <- (sumsqdiffs/2) #sigma squared/n, note that n is 2 for measurement replicates
  
  varspes <- (varspes/averrs)*100 #leave this in if we need to plug in variance of the mean as a percentage of the mean
  vares<- append(vares, varspes)
}

avevarerrs <- mean(as.numeric(vares)) #about 0.028% measurement error
