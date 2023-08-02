#Testing for different allometric intercepts among therapsids
#Note that this method assumes identical slope across regimes, at least the way I'm running it
#see if there's a way to run it with different slopes between regimes.
#Also can't get at differences in adaptive rate or evolutionary step size with this setup

#Re-run the analysis with different distributions of the selective regimes since I'm not sure which groups
#were exposed to a hearing selective regime
#the fit of the model might give a hint at which groups were under which selective regimes

library(ape)
source("http://www.graemetlloyd.com/pubdata/functions_2.r")
library(phytools)
library(readr)
library(tidyverse)
library(geiger)
library(slouch)
#remotes::install_github("wabarr/ggphylomorpho")
library(ggphylomorpho)
#setwd("D:/Work backup/Research Projects/RL morpho/Analyses")

############################

DataSynOU <- read_csv("data/RL_allometry_forR_finalOU.csv") #this dataset doesn't include cynodonts because
#I think they'd make the model too complicated. I'd have to add more regimes to accommodate them and would run into overfitting

#first I need to make a time calibrated tree for this dataset
AllomeTree<-read.nexus("data/RLallometry.nex")
#plotTree(AllomeTree,ftype="i",fsize=0.6,lwd=1) #To gaze upon its beauty

obj<-name.check(AllomeTree,DataSynOU) #check that the tree matches the dataset
AllomeTree<-drop.tip(AllomeTree, obj$tree_not_data)
agesA<-read.table("data/SynapsidTreeAgesA.txt",row.names=1) #needs to be text file with names as in tree
#and age of first occurrence

#calculating branch lengths
ttreeA<-date.phylo(AllomeTree, agesA, rlen=1, method="equal") #uses method from Brusatte et al., 2008
#plotTree(ttreeA,ftype="i",fsize=0.6,lwd=1) #to gaze upon its beauty

## Reorder dataframe so the species labels match the tree tip labels
DataSynOU <- DataSynOU[match(ttreeA$tip.label, DataSynOU$Species), ]
#DataSynOU$Species == ttreeA$tip.label

###########################

## HYPOTHESIS 0 none hearing
## Null hypothesis that all synapsids were in the same selective regime
## Fit the model without a fixed factor for the selective regime because they're all the same regime
tree_height <- max(node.depth.edgelength(ttreeA))

m0 <- slouch.fit(phy = ttreeA,
                 hillclimb = F,
                 vy_values = seq(0.001, 0.15, length.out = 60),
                 hl_values = seq(0.1, 150, length.out = 60),
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logJaw), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m0)
summary(m0)

resid1 <- m0$beta_primary$residuals

## Save the residuals to a nexus file
l <- list()
for (i in 1:length(ttreeA$tip.label)){
  splab <- ttreeA$tip.label[i]
  #splab <- gsub(".", "", splab, fixed = TRUE)
  #splab <- gsub(".", "", splab, fixed = TRUE)
  l[[splab]] <- list(resid1[i])
}

ape::write.nexus.data(l, file = "data/logRL_logJaw_residuals.nex", format = "continuous")

## Save the regimes to a nexusfile
l2 <- list()
for (i in 1:length(ttreeA$tip.label)){
  splab <- ttreeA$tip.label[i]

  dat <- DataSynOU[DataSynOU$Species == splab,c("Regime0", "Regime1", "Regime2", "Regime3", "Regime4")]
  dat1 <- unlist(dat)
  dat2 <- unname(dat1)
  l2[splab] <- list(dat2)
}

ape::write.nexus.data(l2, file = "data/regimes.nex", format = "standard")

## save the tree to a file
ape::write.nexus(ttreeA, file = "data/therapsids_tree.nex")


brown.fit(ttreeA,
          species = DataSynOU$Species,
          response = resid1)


