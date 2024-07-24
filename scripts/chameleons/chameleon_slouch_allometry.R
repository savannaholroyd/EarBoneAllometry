#Testing for different allometric intercepts among chameleons

library(ape)
library(phytools)
library(readr)
library(tidyverse)
library(geiger)
library(slouch)
library(ggphylomorpho)

############################

DataCham <- read_csv("data/chameleons/largest_chameleons.csv")


ChamTree<-read.nexus("data/chameleons/Tolley2013trimmedL.trees")
#ChamTree2<-read.nexus("data/Tolley2013trimmed.trees")
#plotTree(ChamTree,ftype="i",fsize=0.6,lwd=1) #To gaze upon its beauty

## Reorder dataframe so the species labels match the tree tip labels
DataCham <- DataCham[match(ChamTree$tip.label, DataCham$Species), ]
#DataCham$Species == ChamTree$tip.label

###########################

## Fit the model to get residuals
tree_height <- max(node.depth.edgelength(ChamTree))

m0 <- slouch.fit(phy = ChamTree,
                 hillclimb = T,
                 vy_values = seq(0.0001, 0.05, length.out = 30),
                 hl_values = seq(0.001, 100, length.out = 30),
                 species = DataCham$Species,
                 response = DataCham$logpter,
                 direct.cov = DataCham$logbsl,
                 mv.response = rep(0.02 * var(DataCham$logpter), length(ChamTree$tip.label)), #plus 5.5% intraspecific variation
                 mv.direct.cov = rep(0.02 * var(DataCham$logbsl), length(ChamTree$tip.label))) #plus 5.5% intraspecific variation
plot(m0)
summary(m0)

resid1 <- m0$beta_primary$residuals

## Save the residuals to a nexus file
l <- list()
for (i in 1:length(ChamTree$tip.label)){
  splab <- ChamTree$tip.label[i]
  #splab <- gsub(".", "", splab, fixed = TRUE)
  #splab <- gsub(".", "", splab, fixed = TRUE)
  l[[splab]] <- list(resid1[i])
}

ape::write.nexus.data(l, file = "data/chameleons/logpter_logbsl_residuals.nex", format = "continuous")

## Save the regimes to a nexusfile
l2 <- list()
for (i in 1:length(ChamTree$tip.label)){
  splab <- ChamTree$tip.label[i]
  
  dat <- DataCham[DataCham$Species == splab,hyp]
  dat1 <- unlist(dat)
  dat2 <- unname(dat1)
  l2[splab] <- list(dat2)
}
ape::write.nexus.data(l2, file = paste0("data/chameleons/regimes.nex"), format = "standard")





