#Testing for different allometric intercepts among chameleons

library(ape)
library(tidyverse)
library(slouch)

############################

DataCham <- read_csv("data/chameleons/chameleons_speciesdata.csv")

ChamTree<-read.nexus("data/chameleons/Tolley2013trimmedL.trees")

## Reorder dataframe so the species labels match the tree tip labels
DataCham <- DataCham[match(ChamTree$tip.label, DataCham$Species), ]

###########################

## Fit the model to get residuals
tree_height <- max(node.depth.edgelength(ChamTree))

m0 <- slouch.fit(phy = ChamTree,
                 hillclimb = T,
                 vy_values = seq(0.00001, 0.004, length.out = 30),
                 hl_values = seq(0.0001, 140, length.out = 30),
                 species = DataCham$Species,
                 response = DataCham$logpter_mean,
                 mv.response = DataCham$logpter_varmean,
                 direct.cov = DataCham$logbsl_mean,
                 mv.direct.cov = DataCham$logbsl_varmean
)
                 
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





