

library(ggplot2)
library(ape)


DataCham <- read_csv("data/chameleons/largest_chameleons.csv")
ChamTree<-read.nexus("data/chameleons/Tolley2013trimmedL.trees")

DataCham <- DataCham[match(ChamTree$tip.label, DataCham$Species), ]

ggplot(DataCham, aes(x = logbsl, y = logpter, color = as.factor(Regime))) +
  geom_point() +
  geom_smooth(method = "lm")


discrete_trait <- DataCham$Regime
ans <- ape::ace(discrete_trait, ChamTree, type = "discrete")

internal_regimes <- colnames(ans$lik.anc)[apply(ans$lik.anc, 1, which.max)]
regimes <- c(DataCham$Regime, internal_regimes)

edge_regimes <- factor(regimes[ChamTree$edge[,2]])
cls <- c("orange", "black")

plot(ChamTree, edge.color = cls[edge_regimes], edge.width = 3, cex = 1.3)


