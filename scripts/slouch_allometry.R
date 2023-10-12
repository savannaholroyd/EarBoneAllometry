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
DataSynOU$Species[!DataSynOU$Species %in% AllomeTree$tip.label]
not_in_data <- AllomeTree$tip.label[!AllomeTree$tip.label %in% DataSynOU$Species]

AllomeTree<-drop.tip(AllomeTree, not_in_data)
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
                 hillclimb = T,
                 #vy_values = seq(0.001, 0.15, length.out = 30),
                 #hl_values = seq(0.1, 150, length.out = 30),
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw,
                 mv.response = rep(0.1 * var(DataSynOU$logRL), length(ttreeA$tip.label)), #plus 5.5% intraspecific variation
                 mv.direct.cov = rep(0.1 * var(DataSynOU$logJaw), length(ttreeA$tip.label))) #plus 5.5% intraspecific variation
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

hypotheses <- c("Hypothesis0", "Hypothesis1", "Hypothesis2", "Hypothesis3", "Hypothesis4")
for (hyp in hypotheses){
  l2 <- list()
  for (i in 1:length(ttreeA$tip.label)){
    splab <- ttreeA$tip.label[i]
    
    dat <- DataSynOU[DataSynOU$Species == splab,hyp]
    dat1 <- unlist(dat)
    dat2 <- unname(dat1)
    l2[splab] <- list(dat2)
  }
  ape::write.nexus.data(l2, file = paste0("data/regimes_", hyp, ".nex"), format = "standard")
}

DataSynOU$resid1 <- resid1

library(tidytree)
library(ggtree)
library(patchwork)

DataSynOU %>%
  ggplot(aes(x = logJaw, y = logRL, color = factor(Hypothesis2))) +
  geom_point()

tr1 <- as_tibble(ttreeA)
for (i in 1:4){
  ans <- ape::ace(DataSynOU[[paste0("Hypothesis", i)]], ttreeA, type = "discrete")
  internal_states <- c(0, 1)[apply(ans$lik.anc, 1, which.max)]
  tip_states <- DataSynOU[[paste0("Hypothesis", i)]]
  tr1[[paste0("regime",i)]] <- c(tip_states, internal_states)
}
tr2 <- as.treedata(tr1)

tree_plots <- list(
  ggtree(tr2, aes(color = factor(regime1))) + 
    geom_tiplab() + ggtitle("Hypothesis 1") + theme(legend.pos = c(0.1, 0.9)),
  ggtree(tr2, aes(color = factor(regime2))) + 
    geom_tiplab() + ggtitle("Hypothesis 2") + theme(legend.pos = "none"),
  ggtree(tr2, aes(color = factor(regime3))) + 
    geom_tiplab() + ggtitle("Hypothesis 3") + theme(legend.pos = "none"),
  ggtree(tr2, aes(color = factor(regime4))) + 
    geom_tiplab() + ggtitle("Hypothesis 4") + theme(legend.pos = "none")
)


p1 <- Reduce("+", tree_plots) +
  plot_layout(ncol = 1)

intercept <- m0$beta_primary$coefficients[1,1]
slope <- m0$beta_primary$coefficients[2,1]

allometry_plots <- list(
  ggplot(DataSynOU, (aes(x = logJaw, y = logRL, color = factor(Hypothesis1)))) +
    geom_point() + geom_abline(slope = slope, intercept = intercept),
    theme_classic() + theme(legend.pos = c(0.3, 0.8)),
  ggplot(DataSynOU, (aes(x = logJaw, y = logRL, color = factor(Hypothesis2)))) +
    geom_point() + geom_abline(slope = slope, intercept = intercept) +
    theme_classic() + theme(legend.pos = "none"),
  ggplot(DataSynOU, (aes(x = logJaw, y = logRL, color = factor(Hypothesis3)))) +
    geom_point() + geom_abline(slope = slope, intercept = intercept) +
    theme_classic() + theme(legend.pos = "none"),
  ggplot(DataSynOU, (aes(x = logJaw, y = logRL, color = factor(Hypothesis4)))) +
    geom_point() + geom_abline(slope = slope, intercept = intercept) +
    theme_classic() + theme(legend.pos = "none")
)

p2 <- Reduce("+", allometry_plots) +
  theme(legend.pos = "none") +
  plot_layout(ncol = 1)

p1 | p2

dfs <- list()
for (i in 1:4){
  df1 <- read.table(paste0("output/state_dependent_BM_character_", i, ".log"), header = TRUE)
  
  df1 <- df1 %>%
    select(-starts_with("state_branch")) %>%
    select(-starts_with("num_changes")) %>%
    select(-starts_with("branch_rates"))
  
  df1$zeta_diff <- df1$zeta.1. - df1$zeta.2.
  df1$character <- factor(i)
  dfs[[i]] <- df1
}
df2 <- bind_rows(dfs)

zeta_plots <- list()
for (i in 1:4){
  p <- dfs[[i]] %>%
    ggplot(aes(x = zeta_diff)) +
    #geom_density(alpha = 0.5, fill = "gray") +
    geom_histogram() +
    theme_classic() +
    geom_vline(xintercept = 0, color = "red") +
    xlim(c(-1.0, 2.0))
  zeta_plots[[i]] <- p
}

p3 <- Reduce("+", zeta_plots) +
  plot_layout(ncol = 1)

combined_plot <- p1 | p2 | p3

ggsave("figures/hypotheses_may_moore_model.pdf", combined_plot,
       width = 400, height = 400, units = "mm")


summary(lm(DataSynOU$resid1 ~ DataSynOU$Hypothesis1))


## save the tree to a file
ape::write.nexus(ttreeA, file = "data/therapsids_tree.nex")


brown.fit(ttreeA,
          species = DataSynOU$Species,
          response = resid1)


