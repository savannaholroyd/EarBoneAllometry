#Testing for stabilizing selection around optimum in therapsids
library(ape)
source("http://www.graemetlloyd.com/pubdata/functions_2.r")
library(phytools)
library(readr)
library(tidyverse)
library(geiger)
library(OUwie)
library(slouch)
devtools::install_github("wabarr/ggphylomorpho")
setwd("D:/Work backup/Research Projects/RL morpho/Analyses")


DataSynOUR <- as.data.frame(read.csv("DataSynOUready.csv")) #residuals from PGLS regressions
#this dataset doesn't include cynodonts because I think they'd make the model too complicated

#first I need to make a time calibrated tree for this dataset

AllomeTree<-read.nexus("RLallometry.nex") #load the cladogram. This is a tree that I compiled based on several datasets
#plotTree(AllomeTree,ftype="i",fsize=0.6,lwd=1)
rownames(DataSynOUR) <- DataSynOUR$Species
obj<-name.check(AllomeTree,DataSynOUR)
AllomeTree<-drop.tip(AllomeTree, obj$tree_not_data)
agesA<-read.table("SynapsidTreeAgesA.txt",row.names=1) #needs to be text file with names as in tree and age of first occurrence

#calculating branch lengths
ttreeA<-date.phylo(AllomeTree, agesA, rlen=1, method="equal") #uses method from Brusatte et al., 2008
plotTree(ttreeA,ftype="i",fsize=0.6,lwd=1)

# Reorder dataframe so the species labels match the tree tip labels
DataSynOUR <- DataSynOUR[match(ttreeA$tip.label, DataSynOUR$Species), ]

#############################

#I don't feel like making five million files, so just manually copy and paste the hypothesis column you want to use from
#the "RL_allometry_forR_finalOU" file to the "DataSynOUready" excel spreadsheet and re-run code for each hypothesis
#note that the null hypothesis is more or less represented by the OUM models. Same parameters except optimum can be
#different, but we'd expect the optimum to basically be the same (0) for all of these anyway since the data are residuals
DataSynOUR <- as.data.frame(read.csv("DataSynOUready.csv")) #reloading the data file once the hypotheses have been altered as desired
DataSynOUR <- DataSynOUR[match(ttreeA$tip.label, DataSynOUR$Species), ]

#modeling evolution towards optimum with different amounts of drift and adaptive rates
#THIS ALL REQUIRES MOST RECENT VERSION OF R

#map hearing ability as discrete character on tree
earS <- as.numeric(DataSynOUR$hypothesis)
names(earS) <-ttreeA$tip.label
mod = matrix(c(0,0,1,0),2) #this forces is to only allow 0->1 transformations since right now it thinks ear is basal
#I need to not use the above line when doing hypothesis 2 because I need to test for a loss in bidentalians
ttreeAR<-make.simmap(ttreeA,earS,nsim=1, model = mod) #add for when it reconstructs the ear as basal
#check the tree for troubleshooting
#do this before running OU analysis because sometimes it reconstructs the character states weirdly
cols<-setNames(c("blue","red"),c(0,1))
plotSimmap(ttreeAR,cols,ftype="i",fsize=0.7)

fitOUMV<-OUwie(ttreeAR, DataSynOUR,model="OUMV",simmap.tree=TRUE) #uses hypothesized ear evolution to fit an OU model
#of coefficient evolution. OUMV means the model is using a different drift for each selective regime (ear or no)
#this gives an AIC (and sample size corrected AICc) that can be compared with other models

#compare to OU model with same drift and adaptive rate for each selective regime
#this is the only one that can be reported for hypothesis0 since there's only one regime
fitOUM<-OUwie(ttreeAR, DataSynOUR,model="OUM",simmap.tree=TRUE) 

#compare to OU model with same drift but different adaptive rate for each selective regime
fitOUMA<-OUwie(ttreeAR, DataSynOUR,model="OUMA",simmap.tree=TRUE) 

#compare to OU model with different drift and adaptive rate for each selective regime
fitOUMVA<-OUwie(ttreeAR, DataSynOUR,model="OUMVA",simmap.tree=TRUE) 

#see results
fitOUMV
fitOUM
fitOUMA
fitOUMVA

#############################################################################################

#I want to also try running this with a separate PGLS regression for each group because I think there's a phylogenetic impact
#making it look like each clade has strong directional selection when they're actually just moving to their own optima

library(OUwie)
DataSynOURS <- as.data.frame(read.csv("DataSynOUR_sub.csv"))

#map hearing ability as discrete character on tree
earS <- as.numeric(DataSynOURS$hypothesis)
names(earS) <-ttreeA$tip.label
mod = matrix(c(0,0,1,0),2) #this forces is to only allow 0->1 transformations since right now it thinks ear is basal
#I need to not use the above line when doing hypothesis 2 because I need to show a loss in bidentalians
ttreeAR<-make.simmap(ttreeA,earS,nsim=1, model = mod) #add for when distribution of regimes could make it think the ear is basal
#check the tree for troubleshooting
#do this before running OU analysis because sometimes it reconstructs a weird tree
cols<-setNames(c("blue","red"),c(0,1))
plotSimmap(ttreeAR,cols,ftype="i",fsize=0.7)

fitOUMV<-OUwie(ttreeAR, DataSynOURS,model="OUMV",simmap.tree=TRUE) #uses hypothesized ear evolution to fit an OU model
#of coefficient evolution. OUMV means the model is using a different drift for each selective regime (ear or no)
#this gives an AIC (and sample size corrected AICc) that can be compared with other models

#compare to OU model with same drift and adaptive rate for each selective regime
#this is the only one that can be reported for hypothesis0 since there's only one regime
fitOUM<-OUwie(ttreeAR, DataSynOURS,model="OUM",simmap.tree=TRUE) 

#compare to OU model with same drift but different adaptive rate for each selective regime
fitOUMA<-OUwie(ttreeAR, DataSynOURS,model="OUMA",simmap.tree=TRUE) 

#compare to OU model with different drift and adaptive rate for each selective regime
fitOUMVA<-OUwie(ttreeAR, DataSynOURS,model="OUMVA",simmap.tree=TRUE) 

#see results
fitOUMV
fitOUM
fitOUMA
fitOUMVA


####################################################

##CHAMELEONS

setwd("D:/Work backup/Research Projects/Chameleons/Allometry/Data/Chameleon Measurements")

#loading data
coeff <- na.omit(read_csv("coefficientsOU.csv")) #this is using the coefficients from intraspecific allometric regressions
#Can also try using residuals from interspecific allometry so it's more comparable with the synapsid dataset
SpNames<- as.list(t(coeff[,1]))
Tolley2013Trees <- "Tolley2013trimmed.trees"
cham.tree<-read.nexus(Tolley2013Trees)
#To gaze upon its beauty: 
#plotTree(cham.tree,ftype="i",fsize=0.6,lwd=1)

#ancestral state reconstruction
ear <- c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
names(ear) <-SpNames
mod = matrix(c(0,0,1,0),2) #this forces is to only allow 0->1 transformations since right now it thinks ear is basal
pterMap<-make.simmap(cham.tree,ear,model=mod,nsim=1) #will want to change ER if you want a different evolutionary model
plot(pterMap)

#BM model
fitBM <- OUwie(pterMap, chamOU,model="BMS",simmap.tree=TRUE, diagn = T)#I might use this BM model instead because I'm having trouble getting theta
fitBM
fitBM$eigval
fitBM$solution.se

#OU models
chamOU <- as.data.frame(read.csv('coefficientsOU.csv'))

fitOUMV<-OUwie(pterMap, chamOU,model="OUMV",simmap.tree=TRUE, diagn=T) #uses reconstructed ear evolution to fit an OU model
#of coefficient evolution. OUMV means the model is using a different drift for each selective regime (ear or no)
#this gives an AIC (and sample size corrected AICc) that can be compared with other models

#compare to OU model with same drift and adaptive rate for each selective regime
#fitOUM<-OUwie(pterMap, chamOU,model="OUM",simmap.tree=TRUE, diagn=T) 

#compare to OU model with same drift but different adaptive rate for each selective regime
#fitOUMA<-OUwie(pterMap, chamOU,model="OUMA",simmap.tree=TRUE, diagn=T) 

#compare to OU model with different drift and adaptive rate for each selective regime
fitOUMVA<-OUwie(pterMap, chamOU,model="OUMVA",simmap.tree=TRUE, diagn=T) 

fitOUMV
fitOUM
fitOUMA
fitOUMVA
