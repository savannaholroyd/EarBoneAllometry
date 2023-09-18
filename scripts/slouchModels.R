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
devtools::install_github("wabarr/ggphylomorpho")
setwd("D:/Work backup/Research Projects/RL morpho/Analyses")

############################

DataSynOU <- read_csv("RL_allometry_forR_finalOU.csv") #this dataset doesn't include cynodonts because
#I think they'd make the model too complicated. I'd have to add more regimes to accommodate them and would run into overfitting

#first I need to make a time calibrated tree for this dataset
AllomeTree<-read.nexus("RLallometry.nex")
#plotTree(AllomeTree,ftype="i",fsize=0.6,lwd=1) #To gaze upon its beauty 
rownames(DataSynOU) <- DataSynOU$Species
obj<-name.check(AllomeTree,DataSynOU) #check that the tree matches the dataset
AllomeTree<-drop.tip(AllomeTree, obj$tree_not_data)
agesA<-read.table("SynapsidTreeAgesA.txt",row.names=1) #needs to be text file with names as in tree
#and age of first occurrence

#calculating branch lengths
ttreeA<-date.phylo(AllomeTree, agesA, rlen=1, method="equal") #uses method from Brusatte et al., 2008
#plotTree(ttreeA,ftype="i",fsize=0.6,lwd=1) #to gaze upon its beauty

## Reorder dataframe so the species labels match the tree tip labels
DataSynOU <- DataSynOU[match(ttreeA$tip.label, DataSynOU$Species), ]
#DataSynOU$Species == ttreeA$tip.label

###########################

##HYPOTHESIS 0 none hearing
#Null hypothesis that all synapsids were in the same selective regime
## Fit the model without a fixed factor for the selective regime because they're all the same regime

m0 <- slouch.fit(phy = ttreeA,
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60), 
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m0)
summary(m0)

#calculated intraspecific variation by taking variance around the mean for Diictodon, Lystrosaurus, and chameleon samples
#all gave variance under 0.5%, so this should be a reasonable input
#btw I think variance is so small because this is the log transformed data
#and that's the data we're feeding to the model, so the variance should be on this same scale

#rinse and repeat for every hypothesized distribution of selective regimes

######

##HYPOTHESIS 1: All therapsids hearing

#ancestral state reconstruction to add regimes onto the tree
## Reconstruct maximum-likelihood ancestral states, using equal-rates Markov chain model
anc_recon1 <- ace(DataSynOU$hypothesis1, phy = ttreeA, type = "discrete", model = "ER")
MLE_idx1 <- apply(anc_recon1$lik.anc, 1, which.max)
MLE_hypothesis1 <- colnames(anc_recon1$lik.anc)[MLE_idx1]
ttreeA1 <- ttreeA #create a new phylogeny object that I can put regime tip labels onto
ttreeA1$node.label <- MLE_hypothesis1 #asign node labels based on anc recon
plot(ttreeA1, show.node.label = TRUE) #check that it worked

## Fit the model
m1 <- slouch.fit(phy = ttreeA1,
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60),
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw, 
                 fixed.fact = DataSynOU$hypothesis1,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation

plot(m1)
summary(m1) #slightly better than null

######

##HYPOTHESIS 2: Only basal anomodonts and therocephalians hearing
anc_recon2 <- ace(DataSynOU$hypothesis2, phy = ttreeA, type = "discrete", model = "ER")
MLE_idx2 <- apply(anc_recon2$lik.anc, 1, which.max)
MLE_hypothesis2 <- colnames(anc_recon2$lik.anc)[MLE_idx2]

## Create a new phylogeny object
ttreeA2 <- ttreeA
ttreeA2$node.label <- MLE_hypothesis2

## Fit the model
m2 <- slouch.fit(phy = ttreeA2,
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60), 
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw, 
                 fixed.fact = DataSynOU$hypothesis2,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m2)
summary(m2) #slightly worse than the null

#####

##HYPOTHESIS 3: Only basal anomodonts hearing
anc_recon3 <- ace(DataSynOU$hypothesis3, phy = ttreeA, type = "discrete", model = "ER")
MLE_idx3 <- apply(anc_recon3$lik.anc, 1, which.max)
MLE_hypothesis3 <- colnames(anc_recon3$lik.anc)[MLE_idx3]

## Create a new phylogeny object
ttreeA3 <- ttreeA
ttreeA3$node.label <- MLE_hypothesis3

## Fit the model
m3 <- slouch.fit(phy = ttreeA3,
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60), 
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw, 
                 fixed.fact = DataSynOU$hypothesis3,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m3)
summary(m3) #just barely better than H1

#####

##HYPOTHESIS 4: All anomodonts and theorcephalians hearing 
anc_recon4 <- ace(DataSynOU$hypothesis4, phy = ttreeA, type = "discrete", model = "ER")
MLE_idx4 <- apply(anc_recon4$lik.anc, 1, which.max)
MLE_hypothesis4 <- colnames(anc_recon4$lik.anc)[MLE_idx4]

## Create a new phylogeny object
ttreeA4 <- ttreeA
ttreeA4$node.label <- MLE_hypothesis1

## Fit the model
m4 <- slouch.fit(phy = ttreeA4,
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60), 
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw, 
                 fixed.fact = DataSynOU$hypothesis4,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m4)
summary(m4) #performs about as well as the null

summary(m0) #none hearing
summary(m1) #all therapsids hearing
summary(m2) #only therocephalians and basal anomodonts hearing
summary(m3) #only basal anomodonts hearing
summary(m4) #only therocephalians and all anomodonts hearing

#overall models perform similarly, probably because sample size is so low
#m1 and m3 perform slightly better than other models
#every model finds higher intercept for regime 1 than for regime 0
#probably driven by the basal anomodonts

######

##HYPOTHESIS 5: Only therocephalians and basal anomodonts hearing (hypothesis2), but also shift at therapsids (hypothesis1)
anc_recon5 <- ace(DataSynOU$Regime5, phy = ttreeA, type = "discrete", model = "ER")
MLE_idx5 <- apply(anc_recon5$lik.anc, 1, which.max)
MLE_hypothesis5 <- colnames(anc_recon5$lik.anc)[MLE_idx5]

## Create a new phylogeny object
ttreeA5 <- ttreeA
ttreeA5$node.label <- MLE_hypothesis5

## Fit the model
m5 <- slouch.fit(phy = ttreeA5,
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60), 
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw, 
                 fixed.fact = DataSynOU$Regime5,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m5)
summary(m5) #not great


##HYPOTHESIS 6: Only basal anomodonts hearing, but also shift at therapsids
anc_recon6 <- ace(DataSynOU$Regime6, phy = ttreeA, type = "discrete", model = "ER")
MLE_idx6 <- apply(anc_recon6$lik.anc, 1, which.max)
MLE_hypothesis6 <- colnames(anc_recon6$lik.anc)[MLE_idx6]

## Create a new phylogeny object
ttreeA6 <- ttreeA
ttreeA6$node.label <- MLE_hypothesis6

## Fit the model
m6 <- slouch.fit(phy = ttreeA6,
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60), 
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw,
                 fixed.fact = DataSynOU$Regime6,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m6)
summary(m6) #about as good as H1 and H3. Intercept increases in therapsids and again in anomodonts

######

##HYPOTHESIS 7: Only all anomodonts hearing, but also shift at therapsids
anc_recon7 <- ace(DataSynOU$Regime7, phy = ttreeA, type = "discrete", model = "ER")
MLE_idx7 <- apply(anc_recon7$lik.anc, 1, which.max)
MLE_hypothesis7 <- colnames(anc_recon7$lik.anc)[MLE_idx7]

## Create a new phylogeny object
ttreeA7 <- ttreeA
ttreeA7$node.label <- MLE_hypothesis7

## Fit the model
m7 <- slouch.fit(phy = ttreeA7, 
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60), 
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw, 
                 fixed.fact = DataSynOU$Regime7,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m7)
summary(m7) #slightly better fit that H6

#################

#just gonna play around with changing a hypothesis 8 in the csv
#try different things to see what patterns arise
DataSynOU <- read_csv("RL_allometry_forR_finalOU.csv")

anc_recon8 <- ace(DataSynOU$Regime8, phy = ttreeA, type = "discrete", model = "ER")
MLE_idx8 <- apply(anc_recon8$lik.anc, 1, which.max)
MLE_hypothesis8 <- colnames(anc_recon8$lik.anc)[MLE_idx8]

## Create a new phylogeny object
ttreeA8 <- ttreeA
ttreeA8$node.label <- MLE_hypothesis8

## Fit the model
m8 <- slouch.fit(phy = ttreeA8, #use the parameters from this run for everything
                 hillclimb = F,
                 vy_values = seq(0.001, 1, length.out = 60),
                 hl_values = seq(0.1, 50, length.out = 60), 
                 species = DataSynOU$Species,
                 response = DataSynOU$logRL,
                 direct.cov = DataSynOU$logJaw, #I was doing direct cov but I think it should be random since it's also evolving
                 fixed.fact = DataSynOU$Regime8,
                 mv.response = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label)),#accounting for 5% measurement error
                 mv.direct.cov = rep(0.055 * var(DataSynOU$logRL), length(ttreeA$tip.label))) #plus 0.5% intraspecific variation
plot(m8)
summary(m8)

#The following are notes from my initial runs with this dataset
#Probably not accurate now but a good reminder of things to try if we get the model tweaked to our liking

#moderately high support when each therapsid subclade got its own optimum
#pretty high support for ALL NON-BIARMOSUCHIAN THERAPSIDS AS ONE REGIME (got slightly worse when I also made dinocephalians 0)
#didn't get as high of support when excluding biarmosuchians from therapsid regime in nb only with therapsid shift setup
#pretty good support for bidentalians going to own regime after nb shift, but not as high as going back to therapsid state
#similar support for bidentalians going to a new regime that's the same as dinocephalians

#modeling Kannemeyeriids and cryptodonts as reversions to therapsid regime gets even higher support that a reversion at Bidentalia.
#even higher when Ulemica is lumped in with other therapsids
#this implies that it's a change at the base of Dicynodontia and carries through to the smaller bidentalians
#which implies that it doesn't match up with the surface topology of the RL
#other than the small cryptodonts this does still make sense for thickness of the RL
#maybe I should try setting regimes based on RL thickness

#ok nope figured it out. If I have a shift at the base of Dicynodontia and then a reversion in Endothiodon and Bidentalia
#it has the highest support of any model. It still wants Ulemica to be with the other therapsids which is fine with me
#but this really supports the idea that larger dicynodonts lost their sense of hearing
#assuming basal bidentalians actually were large and therefore subjected to this shift

#brownian motion model instead of OU

mb <- brown.fit(phy = ttreeA,
                hillclimb = TRUE,
                species = DataSynOU$Species,
                response = DataSynOU$logRL,
                direct.cov = DataSynOU$logJaw)
plot(mb)
summary(mb)
#real bad lol


#######################################################################################################################

#Run with largest chameleons of each species
#note that this dataset is slightly expanded from the intraspecific allometry dataset, so it uses a different trimmed tree
setwd("D:/Work backup/Research Projects/Chameleons/Measurements")

Tolley2013Trees <- "Tolley2013trimmedL.trees" #tree trimmed to just species with largest specimen data
cham.tree<-read.nexus(Tolley2013Trees)
#To gaze upon its beauty: 
#plotTree(cham.tree,ftype="i",fsize=0.6,lwd=1)


DataChamL <- read_csv("MasterLargest.csv")

rownames(DataChamL) <- DataChamL$Species #check that dataset and tree names are identical
obj<-name.check(cham.tree,DataChamL)
cham.tree<-drop.tip(cham.tree, obj$tree_not_data) #this line is bein wierd idk the dataset is fine though

#plotTree(cham.tree,ftype="i",fsize=0.6,lwd=1)

## Reorder dataframe so the species labels match the tree tip labels
DataChamL <- DataChamL[match(cham.tree$tip.label, DataChamL$Species), ]

##HYPOTHESIS 0 all the same
MLE_hypothesis0 <- as.character(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) #just set this for everything I guess lmao
## Create a new phylogeny object
CT0 <- cham.tree
CT0$node.label <- MLE_hypothesis0
## Fit the model without a fixed factor for the selective regime because Regime0 is all non hearing
mC0 <- slouch.fit(phy = CT0,
                  hillclimb = F,
                  vy_values = seq(0.001, 1, length.out = 60), #this is giving negative slope even when I increase hl
                  hl_values = seq(2, 30, length.out = 60), 
                  species = DataChamL$Species,
                  response = DataChamL$logpter,
                  direct.cov = DataChamL$logbsl,
                  mv.response = rep(0.055 * var(DataChamL$logpter), length(cham.tree$tip.label)),#accounting for 5% measurement error
                  mv.direct.cov = rep(0.055 * var(DataChamL$logbsl), length(cham.tree$tip.label))) #plus 0.5% intraspecific variation
plot(mC0)
summary(mC0)

#######

##REGIME 1 Hearing different regime
MLE_hypothesis1 <- as.character(c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1)) #just set this for everything I guess lmao
## Create a new phylogeny object
CT1 <- cham.tree
CT1$node.label <- MLE_hypothesis1
plot(CT1, show.node.label = TRUE) #make sure internal nodes were assigned correctly
mC1 <- slouch.fit(phy = CT1,
                  hillclimb = F,
                  vy_values = seq(0.001, 1, length.out = 60),
                  hl_values = seq(2, 30, length.out = 60), 
                  species = DataChamL$Species,
                  response = DataChamL$logpter,
                  direct.cov = DataChamL$logbsl,
                  mv.response = rep(0.055 * var(DataChamL$logpter), length(cham.tree$tip.label)),#accounting for 5% measurement error
                  mv.direct.cov = rep(0.055 * var(DataChamL$logbsl), length(cham.tree$tip.label)), #plus 0.5% intraspecific variation
                  fixed.fact = DataChamL$Regime1)
plot(mC1)
summary(mC1) #worse than null
regimeplot(mC1)

#this doesn't perform as well, but I think that's because it's assuming the same optimum for all the non-hearing chameleons
#when really the point is that they're exploring lots of different optima

#brownian motion model instead of OU

mbC <- brown.fit(phy = cham.tree,
                 hillclimb = TRUE,
                 species = DataChamL$Species,
                 response = DataChamL$logpter,
                 direct.cov = DataChamL$logbsl)
plot(mbC)
summary(mbC) #this one performs pretty badly

#######

##HYPOTHESIS 2: hearing group is its own regime, the genera Furcifer and Trioceros as separate regimes (except hearing taxa)
MLE_hypothesis2 <- as.character(c(0,0,0,0,2,2,2,0,3,3,3,3,3,1,1,1,1,1,1,1)) #just set this for everything I guess lmao
## Create a new phylogeny object
CT2 <- cham.tree
CT2$node.label <- MLE_hypothesis2
plot(CT2, show.node.label = TRUE) #make sure internal nodes were assigned correctly
mC2 <- slouch.fit(phy = CT2,
                  hillclimb = TRUE,
                  species = DataChamL$Species,
                  response = DataChamL$logpter,
                  direct.cov = DataChamL$logbsl,
                  fixed.fact = DataChamL$Regime2,
                  mv.response = rep(0.055 * var(DataChamL$logpter), length(cham.tree$tip.label)),#accounting for 5% measurement error
                  mv.direct.cov = rep(0.055 * var(DataChamL$logbsl), length(cham.tree$tip.label))) #plus 0.5% intraspecific variation
plot(mC2)
summary(mC2) #better than null
regimeplot(mC2)

#played around and got very good support for Trioceros and Furcifer+Calumma having their own non hearing regimes and
#hearing group being its own regime
#for some reason it likes Calumma to be in the same regime as Furcifer but not Kinyongia with Trioceros idk

########

#the following is just to mess around with a Hypothesis 3

DataChamL <- read_csv("MasterLargest.csv")
#Reorder dataframe so the species labels match the tree tip labels
DataChamL <- DataChamL[match(cham.tree$tip.label, DataChamL$Species), ]

MLE_hypothesis3 <- as.character(c(0,0,0,0,2,2,2,0,3,3,3,3,3,3,1,1,1,1,1,1)) #just set this for everything I guess lmao
## Create a new phylogeny object
CT3 <- cham.tree
CT3$node.label <- MLE_hypothesis3
plot(CT3, show.node.label = TRUE) #make sure internal nodes were assigned correctly
mC3 <- slouch.fit(phy = CT3,
                  hillclimb = TRUE,
                  species = DataChamL$Species,
                  response = DataChamL$logpter,
                  direct.cov = DataChamL$logbsl,
                  fixed.fact = DataChamL$Regime3,
                  mv.response = rep(0.055 * var(DataChamL$logpter), length(cham.tree$tip.label)),#accounting for 5% measurement error
                  mv.direct.cov = rep(0.055 * var(DataChamL$logbsl), length(cham.tree$tip.label))) #plus 0.5% intraspecific variation
plot(mC3)
summary(mC3)
regimeplot(mC3)
#old notes from initial runs, might not still be accurate but something to play with in the future
#way less support when the hearing taxa in these genera are lumped into the non hearing regimes for their genus
#best support when T. hoehnelii is its own regime
