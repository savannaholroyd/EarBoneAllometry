library(phylolm)
library(readr)
library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(phytools)
require(devtools)
install_github("JeroenSmaers/evomap")
library(evomap)

#######################################################

DataANC <- read_csv("RL_allometry_ANCOVA.csv") #copy and paste the groups you want to test into this spreadsheet from RL_allometry_forR_final.csv


#pANCOVA sensu Smaers and Rohlf, 2016

data<-as.data.frame(DataANC); tree<-AllomeTree
plotTree(tree)
rownames(data) <- DataANC$Species
Y<-"logRL"; X<-"logJaw"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]

tree<-treedata(tree,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
colnames(data)<-c("Dependent","Independent") 

Therapsids<-getTips(tree,findMRCA(tree,c("Procynosuchus", "Trirachodon"))) #add in manually. Sorry. #"Procynosuchus", "Trirachodon"; "Ulemica", "Placerias";"Sinophoneus", "Moschops"; "Biarmosuchus_tener", "Paraburnetia"; "Aelurognathus", "Aloposaurus"
Pely<-getTips(tree,findMRCA(tree,c("Palaeohatteria", "Dimetrodon_limbatus")))
Pely <- as.integer(c(25, 26, 3, 24, 1, 2)) #set whichever weren't in Therapsids

plot(data$Dependent~data$Independent,col="white",xlab="", ylab="") 
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=Therapsids,col="green") 
pGLS.plotGrade("Dependent","Independent",data,tree,model="BM",group=Pely,col="grey")

grpS<-rep("A",length(rownames(data))) 
grpS[Therapsids]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(data)

#For differences in intercept: 
grpI<-rep("A",length(rownames(data))) 
grpI[Therapsids]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

Model<-model.matrix(as.formula(Dependent~Independent),data)
#(1) Differences in slopes, holding intercept constant: 
Model_S<-model.matrix(as.formula(Dependent~grpS:Independent),data) 

#(2) Differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Dependent~grpI + Independent),data) 

#(3) Differences in slopes and differences in intercept: 
Model_SI<-model.matrix(as.formula(Dependent~grpI + grpS:Independent),data)


#(1) Differences in slopes, holding intercept constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_S)

#(2) Differences in intercept, holding slopes constant:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_I)

#(3) Differences in slopes and differences in intercept:
gls.ancova(Dependent~Independent,vcv(tree),Model,Model_SI)