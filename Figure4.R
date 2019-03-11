##FlowSOM based analysis of memory cells#MYL
##Figure 4: Expression of surface proteins on acute, latent and chronic memory
## Stain list:
# CD8-V500, CD90.1-Pacific Blue, CD27-PECy7, CD62L-APCCy7, CD127-FITC, KLRG1-PE, CD90.1-Pacific Blue, GP33-APC

#Install packages
library(flowCore)
library(flowTrans)
library(ggplot2) 
library(tidyr)
library(plyr)
library(RColorBrewer)
library(FlowSOM)

##import data
#set up working directory
setwd("~/Dropbox/Analysis of T cell fate/Publication/R code/surface_population1") 
#read the written Rda file
load("s5_memory.Rda") 
#take out mouse label and infection type for data transformation
ALL<-s5_with_label[,c(-5,-6)]

##ArcSinh transformation using flowTrans package
ALL<-as.matrix(ALL)
#convert matrix to flow frame
ALL<-flowFrame(ALL)
#ArcSinh Transformation
ALL_trans<-flowTrans(dat=ALL,fun="mclMultivArcSinh",colnames(ALL), n2f=FALSE,parameters.only = FALSE)
#extract the transformed data
ALL_after_transformation<-exprs(ALL_trans$result) 
#scale the transformed data
ALL<-as.matrix(scale(ALL_after_transformation))

##FCS conversion for FlowSOM
# convert the pooled cells to FCS
ALL_FCS <- flowFrame(ALL)
#add infection label back
ALL_with_label <- cbind(as.data.frame(ALL),infection=as.data.frame(s5_with_label)[,6]) 
#extract acute memory data
ARM <- ALL_with_label[ALL_with_label$infection == "ARM",][,-5]
#extract latent memory data
MHV <- ALL_with_label[ALL_with_label$infection == "MHV",][,-5]
#extract chronic memory data
CL13 <- ALL_with_label[ALL_with_label$infection == "CL13",][,-5] 
#extract naive data
Naive <- ALL_with_label[ALL_with_label$infection == "Naive",][,-5] 

#convert data to FCS
ARM_FCS <- flowFrame(as.matrix(ARM))
MHV_FCS <- flowFrame(as.matrix(MHV)) 
CL13_FCS <- flowFrame(as.matrix(CL13)) 
Naive_FCS <- flowFrame(as.matrix(Naive)) 

##FlowSOM based analysis
set.seed(1234)
#Read the pooled population into FlowSom
fSOM_ALL <- ReadInput(ALL_FCS,compensate = FALSE, transform = FALSE) 
#Build Self-Organizing Map
fSOM_ALL <- BuildSOM(fSOM_ALL, xdim=5,ydim=5,rlen=100) 
#Build Minimal Spanning Tree
fSOM_ALL <- BuildMST(fSOM_ALL) 
#Plot MST of the pooled population
PlotStars(fSOM_ALL,legend = FALSE) 

##Plot MST by infection
#Map acute memory data to the SOM grid 
fSOM_ARM <- NewData(fSOM_ALL,ARM_FCS) 
PlotStars(fSOM_ARM, legend = FALSE)
#Map latent memory data to the SOM grid (standard formed by total population)
fSOM_MHV <- NewData(fSOM_ALL,MHV_FCS)
PlotStars(fSOM_MHV, legend = FALSE)
#Map chronic memory data to the SOM grid (standard formed by total population)
fSOM_CL13 <- NewData(fSOM_ALL,CL13_FCS) 
PlotStars(fSOM_CL13, legend = FALSE)
#Map naive data to the SOM grid (standard formed by total population)
fSOM_Naive <- NewData(fSOM_ALL,Naive_FCS) 
PlotStars(fSOM_Naive,legend = FALSE)




