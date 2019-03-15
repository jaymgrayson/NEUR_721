# Dimensionality reduction comparison
# Jason Grayson
# Started 9-17-18

rm(list=ls())


load("file=s5_memory.Rda")
Memory_With_Label<-s5_with_label
rm(s5_with_label)
str(Memory_With_Label)
Memory_No_Label<-Memory_With_Label[,-c(5:6)]

set.seed(1234)
library(flowCore)
library(flowTrans)
ALL<-Memory_No_Label
ALL<-as.matrix(ALL)
ALL<-flowFrame(ALL) #convert matrix to flow frame
ALL_trans<-flowTrans(dat=ALL,fun="mclMultivArcSinh",colnames(ALL), n2f=FALSE,parameters.only = FALSE)
ALL_after_transformation<-exprs(ALL_trans$result)
ALL<-as.data.frame(scale(ALL_after_transformation)) #scale the transformed data

library(psych)
library(ggplot2)

fa.parallel(ALL,fa="pc",n.iter=100,show.legend = FALSE,
            main="screen plot with parallel analysis")
ALL_PCA<-principal(ALL,2)
print(ALL_PCA)
ALL_pred <- as.data.frame(ALL_PCA$scores)
ALL_pred <- cbind(ALL_pred, infection = Memory_With_Label$infection)
colnames(ALL_pred) <- c("PC1", "PC2","infection")

PCA_Plot<-ggplot(ALL_pred, aes(x=PC1,y=PC2,color=infection))+geom_point(size=0.5)+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(colour = guide_legend(override.aes = list(size=4)))+ggtitle("PCA of Memory Surface Data")
PCA_Plot

library(Rtsne)
All_tsne<-Rtsne(ALL,dims = 2, perplexity=30, verbose=TRUE, max_iter = 5000)
tsne_time<-1596
save(All_tsne,file="All_tsne.Rda")
load(file = "All_tsne.Rda")
Flow_tsne<-All_tsne$Y
Flow_tsne<-as.data.frame(Flow_tsne)
Flow_tsne<-cbind(Flow_tsne,Memory_With_Label$infection)
colnames(Flow_tsne)<-c("V1", "V2","infection")

tsne_Plot<-ggplot(Flow_tsne, aes(x=V1,y=V2,color=infection))+geom_point(size=0.5)+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(colour = guide_legend(override.aes = list(size=4)))+ggtitle("tsne of Memory Surface Data")
tsne_Plot

library(umap)
Flow_umap<-umap(ALL)
x<-120
Flow_umap2<-Flow_umap$layout
Flow_umap2<-as.data.frame(Flow_umap2)
Flow_umap2<-cbind(Flow_umap2,Memory_With_Label$infection)
colnames(Flow_umap2)<-c("V1", "V2","infection")
save(Flow_umap2,file="Flow_umap.Rda")
load("~/Desktop/Critical Docs/Meetup Talk/Flow_umap.Rda")
umap_Plot<-ggplot(Flow_umap2, aes(x=V1,y=V2,color=infection))+geom_point(size=0.5)+theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ guides(colour = guide_legend(override.aes = list(size=4)))+ggtitle("umap of Memory Surface Data")
umap_Plot
