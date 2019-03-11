##Cluster Analysis of memory cells#MYL
##Figure 3: Expression of surface proteins on acute, latent and chronic memory
## Stain list:
# CD8-V500, CD90.1-Pacific Blue, CD27-PECy7, CD62L-APCCy7, CD127-FITC, KLRG1-PE, CD90.1-Pacific Blue, GP33-APC
## Mouse label:
# Acute memory:M01-M13; Latent Memory:M14-M19; Chronic Memory:M20-M28; Naive P14:M29-M31


#Install packages
library(flowCore)
library(flowTrans)
library(ggplot2)
library(tidyr)
library(plyr)
library(RColorBrewer)

##import data
#set up working directory
setwd("~/Dropbox/Analysis of T cell fate/Publication/R code/surface_population1")
#read the written Rda file
load("s5_memory.Rda")
#take out mouse label and infection type for data transformation
ALL<-s5_with_label[,c(-5,-6)]

##ArcSinh transformation using flowTrans package
#convert dataframe to matrix
ALL<-as.matrix(ALL)
#convert matrix to flow frame
ALL<-flowFrame(ALL)
#ArcSinh Transformation
ALL_trans<-flowTrans(dat=ALL,fun="mclMultivArcSinh",colnames(ALL), n2f=FALSE,parameters.only = FALSE) 
plot(ALL_trans$result,c("CD127","KLRG1"))
#contour plot of CD127 by KLRG1
contour(ALL_trans$result,c("CD127","KLRG1"),add=TRUE)
#extract the transformed data
ALL_after_transformation<-exprs(ALL_trans$result)  
#scale the transformed data
ALL<-as.data.frame(scale(ALL_after_transformation)) 

##determine optimal number of clusters by WSS plot
set.seed(2L)
wss <- (nrow(ALL)-1)*sum(apply(ALL,2,var))
#plot k by wss to find the elbow of the plot
for (i in 2:15) wss[i] <- sum(kmeans(ALL,algorithm = "Lloyd", centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") 
k <- 6

##K-means clustering
set.seed(1L) 
ALL_cluster <- kmeans(ALL,k,algorithm = "Lloyd", nstart = 100) 
#Add cluster number to the dataset
ALL_after_cluster <- cbind(ALL,cluster=ALL_cluster$cluster) 

##Visualize the characteristics of clusters
#stacked density plot to show the characteristic of each cluster
#calculate the density of fluorescent molecules depending on the cluster
Density=NULL
for (i in 1:5){
  #take out cluster and combine with each fluorscence in each round
  my.data <- as.data.frame(cbind( V1 = ALL_after_cluster[,i], cluster=ALL_after_cluster[,ncol(ALL_after_cluster)])) 
  #calculate density of chosen molecules depending on cluster
  res <- dlply(my.data, .(cluster), function(x) density(x$V1))
  # take out the density and value from the result list
  dd <- ldply(res, function(z){data.frame(Values = z[["x"]], 
                                          V1_density = z[["y"]],
                                          V1_count = z[["y"]]*z[["n"]])}) 
  # adjust an offset
  dd$offest <- -0.4*dd$cluster 
  # add the offset to the density to avoid verical overlap
  dd$V1_density_offest=dd$V1_density+dd$offest 
  Fluor <- colnames(ALL_after_cluster)[i]
  #generate a dataset combining offset, density and original value
  DD <- cbind(dd$cluster,dd$Values,dd$V1_density_offest,dd$offest) 
  colnames(DD) <- c("cluster",paste(Fluor,"Value"),paste(Fluor,"density_offset"),paste(Fluor,"offset"))
  #combine results of different molecules
  Density <- cbind(Density,DD) 
}


##separate the dataset by molecules and add label
CD127_density <- as.data.frame(Density[,1:4])
colnames(CD127_density)<-c("cluster","FI","density","offset")

KLRG1_density <- as.data.frame(Density[,5:8])
colnames(KLRG1_density)<-c("cluster","FI","density","offset")

CD27_density <- as.data.frame(Density[,9:12])
colnames(CD27_density)<-c("cluster","FI","density","offset")

CD62L_density <- as.data.frame(Density[,13:16])
colnames(CD62L_density)<-c("cluster","FI","density","offset")

#plot the characteristics of 6 clusters (Figure 3C)
ggplot(CD127_density, aes(x=FI, y=density,group=cluster)) + 
  geom_line(color="white") + geom_ribbon(aes(FI, ymin=offset,ymax=density),alpha=0.5) + 
  scale_y_continuous(breaks=NULL) +
  ggtitle("") + xlab("") + ylab("") +
  theme(panel.background = element_rect(fill="white",colour = "black", size = 0.75),axis.text.x=element_blank(), 
        axis.text.y=element_blank(), axis.title.x=element_blank(),
        axis.ticks=element_blank()) + theme_bw() 

ggplot(KLRG1_density, aes(x=FI, y=density,group=cluster)) + 
  geom_line(color="white") + geom_ribbon(aes(FI, ymin=offset,ymax=density),alpha=0.5) + 
  scale_y_continuous(breaks=NULL) +
  ggtitle("") + xlab("") + ylab("") +
  theme(panel.background = element_rect(fill="white",colour = "black", size = 0.75),axis.text.x=element_blank(), 
        axis.text.y=element_blank(), axis.title.x=element_blank(),
        axis.ticks=element_blank()) + theme_bw()

ggplot(CD27_density, aes(x=FI, y=density,group=cluster)) + 
  geom_line(color="white") + geom_ribbon(aes(FI, ymin=offset,ymax=density),alpha=0.5) + 
  scale_y_continuous(breaks=NULL) +
  ggtitle("") + xlab("") + ylab("") +
  theme(panel.background = element_rect(fill="white",colour = "black", size = 0.75),axis.text.x=element_blank(), 
        axis.text.y=element_blank(), axis.title.x=element_blank(),
        axis.ticks=element_blank()) + theme_bw()

ggplot(CD62L_density, aes(x=FI, y=density,group=cluster)) + 
  geom_line(color="white") + geom_ribbon(aes(FI, ymin=offset,ymax=density),alpha=0.5) + 
  scale_y_continuous(breaks=NULL) +
  ggtitle("") + xlab("") + ylab("") +
  theme(panel.background = element_rect(fill="white",colour = "black", size = 0.75), axis.text.x=element_blank(), 
        axis.text.y=element_blank(), axis.title.x=element_blank(),
        axis.ticks=element_blank()) + theme_bw()

##add mouse label to the dataset
ALL_after_cluster<- cbind(ALL_after_cluster, mouse=as.matrix(s5_with_label)[,ncol(s5_with_label)-1])

##create a function to calculate percentage of cells in each cluster
r=0
cluster_ratio<-function(x){
  for (i in 1:k){
    count<-length(x[x==i])
    ratio<-round(count/length(x),4)*100
    r[i]=ratio
  }
  return(r)
}

#apply the cluster_ratio function to calculate the percentage of cell in each cluster
result_of_all<-as.matrix(aggregate(ALL_after_cluster$cluster, by=list(ALL_after_cluster$mouse), FUN=cluster_ratio))
#rename the dataset
colnames(result_of_all)<-c("mouse","Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6") 
print(result_of_all)

##create infection label
infection<-c(rep("ARM",13),rep("MHV",6),rep("CL13",9),rep("Naive",3)) #"ARM",LCMV-ARM;"MHV",MHV68-GP33;"CL13",LCMV-CLONE 13
#add the label to result_of_all
result_of_all<-as.data.frame(cbind(result_of_all,infection=infection)) 

##remove infection label
result_of_all_1 <- result_of_all[,c(-1,-(k+2))]

##create a convert_to_numeric function
numeric_conversion<- function(x){
  for (i in 1:ncol(x)){
    x[,i] = as.numeric(as.character(x[,i]))
  }
  return(x)
}

##convert factor in result_of_all_1 to numeric variables
result_of_all <- cbind(numeric_conversion(result_of_all_1),infection = result_of_all[,k+2])

#Aggregate the result by infections for analysis
#gather the data using tidyr package
ALL_gather <- gather(result_of_all,"cluster","percentage_of_cluster", -infection) 
#calculate the mean percentage of cells in each cluster
mean_for_rose<-as.data.frame(aggregate(ALL_gather$percentage_of_cluster, by=list(ALL_gather$infection,ALL_gather$cluster), FUN=mean))
#calculate the standard deviation
sd_for_rose<-as.data.frame(aggregate(ALL_gather$percentage_of_cluster, by=list(ALL_gather$infection,ALL_gather$cluster), FUN=sd))
#generate a dataset that have mean and sd
ALL_for_rose <-cbind(mean_for_rose,sd_for_rose[,3])
colnames(ALL_for_rose) <- c("infection","cluster","percentage_of_cluster","sd")


#generate a barplot showing the distribution of CD90.1 cells to the 6 clusters
#plot all infections
ALL_bar <- ggplot(ALL_for_rose, aes(x=cluster, y=percentage_of_cluster,group=cluster, fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") +
  geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd),colour="black") +facet_grid(.~infection)+ scale_fill_brewer() + scale_color_brewer()
#plot each infection individually
ARM_bar <- ggplot(ALL_for_rose[ALL_for_rose$infection=="ARM",], aes(x=cluster, y=percentage_of_cluster,group=cluster, fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + 
  geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd),colour = "black")+ scale_fill_brewer() + scale_color_brewer()
CL13_bar <- ggplot(ALL_for_rose[ALL_for_rose$infection=="CL13",], aes(x=cluster, y=percentage_of_cluster, group= cluster,fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + 
  geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd), colour = "black")+ scale_fill_brewer() + scale_color_brewer()
MHV_bar <- ggplot(ALL_for_rose[ALL_for_rose$infection=="MHV",], aes(x=cluster, y=percentage_of_cluster, grouop=cluster ,fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + 
  geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd), colour = "black")+ scale_fill_brewer() + scale_color_brewer()
Naive_bar <- ggplot(ALL_for_rose[ALL_for_rose$infection=="Naive",], aes(x=cluster, y=percentage_of_cluster, grouop=cluster ,fill = cluster)) + geom_bar(width = 1, stat="identity",colour="black") + 
  geom_errorbar(aes(ymin= percentage_of_cluster, ymax=percentage_of_cluster + sd), colour = "black")+ scale_fill_brewer() + scale_color_brewer()

#convert the barplot to rose diagram
ALL_bar + scale_y_continuous(breaks = c(0,10,20,30,40,50,60)) + coord_polar() + labs(x = "", y = "") + ggtitle("") +
  facet_wrap(~ infection, nrow = 2) 
#plot each infection individually (Figure 2D)
ARM_bar + scale_y_continuous(breaks=c(0,10,20,30,40,50),limits = c(0,55)) + coord_polar() + labs(x = "", y = "") + ggtitle("") + theme_bw() + theme(legend.position="none",panel.background = element_rect(fill="white",colour = "black", size = 0.75), 
                                                                                                                                                 axis.text.x=element_blank(), 
                                                                                                                                                 axis.text.y=element_blank(),
                                                                                                                                                 axis.ticks=element_blank(),
                                                                                                                                                 panel.grid.major = element_line(colour = "gray",size=0.75),
                                                                                                                                                 axis.title.x=element_blank())
CL13_bar + scale_y_continuous(breaks=c(0,10,20,30,40,50),limits = c(0,55)) + coord_polar() + labs(x = "", y = "") + ggtitle("") + theme_bw() + theme(legend.position="none",panel.background = element_rect(fill="white",colour = "black", size = 0.75), 
                                                                                                                                                  axis.text.x=element_blank(), 
                                                                                                                                                  axis.text.y=element_blank(),
                                                                                                                                                  axis.ticks=element_blank(),
                                                                                                                                                  panel.grid.major = element_line(colour = "gray",size=0.75),
                                                                                                                                                  axis.title.x=element_blank())
MHV_bar + scale_y_continuous(breaks=c(0,10,20,30,40,50),limits = c(0,55)) + coord_polar() + labs(x = "", y = "") + ggtitle("")  + theme_bw() + theme(legend.position="none",panel.background = element_rect(fill="white",colour = "black", size = 0.75), 
                                                                                                                                                  axis.text.x=element_blank(), 
                                                                                                                                                  axis.text.y=element_blank(),
                                                                                                                                                  axis.ticks=element_blank(),
                                                                                                                                                  panel.grid.major = element_line(colour = "gray",size=0.75),
                                                                                                                                                  axis.title.x=element_blank())

Naive_bar + scale_y_continuous(limits = c(0,65)) + coord_polar() + labs(x = "", y = "") + ggtitle("")  + theme_bw() + theme(legend.position="none",panel.background = element_rect(fill="white",colour = "black", size = 0.5), 
                                                                                                                            axis.text.x=element_blank(), 
                                                                                                                            axis.text.y=element_blank(),
                                                                                                                            axis.ticks=element_blank(),
                                                                                                                            panel.grid.major = element_line(colour = "gray",size=0.5),
                                                                                                                            axis.title.x=element_blank())
##ANOVA to compare the effects of various infections on expression of surface markers on CD8
#Cluster 1
AOV_C1<- aov (Cluster1~infection, data = result_of_all)
summary(AOV_C1)
#multiple comparison
TukeyHSD(AOV_C1) 
#Cluster 2
AOV_C2<- aov (Cluster2~infection, data = result_of_all)
summary(AOV_C2)
TukeyHSD(AOV_C2)
#Cluster 3
AOV_C3<- aov (Cluster3~infection, data = result_of_all)
summary(AOV_C3)
TukeyHSD(AOV_C3)
#Cluster 4
AOV_C4<- aov (Cluster4~infection, data = result_of_all)
summary(AOV_C4)
TukeyHSD(AOV_C4)
#Cluster 5
AOV_C5<- aov (Cluster5~infection, data = result_of_all)
summary(AOV_C5)
TukeyHSD(AOV_C5)
#Cluster 6
AOV_C6<- aov (Cluster6~infection, data = result_of_all)
summary(AOV_C6)
TukeyHSD(AOV_C6)

                                                                                                                      
                                                                                                                            