#By: Theresa Christiansen
#Note: some debugging done with ChatGPT bc I am not great at R syntax
#Citation: Wilkerson, M. D., & Hayes, D. N. (2010). ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics (Oxford, England), 26(12), 1572–1573. https://doi.org/10.1093/bioinformatics/btq170
#This file is consensus clustering on each of the two datasets and the combined dataset
#It also does Kmeans clustering on the entire dataset ordered by Differential Gene expression (case vs. ctrl)
#The top 10 most differentially expressed genes per cluster are labeled

#Install Consensus Cluster Plus
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

#Import GSE20680.csv
d <- GSE20680
d <- d[-1, ]
d <- d[, -c(1, 2)]
print(d[1:5, 1:5])
d <- apply(d, 2, as.numeric)

#From the tutorial: calculate MAD of each gene and pick the top 5000 w/most variability
mads <- apply(d, 2, mad, na.rm = TRUE)
d <- d[, rev(order(mads))[1:5000]]
d <- sweep(d, 1, apply(d, 1, median, na.rm = TRUE))
print(d[1:5, 1:5])

#Do Consensus Clustering with Hierarchical and Kmeans
library(ConsensusClusterPlus) 
title=tempdir() 
results = ConsensusClusterPlus(d,maxK=5,reps=10,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot=NULL)
icl = calcICL(results,title=title,plot=NULL)

resultskm = ConsensusClusterPlus(d,maxK=5,reps=100,pItem=0.8,pFeature=1,title=title,clusterAlg="km",distance="pearson",seed=1262118388.71279,plot=NULL)
icl = calcICL(resultskm,title=title,plot=NULL)

#Consensus Clustering with GSE20681_preprocessed.csv
k <- GSE20681_preprocessed
k <- k[-1, ]
k <- k[, -c(1)]
k <- apply(k, 2, as.numeric)
tk <- t(k)
print(tk[1:5, 1:5])

#From the tutorial: calculate MAD of each gene and pick the top 5000 w/most variability
mads <- apply(tk, 2, mad, na.rm = TRUE)
tk <- tk[, rev(order(mads))[1:1000]]
tk <- sweep(tk, 1, apply(tk, 1, median, na.rm = TRUE))
print(tk[1:5, 1:5])

#Do Consensus Clustering with Hierarchical and Kmeans
library(ConsensusClusterPlus) 
title=tempdir() 
results = ConsensusClusterPlus(tk,maxK=5,reps=10,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot=NULL)
icl = calcICL(results,title=title,plot=NULL)

resultskm = ConsensusClusterPlus(tk,maxK=5,reps=100,pItem=0.8,pFeature=1,title=title,clusterAlg="km",distance="pearson",seed=1262118388.71279,plot=NULL)
icl = calcICL(resultskm,title=title,plot=NULL)

#Make Combined Dataset
k <- GSE20681_preprocessed
k <- k[-1, ]
k <- k[, -c(1)]
k <- apply(k, 2, as.numeric)
tk <- t(k)
d <- GSE20680
d <- d[-1, ]
d <- d[, -c(1, 2)]
d <- apply(d, 2, as.numeric)
print(d[1:5, 1:5])

combined <- rbind(tk, d)
print(combined[1:5, 1:5])
write.csv(data, file = "combinedGSE.csv", row.names = FALSE)

#Consensus Clustering Results for Combined Dataset 
combined <- combinedGSE
combined <- combined[-1, ]
combined <- as.data.frame(lapply(combined, as.numeric))
mads <- apply(combined, 2, mad, na.rm = TRUE)
combined <- combined[, rev(order(mads))[1:5000]]
combined <- sweep(combined, 1, apply(combined, 1, median, na.rm = TRUE))
print(combined[1:5, 1:5])
combined <- as.matrix(combined)
combinedt = as.dist(1-cor(combined,method="pearson"))

#Call Consensus Clustering for Hierarchical and Kmeans
library(ConsensusClusterPlus) 
title=tempdir() 
results = ConsensusClusterPlus(combinedt,maxK=5,reps=10,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot=NULL)
icl = calcICL(results,title=title,plot=NULL)

resultskm = ConsensusClusterPlus(combinedt,maxK=5,reps=100,pItem=0.8,pFeature=1,title=title,clusterAlg="km",distance="pearson",seed=1262118388.71279,plot=NULL)
icl = calcICL(resultskm,title=title,plot=NULL)

#Clustering Assignments for top 5000 differentially expressed genes
#Import ctrl_vs_case.csv
genes = ctrl_vs_case[1:20000,2] #from Romy's analysis
combined <- combinedGSE
combined <- combined[-1, ]

#adjust column names so they start at V1
new_column_names <- paste0("V", seq_along(combined)) 
colnames(combined) <- new_column_names

#create distribution for clustering
diffgenes <- combined[,c(genes)]
diffgenes <- sweep(diffgenes, 1, apply(combined, 1, median, na.rm = TRUE))
print(diffgenes[1:5, 1:5])
diffgenes <- as.matrix(diffgenes)
diffgenest = as.dist(1-cor(diffgenes,method="pearson"))

#Run Hierarchical Clustering
hc1 <- hclust(diffgenest, method = "complete") 
plot(hc1)
clusters <- cutree(hc1, k = 3)  #k is from consensus clustering
gene_clusters <- data.frame(Gene = genes, Cluster = clusters)
write.csv(gene_clusters, file = "ClusteredDEGs2.csv", row.names = FALSE)

#Print top 10 gene names for each cluster
cluster1 = gene_clusters$Gene[gene_clusters$Cluster == 1]
cluster2 = gene_clusters$Gene[gene_clusters$Cluster == 2]
cluster3 = gene_clusters$Gene[gene_clusters$Cluster == 3]
print(cluster1[1:10])
print(cluster2[1:10])
print(cluster3[1:10])

#Make distance matrix for Kmeans clustering
diffgeneskm <- combined[,genes]
diffgeneskm <- sweep(diffgenes, 1, apply(combined, 1, median, na.rm = TRUE))
print(diffgeneskm[1:5, 1:5])
diffgeneskm <- as.matrix(diffgeneskm)
diffgeneskm = as.dist(diffgeneskm)

#Run Kmeans Clustering
k <- 3
diffgeneskmt = t(diffgenes)
kmeans_result <- kmeans(diffgenest, centers = k)
print(kmeans_result)
gene_clusters <- data.frame(Gene = genes, Cluster = kmeans_result[["cluster"]])

plot(diffgenest, col = kmeans_result$cluster, pch = 16, main = "K-means Clustering of k=3 subgroups by gene", xlab = "X", ylab = "Y")
points(kmeans_result$centers, col = 1:k, pch = 8, cex = 2)  # Plot cluster centers
legend("topright", legend = paste("Cluster", 1:k), col = 1:k, pch = 16) 

#Print top 10 gene names for each cluster
cluster1 = gene_clusters$Gene[gene_clusters$Cluster == 1]
cluster2 = gene_clusters$Gene[gene_clusters$Cluster == 2]
cluster3 = gene_clusters$Gene[gene_clusters$Cluster == 3]
print(cluster1[1:10])
print(cluster2[1:10])
print(cluster3[1:10])
