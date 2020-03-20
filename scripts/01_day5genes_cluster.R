#

rm(list=ls())

setwd("/Users/ap6/Documents/Results/Bioinfo_results/Transcriptome/Bpahangi/")


# day5genes ---------------------------------------------------------------

# mbclusterseq

load("ANALYSIS/MBCluster/mbcluster.RData")

# clusters that present incresed exp at day5
# cls = 8, 20, 1 and 13

keep = c(1,8,13,20)

cls_ls = list()

for (i in 1:length(keep)) {
  cls_ls[[i]] = rownames(mydata$logFC[which(cls$cluster == keep[i] & VeryGoodP),])
}

names(cls_ls) = paste("cls_",keep,sep="")

h.dir = "/Users/ap6/Documents/Google_Drive_avp25/Bioinformatics_coll/Collette/Rproj/"

for (i in 1:length(cls_ls)) {
  write.csv(cls_ls[i]
            , file = paste(h.dir,"tables/mbc/mbc_day5genes_",names(cls_ls[i]),".csv",sep="")
            , quote = F
            , row.names = F
            )
}

# kohonen

load("DESeq_ANALYSIS/Kohonen_Bpah_larva.RData")

# List of cluster number, the ordering is the same as input gene table
# Cluster numbers in plot go like this, starting from bottom left
# ...
# 11 12 13 14 15 
# 6 7 8 9 10 
# 1 2 3 4 5 

# all 40 clusters
keep = c(1,25,33,40)

cluster.res=list()
for (i in 1:length(keep)) {
  cluster.res[[i]]=row.names(rld_larva_meanNorm)[which(data_model$unit.classif %in% c(i))]
}

names(cluster.res) = paste("cls_",keep,sep="")

for (i in 1:length(cluster.res)) {
  write.csv(cluster.res[i]
            , file = paste(h.dir,"tables/koh/koh_day5genes_",names(cluster.res[i]),".csv",sep="")
            , quote = F
            , row.names = F
  )
}


