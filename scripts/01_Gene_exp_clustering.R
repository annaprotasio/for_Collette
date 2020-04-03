#

rm(list=ls())

# mbclusterseq ------------------------------------------------------------

load("/Users/ap6/Documents/Results/Bioinfo_results/Transcriptome/Bpahangi/ANALYSIS/MBCluster/mbcluster.RData")

# clusters that present incresed exp at day5
# cls = 8, 20, 1 and 13

keep = c(1:20)

mbc_cls_genes = list()

for (i in 1:length(keep)) {
  mbc_cls_genes[[i]] = rownames(mydata$logFC[which(cls$cluster == keep[i] & VeryGoodP),])
}

names(mbc_cls_genes) = paste("cls_",keep,sep="")

h.dir = "/Users/ap6/Documents/Google_Drive_avp25/Bioinformatics_coll/Collette/for_Collette/output/"

WriteXLS(lapply(mbc_cls_genes, function (x) as.data.frame(x)) 
         , ExcelFileName = paste(h.dir,"tables/mbc/mbc_exp_clusters.xls",sep="")
         , SheetNames = names(mbc_cls_genes), AdjWidth = T
         )

saveRDS(mbc_cls_genes,"/Users/ap6/Documents/Google_Drive_avp25/Bioinformatics_coll/Collette/for_Collette/RData/01_mbc_cls.RData")


# kohonen clustering ------------------------------------------------------


load("DESeq_ANALYSIS/Kohonen_Bpah_larva.RData")

# List of cluster number, the ordering is the same as input gene table
# Cluster numbers in plot go like this, starting from bottom left
# ...
# 1 2 3 4 5 
# 6 7 8 9 10 
# 11 12 13 14 15 

# all 40 clusters
keep = c(1:40)

koh_cls_genes=list()
for (i in 1:length(keep)) {
  koh_cls_genes[[i]]=row.names(rld_larva_meanNorm)[which(data_model$unit.classif %in% c(i))]
}

names(koh_cls_genes) = paste("cls_",keep,sep="")

WriteXLS(lapply(koh_cls_genes, function (x) as.data.frame(x)) 
         , ExcelFileName = paste(h.dir,"tables/koh/koh_exp_clusters.xls",sep="")
         , SheetNames = names(koh_cls_genes), AdjWidth = T
)


saveRDS(koh_cls_genes,"/Users/ap6/Documents/Google_Drive_avp25/Bioinformatics_coll/Collette/for_Collette/RData/01_koh_cls.RData")
