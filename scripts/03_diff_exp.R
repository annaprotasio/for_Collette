# differential expression

library(DESeq2)

# metadata
sample_map = read.table("/Users/ap6/Documents/Results/Bioinfo_results/Transcriptome/Bpahangi/DATA/sample_lane_matching_sorted")

# define coldata
coldata = data.frame( condition = sub("\\w.(*)","\\1",sample_map$V1)
                      , replicate = sub("(\\w).*","\\1",sample_map$V1)
                      , sample_name = sample_map$V1)
rownames(coldata) = sample_map$V2

# count data
all.data = readRDS("/Users/ap6/Documents/Results/Bioinfo_results/Transcriptome/Bpahangi/DATA/Counts/all_data.RData")
colnames(all.data) = gsub("\\#","_", colnames(all.data))
counts = all.data
counts = counts[,rownames(coldata)]

all(rownames(coldata) == colnames(counts)) # test that all the samples in coldata are in the count table

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "L3")

# keep <- rowSums(counts(dds)) >= 5 # remove low counts
# dds <- dds[keep,]

# QC ----------------------------------------------------------------------

vsd <- vst(dds, blind=FALSE)

# pdf(file = "output/figures/PCA_condition.pdf")
# plotPCA(vsd, intgroup=c("condition"))
# dev.off()


# DE  ---------------------------------------------------------------------

dds <- DESeq(dds)

# res = list(
#   res_RNF145 = results(dds, contrast  = c('condition', 'gRNF145', 'gCTR'))
#   , res_ADIPOR2 = results(dds, contrast  = c('condition', 'gADIPOR2', 'gCTR'))
#   , res_both = results(dds, contrast  = c('condition', 'gRNF145_gADIPOR2', 'gCTR'))
# )

# x24hr vs x5day ----------------------------------------------------------

res_24h_vs_5d = results(dds, contrast  = c('condition', '24hrPI', '5dayPI'), alpha = 0.05)

DESeq2::plotMA(res_24h_vs_5d, alpha = 0.05, main = "res_24h_vs_5d, adj.P-val < 0.05", cex = 1)

up_24h_vs_5d = res_24h_vs_5d[which(res_24h_vs_5d$log2FoldChange > 0 & res_24h_vs_5d$padj < 0.05),]
dw_24h_vs_5d = res_24h_vs_5d[which(res_24h_vs_5d$log2FoldChange < 0 & res_24h_vs_5d$padj < 0.05),]

up_24h_vs_5d.TG = runtopGO.simple("BP", rownames(up_24h_vs_5d), "data/GO_mart_export_topGO_ready.txt")
dw_24h_vs_5d.TG = runtopGO.simple("BP", rownames(dw_24h_vs_5d), "data/GO_mart_export_topGO_ready.txt")

comp_24h_vs_5d = list(
  up_24h_vs_5d_genes = as.data.frame(up_24h_vs_5d)
  , dw_24h_vs_5d_genes = as.data.frame(dw_24h_vs_5d)
  , up_24h_vs_5d.topGO = up_24h_vs_5d.TG
  , dw_24h_vs_5d.topGO = dw_24h_vs_5d.TG
)

WriteXLS(comp_24h_vs_5d, ExcelFileName = "output/tables/diff_exp/comp_24h_vs_5d.xls", row.names = T, AdjWidth = T)
