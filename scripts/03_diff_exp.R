# differential expression

# libs
library(DESeq2)
library(topGO)
library(WriteXLS)

# functions ---------------------------------------------------------------

runtopGO.simple = function(ontol, genes, reference) { 
  
  #ontol-> "BP", "CC" or "MF", genes is vector, "reference" is path to file
  
  ref = read.table(file=reference, sep=" ", stringsAsFactor=F)
  
  ### COMPARISON 
  names(ref) = c('id', 'go')
  ref.vec = strsplit(ref$go, split=',', fixed=T)
  names(ref.vec) <- ref$id
  all.ids <- ref$id
  
  
  #for up reg genes
  scores <- rep(0, nrow(ref)) # list of scores
  names(scores) <- ref$id
  scores[ref$id %in% genes] <- 1
  
  # Just selects items w/ a score of 1
  geneSelectionFun <- function(score){
    return(score >= 1)
  }
  
  GOdata <- new("topGOdata",
                ontology = ontol,
                allGenes = scores,
                annot = annFUN.gene2GO,
                gene2GO = ref.vec,
                geneSelectionFun = geneSelectionFun,
                nodeSize = 3, # can change this!!
                description = ''
  )
  
  resultFisher <- runTest(GOdata,algorithm="classic",statistic="Fisher")
  resultTopgo <- runTest(GOdata,algorithm="weight01",statistic="Fisher")
  resultElim <- runTest(GOdata,algorithm="elim",statistic="Fisher")
  
  resultFis=resultFisher
  resultWeight=resultTopgo
  
  
  result <- GenTable(
    GOdata,
    #fisher = resultFisher,
    topGO = resultTopgo,
    #elim = resultElim,
    orderBy = "topGO",
    ranksOf = "fisher",
    topNodes = 100
  )
  
  colnames(result)=c("GO.id","Term","Ann","Sig","Exp","p.value")
  mode(result$p.value)=c("numeric")
  result=subset(result, p.value < 5e-02)
  return(result)
} 



# metadata
sample_map = read.table("data/sample_lane_matching_sorted")

# define coldata
coldata = data.frame( condition = sub("\\w.(*)","\\1",sample_map$V1)
                      , replicate = sub("(\\w).*","\\1",sample_map$V1)
                      , sample_name = sample_map$V1)
rownames(coldata) = sample_map$V2

# count data
all.data = readRDS("RData/all_data.RData")
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

# x5d vs x10d ----------------------------------------------------------

res_5d_vs_10d = results(dds, contrast  = c('condition', '5dayPI', '10dayPI'), alpha = 0.05)

DESeq2::plotMA(res_5d_vs_10d, alpha = 0.05, main = "res_5d_vs_10d, adj.P-val < 0.05", cex = 1)

up_5d_vs_10d = res_5d_vs_10d[which(res_5d_vs_10d$log2FoldChange > 0 & res_5d_vs_10d$padj < 0.05),]
dw_5d_vs_10d = res_5d_vs_10d[which(res_5d_vs_10d$log2FoldChange < 0 & res_5d_vs_10d$padj < 0.05),]

up_5d_vs_10d.TG = runtopGO.simple("BP", rownames(up_5d_vs_10d), "data/GO_mart_export_topGO_ready.txt")
dw_5d_vs_10d.TG = runtopGO.simple("BP", rownames(dw_5d_vs_10d), "data/GO_mart_export_topGO_ready.txt")

comp_5d_vs_10d = list(
  up_5d_vs_10d_genes = as.data.frame(up_5d_vs_10d)
  , dw_5d_vs_10d_genes = as.data.frame(dw_5d_vs_10d)
  , up_5d_vs_10d.topGO = up_5d_vs_10d.TG
  , dw_5d_vs_10d.topGO = dw_5d_vs_10d.TG
)

WriteXLS(comp_5d_vs_10d, ExcelFileName = "output/tables/diff_exp/comp_5d_vs_10d.xls", row.names = T, AdjWidth = T)

