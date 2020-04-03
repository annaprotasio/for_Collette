# topgo

rm(list=ls())


out.dir = "/Users/ap6/Documents/Google_Drive_avp25/Bioinformatics_coll/Collette/for_Collette/output/"

# libraries ---------------------------------------------------------------

library(topGO)

# function ----------------------------------------------------------------

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


# topGO mbcluster ----------------------------------------------------

# runtopGO.simple(ontol = "BP",genes =  mbc_cls_genes$cls_1, reference = "data/GO_mart_export_topGO_ready.txt")

mbc_cls_genes = readRDS("RData/01_day5genes_cluster_mbc_cls_genes.RData")

mbc_go = list()

for (i in 1:length(mbc_cls_genes)) {
  if (length(mbc_cls_genes[[i]]) > 99) {
    mbc_go[[i]] = runtopGO.simple("BP", mbc_cls_genes[[i]], "data/GO_mart_export_topGO_ready.txt")
  } else {
    mbc_go[[i]] = "gene count less than 100"
  }
}

WriteXLS(lapply(mbc_go, function (x) as.data.frame(x)) 
         , ExcelFileName = paste(out.dir,"tables/mbc/mbc_cls_TopGO.xls",sep="")
         , SheetNames = names(mbc_cls_genes), AdjWidth = T
)


# group MBC clusters ------------------------------------------------------

# group those increased at 5dayPI

mbc.keep = c(1,8,13,20)

day5.mbc = unique(as.character(unlist(mbc_cls_genes[mbc.keep])))

day5.mbc.topGO = runtopGO.simple("BP", day5.mbc, "data/GO_mart_export_topGO_ready.txt")

WriteXLS(day5.mbc.topGO, ExcelFileName = paste(out.dir,"/tables/mbc/day5.mbc.topGO.xls",sep=""), AdjWidth = T)

# topGO kohonen ------------------------------------------------------

koh_cls_genes = readRDS("RData/01_day5genes_cluster_koh_cls_genes.RData")

koh_go = list()

for (i in 1:length(koh_cls_genes)) {
  if (length(koh_cls_genes[[i]]) > 99) {
    koh_go[[i]] = runtopGO.simple("BP", koh_cls_genes[[i]], "data/GO_mart_export_topGO_ready.txt")
  } else {
    koh_go[[i]] = "gene count less than 100"
  }
}

WriteXLS(lapply(koh_go, function (x) as.data.frame(x)) 
         , ExcelFileName = paste(out.dir,"tables/koh/koh_cls_TopGO.xls",sep="")
         , SheetNames = names(koh_cls_genes), AdjWidth = T
)

# group KOH clusters ------------------------------------------------------

# group those increased at 5dayPI

koh.keep = c(1,8,13,20)

day5.koh = unique(as.character(unlist(koh_cls_genes[koh.keep])))

day5.koh.topGO = runtopGO.simple("BP", day5.koh, "data/GO_mart_export_topGO_ready.txt")

WriteXLS(day5.koh.topGO, ExcelFileName = paste(out.dir,"/tables/koh/day5.koh.topGO.xls",sep=""), AdjWidth = T)

