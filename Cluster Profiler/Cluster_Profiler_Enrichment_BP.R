lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library("rjson")
setwd("~/Documents/Courses /2/2.1/BI/Cluster Profiler")

#--------------------------------------------
genes_to_enr <- c('ROGDI',
 'GABRB3',
 'TSHZ3',
 'CNTNAP2',
 'DNMT3B',
 'DPP6',
 'SNTG2',
 'DHCR24',
 'USH2A',
 'NTSR1',
 'TRIP12',
 'SOX5',
 'CDH10',
 'GABRQ',
 'CIRBP',
 'RYR2',
 'FOXP1',
 'HFE',
 'CHD2',
 'SLC1A1',
 'SOX9',
 'SLC16A2',
 'BDNF',
 'DIO3',
 'DLX1',
 'AVPR1A',
 'MTNR1A',
 'EPHX2',
 'CIC',
 'FEV',
 'TET1',
 'NEMF',
 'DNMT3A',
 'CTTNBP2',
 'CHD8',
 'MYT1L',
 'NRXN1',
 'DDHD1',
 'RAI1',
 'NSUN2',
 'DIP2A',
 'MEF2C',
 'LRRTM3',
 'PCDH9',
 'UNC80',
 'GLRA2',
 'TBL1X',
 'JARID2',
 'CA2',
 'DLG4',
 'CEP41',
 'HEY1',
 'RELN',
 'ANKRD11',
 'MTNR1B',
 'NRXN2',
 'IL1RAPL1',
 'RARB',
 'GTF2I',
 'LAMC3',
 'GPR50',
 'NRXN3',
 'SFSWAP',
 'LRRN3',
 'SIN3A',
 'IQGAP3',
 'DIO2')
#----------- Convert Genes into Symbols for the enrichment
library(org.Hs.eg.db)

gene_to_ID <- function(my.symbols){ 
  
  hs <- org.Hs.eg.db
  symbol_to_id <- select(hs, 
                         keys = my.symbols,
                         columns = c("ENTREZID", "SYMBOL"),
                         keytype = "SYMBOL")
  
  x2 <- na.omit(symbol_to_id$ENTREZID)
  return(x2)
}


gene_ID        <- gene_to_ID(genes_to_enr)

library(clusterProfiler)
library(ggplot2)

organism = 'org.Hs.eg.db'
setwd("~/Documents/Courses /2/2.1/BI/Cluster Profiler/outputs")

dir.create("BP")
ggo <- groupGO(gene = gene_ID,OrgDb = organism,ont="BP", level=3, readable=TRUE)
ego <- enrichGO(gene = gene_ID,OrgDb = organism, ont="BP", pAdjustMethod="BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable=TRUE)

write.csv(as.data.frame(ggo),"BP/ggo_BP.csv")
write.csv(as.data.frame(ego),"BP/ego_BP.csv")

png("BP/barplot_BP_ggo.png",width = 780, height = 480)
barplot(ggo, drop=TRUE, showCategory=12)
dev.off()

png("BP/barplot_BP_ego.png",width = 780, height = 480)
barplot(ego, showCategory=8)
dev.off()

png("BP/dotplot_BP.png",width = 780, height = 480)
dotplot(ego)
dev.off()

png("BP/cnet_plot_BP.png",width = 780, height = 480)
cnetplot(ego, categorySize="pvalue")
dev.off()
