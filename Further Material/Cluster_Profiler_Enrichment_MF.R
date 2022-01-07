lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library("rjson")
setwd("~/Documents/Courses /2/2.1/BI/Cluster Profiler")

genes_to_enr <- c('ABLIM1',
                  'ADRA1D',
                  'ADRB1',
                  'ANKRD26',
                  'APBA1',
                  'ARHGEF26',
                  'BAI1',
                  'C11orf52',
                  'C15orf59',
                  'CASK',
                  'CAV1',
                  'CDC42BPA',
                  'CRHR1',
                  'CTNNA1',
                  'CTNNAL1',
                  'CXADR',
                  'DIRAS3',
                  'DLG1',
                  'DLG2',
                  'DLG3',
                  'DLG5',
                  'DLGAP1',
                  'DLGAP2',
                  'DLGAP3',
                  'DLGAP4',
                  'DMD',
                  'DNAJC5',
                  'DTNA',
                  'DTNB',
                  'DUSP10',
                  'EFR3A',
                  'EFR3B',
                  'EPB41',
                  'EPB41L4A',
                  'ERBB2IP',
                  'ERBB4',
                  'F8A1',
                  'FAM171A1',
                  'FAM171A2',
                  'FAM171B',
                  'FLOT1',
                  'FRS2',
                  'GAB1',
                  'GJA1',
                  'GPRIN3',
                  'GRIN1',
                  'GRIN2A',
                  'GRIN2C',
                  'GRIN3A',
                  'GUCY1A2',
                  'INADL',
                  'KCNA4',
                  'KCNJ12',
                  'KCNJ4',
                  'KHDRBS1',
                  'KIAA0754',
                  'KIF26B',
                  'KRAS',
                  'LCK',
                  'LIN7A',
                  'LIN7B',
                  'LIN7C',
                  'LLGL1',
                  'LPHN2',
                  'LYN',
                  'MAGI2',
                  'MARCKS',
                  'MARK2',
                  'MARK3',
                  'MLLT4',
                  'MPDZ',
                  'MPP2',
                  'MPP3',
                  'MPP5',
                  'MPP6',
                  'MPP7',
                  'MTMR2',
                  'NETO1',
                  'NLGN1',
                  'NLGN2',
                  'NOS1',
                  'OCLN',
                  'PALM',
                  'PARD3',
                  'PHACTR4',
                  'PKP4',
                  'PLCH1',
                  'PLEKHA1',
                  'PLEKHA2',
                  'PSD3',
                  'PTPN13',
                  'PXDC1',
                  'RAB35',
                  'RHOB',
                  'RHPN1',
                  'RPS6KA1',
                  'SCRIB',
                  'SHANK2',
                  'SNTA1',
                  'SNTB1',
                  'SNTB2',
                  'SNTG1',
                  'SPZ1',
                  'STX7',
                  'TAGAP',
                  'TENC1',
                  'TNS1',
                  'TNS3',
                  'USP6NL',
                  'UTRN',
                  'WWC1',
                  'ZBTB7A',
                  'ZDHHC5',
                  'ZFPL1',
                  'ZGPAT')

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

dir.create("MF")
ggo <- groupGO(gene = gene_ID,OrgDb = organism,ont="MF", level=3, readable=TRUE)
ego <- enrichGO(gene =gene_ID,OrgDb = organism, ont="MF", pAdjustMethod="BH", pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable=TRUE)

write.csv(as.data.frame(ggo),"MF/ggo_MF.csv")
write.csv(as.data.frame(ego),"MF/ego_MF.csv")

png("MF/barplot_MF_ggo.png",width = 780, height = 480)
barplot(ggo, drop=TRUE, showCategory=12)
dev.off()

png("MF/barplot_MF_ego.png",width = 780, height = 480)
barplot(ego, showCategory=8)
dev.off()

png("MF/dotplot_MF.png",width = 780, height = 480)
dotplot(ego)
dev.off()

png("MF/cnnet_plot_MF.png",width = 780, height = 480)
cnetplot(ego, categorySize="pvalue")
dev.off()

