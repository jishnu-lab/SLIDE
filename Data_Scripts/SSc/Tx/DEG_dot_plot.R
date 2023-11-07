# plot the expression of the non-linear genes with MRSS

library(Seurat)

load("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/RawData/harmonskinV6.RData")

cluster_name <- mixedsort(as.character(unique(harmonskin@active.ident)))
data <- subset(x = harmonskin, idents = cluster_name[1:18])

cd_genes <- c("KRT1", "KRT2", "KRT10", "KRT14", "KRT15", "S100A2", "ACKR1", "FABP4", 
              "SELE", "SFRP2", "COL1A1", "COMP", "ACTA2", "RERGL", "TAGLN", "RGS5", 
              "RGS16", "CCL2", "CXCL12", "CCL19", "PLA2G2A", "CCL3", "RNASE1", "CXCL8", 
              "COCH", "ASPN", "MFAP4", "LYZ", "HLA-DRA", "G0S2", "CXCR4", "CD52", "LTB", 
              "HIST1H4C", "KRT5", "KRT14", "TPSAB1", "TPSB2", "CTSG", "DES", "ACTG2", 
              "MYLK", "S100B", "GPM6B", "PLP1", "AQP5", "KRT19", "SNORC", "KRT17", 
              "KRT6B", "S100A2", "CCL5",  "NKG7")

cd_genes <- unique(cd_genes)

DotPlot(object = data, features = cd_genes) + 
  RotatedAxis()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))


 d