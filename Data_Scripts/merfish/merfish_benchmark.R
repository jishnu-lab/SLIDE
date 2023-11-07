############################################################################
############################################################################
###                                                                      ###
###                           SLIDE ON MERFISH                           ###
###                                                                      ###
############################################################################
############################################################################

# https://www.nature.com/articles/s41586-021-03705-x#code-availability


##---------------------------------------------------------------
##                      load merFISH data                      --
##---------------------------------------------------------------
library(anndata)
library(dplyr)

merFISH <- read_h5ad("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/merFISH/counts.h5ad")
#toy <- merFISH[1:100, 1:100]

# load gene names that are through combinatorial sequencing imaging of merFISH
gene_names <- as.matrix(read.table("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/merFISH/genes_combinatorial.txt"))
merFISH <- as.data.frame(merFISH)
idx <- which(colnames(merFISH) %in% gene_names)  # only use the genes that are used in combinatorial sequencing imaging
merFISH <- merFISH[ ,idx]

# MARCH1 is not in the count matrix...
# check <- colnames(merFISH)
# !(gene_names %in% check)
# identical(check, gene_names)
# identical(check[110:120], gene_names[110:120])

#make sure that the orders are the same in the count file and the cell_label file
new_label <- c(1: nrow(merFISH))
merFISH['cell_label'] <- new_label # add the new labels since R can't read more than 22 digits

# load the python processed cell label file
cell_names <- read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/merFISH/cell_labels_V2.csv", row.names = 1)
if (nrow(cell_names) != nrow(merFISH)) {stop( "loaded wrong data...")}

##---------------------------------------------------------------
##               ordinal analysis pre-processing               --
##---------------------------------------------------------------

# filter out useless cells for this analysis
#y <- subset(cell_names, grepl("^L\\d", cell_names$subclass))
y <- subset(cell_names, (grepl("^L\\d", cell_names$subclass)) | (cell_names$subclass %in% c("Lamp5", "Sncg", "Vip", "Sst", "Pvalb")))
dim(y)
subset_cells <- as.matrix(y["new_ident"])

subset_mat <- merFISH[merFISH$cell_label %in% subset_cells, ]
dim(subset_mat)
if (nrow(subset_mat) != length(subset_cells)){stop("wrong dimensions...")}

# encode Y to numerical labels
# the ordinal relationship is extracted from figure 1 panel b
y_num <- as.matrix((as.numeric(recode(as.matrix(y['subclass']), "L2/3 IT" = "0", "L4/5 IT" = "0", "L5 IT" = "0", "L6 IT" = "0", "L6 IT Car3" = "0", "L5 ET" = "0", "L5/6 NP" = "0",
       "L6 CT" = "0", "L6b" = "0", "Lamp5" = "1", "Sncg" = "1", "Vip" = "1", "Sst"= "1", "Pvalb" = "1"))))
y["new_y"] <- y_num[, 1]

#final x
row.names(subset_mat) <- subset_mat$cell_label
subset_mat$cell_label <- NULL

#final y
row.names(y_num) <- y$new_ident

#write.csv(subset_mat, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/merFISH/Catagorical_Analysis/Data/x.csv")
#write.csv(y_num, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/merFISH/Catagorical_Analysis/Data/y.csv")

##---------------------------------------------------------------
##                    Binary Pre-processing                    --
##---------------------------------------------------------------
