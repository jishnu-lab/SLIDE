#library(Seurat)
library(ggplot2)

# replicate and sort the feature vecotrs, the x axis of the dot plot
GetRepFeatures <- function(vec_features, rep_num){
  for (i in 1: length(vec_features)){
    tmp_vec <- rep(vec_features[i], rep_num)
    if (i == 1) {final_vec = tmp_vec} else { final_vec <- c(final_vec, tmp_vec)}
  }
  return(final_vec)
}

# The function get the DF for the 2 row dot plot (low, medium, hi)
FillDataFrame <- function(df1, features, list_CD4){
  # for each feature
  for (j in 1:length(features)){
    gene <- features[[j]]
    print(gene)
    list_per <- NULL
    # for that one gene, loop through all stages
    for (i in 1:length(list_CD4)){
      obj <- list_CD4[[i]]
      count <- obj@assays$RNA@data[gene, ]
      med_gene <- median(count[count > 0 ])  # get the count
      list_per <- append(list_per, length(count[count > 0 ])) # get the count of nonzero cells for that one stage for that one gene
      idx <- ((j - 1) * length(list_CD4)) + i
      df1["count"][idx, ] <- med_gene
    }
    if (length(list_per) != length(list_CD4)){stop("percentage list length is off ... ")}
    norm_per <- (list_per / sum(list_per)) * 100 # normalize the percentage. What we want is in the non-zero cells, what are the percentage of each stage.
    for (k in 1:length(norm_per)){
      idx <- ((j - 1) * length(list_CD4)) + k
      df1["percent"][idx, ] <- norm_per[[k]]
    }
  }
  return(df1)
}


Cast2DF <- function(sig_z, index){
  col <-sig_z[ ,i][idx_small]
  tmp <- as.data.frame(col)
  tmp["group"] <- c(rep("small", length(sig_z[ ,i][idx_small])))
  return(tmp)
}

##==================
##  load objects   =
##==================
CD4_obj <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Raw Data/AJ42_cd4.RDS")

#import the CD4 object from Alok, and subset it since there are NAs in the cloneType
CD4_obj <- subset(x = CD4_obj,
                  subset = (cloneType == "Single (0 < X <= 1)" | cloneType == "Small (1 < X <= 9)" | cloneType == "Medium (9 < X <= 100)"))

DefaultAssay(object = CD4_obj) <- "RNA"

#split the object into 3 different objects
CD4_single <- subset(x = CD4_obj, subset = cloneType == "Single (0 < X <= 1)")
CD4_small <- subset(x = CD4_obj, subset = cloneType == "Small (1 < X <= 9)")
CD4_medium <- subset(x = CD4_obj, subset = cloneType == "Medium (9 < X <= 100)")

#################################################################
##                3 row plot (low, miedum,  hi)                ##
#################################################################
list_CD4 <- c(CD4_single, CD4_small, CD4_medium)
# the genes we are interested in plotting
# features <- c("Lag3", "Tigit", "Sell", "Ccr7", "Cd200", "Pdcd1", "Ndfip1", "Tbc1d4", "Cd82",
#               "Ptpn11", "Prr13", "Mllt3", "Ankrd12", "Gimap5", "Anxa5")

features <- c("Lag3", "Tigit", "Pdcd1", "Cd200", "Sell", "Ccr7")

# The y
stages <- c("Low", "Medium", "High")

r_features <- GetRepFeatures(features, length(stages))
r_stages <- rep(stages , length(features))
# sanity check for the feature names
for ( f in features){
  if(f %in% row.names(CD4_obj) == FALSE){stop("wrong feature names")}
}

df1 <- as.data.frame(r_features)
df1["r_stages"] <- r_stages
df1["count"] <- rep(0, length(r_stages))
df1["percent"] <- rep(0, length(r_stages))

df1 <- FillDataFrame(df1, features, list_CD4)
write.csv(df1, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/5_h.csv")
y_order <- stages
x_order <- features
# now move on to plotting df1
ggplot(data = df1, aes(x=factor(r_features, level = x_order), y = factor(r_stages, level = y_order), color = -count, size = percent)) +
  geom_point()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))


##################################################################################

nod_combined <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/Unanue_JEM/Combined_Obj.RDS")

nod_4w_inter <- subset(x = nod_combined, subset = (orig.ident == "nod-4w-1" | orig.ident == "nod-4w-2"))
nod_8w_inter <- subset(x = nod_combined, subset = (orig.ident == "nod-8w-1" | orig.ident == "nod-8w-2"))
nod_15w_inter <- subset(x = nod_combined, subset = (orig.ident == "nod-15w-1"))

nod_list_inter <- list(nod_4w_inter, nod_8w_inter, nod_15w_inter)
# selected_genes <- c("Lag3", "Tigit", "Sell", "Ccr7", "Cd200", "Pdcd1", "Tbc1d4", "Cd82",
#               "Ptpn11", "Gimap5")

selected_genes <- c("Lag3", "Tigit", "Pdcd1", "Cd200", "Sell", "Ccr7")

stages <- c("4w", "8w", "15w")

r_features <- GetRepFeatures(selected_genes, 3)
r_stages <- rep(stages , length(selected_genes))
df1 <- as.data.frame(r_features)
df1["r_stages"] <- r_stages
df1["count"] <- rep(0, length(r_stages))
df1["percent"] <- rep(0, length(r_stages))

df1 <- FillDataFrame(df1, selected_genes, nod_list_inter)
write.csv(df1, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/5_k.csv")
y_order <- stages
x_order <- selected_genes
ggplot(data = df1, aes(x=factor(r_features, level = x_order), y = factor(r_stages, level = y_order), color = -count, size = percent)) +
  geom_point()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))


