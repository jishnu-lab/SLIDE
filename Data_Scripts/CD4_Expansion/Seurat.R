library(Seurat)
library(dplyr)
#read RDS
data <- readRDS('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Raw Data/AJ42_cd4.RDS')

#extract the count matrix
count <- t(data.frame(as.matrix(data@assays[["RNA"]]@data)))
y <- data@meta.data[["cloneType"]]

# delete the cells that doesn't have a label
idx <- which(is.na(y))
y <- y[-idx]
y <- as.data.frame(as.numeric(as.character(recode(y, "Single (0 < X <= 1)" = "0", "Small (1 < X <= 9)" = "1","Medium (9 < X <= 100)" = "2"))))
count <- count[-idx, ]
if (dim(count)[[1]] != dim(y)[1]) {stop("dimension is wrong...\n")}
rownames(y) <- row.names(count)

ZeroFiltering <- function(data, g_thresh, c_thresh){
  cat("Original dataframe dimension is ", dim(data)[[1]], " by ", dim(data)[[2]], "\n")

  #filtered out cells with num zeros greater than a certain threshold
  thresh <- dim(data)[[2]] - g_thresh
  i <- rowSums(data == 0, na.rm=TRUE) < thresh
  filtered <- data[i, ]

  #filtered out genes with num zeros greater than a certain threshold
  thresh <- dim(data)[[1]] - c_thresh
  i <- colSums(data == 0, na.rm=TRUE) < thresh
  filtered <- filtered[ ,i]

  cat("Filtered dataframe dimension is ", nrow(filtered), " by ", ncol(filtered), "\n")
  return(filtered)
}


filtered_count <- ZeroFiltering(count, 1200, 1200)
filtered_y <- y %>% filter(row.names(y) %in% row.names(filtered_count))
if (dim(filtered_count)[[1]] != dim(filtered_y)[1]) {stop("dimension is wrong...\n")}


#write to csv
write.csv(filtered_count, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/Data/expansion_x.csv')
write.csv(filtered_y, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/Data/expansion_y.csv')

