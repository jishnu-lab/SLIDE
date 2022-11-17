rm(list = ls())
cat("\014")
library(Seurat)
library(ggplot2)

##################################################################
##                Load All Seurat Objs to a List                ##
##################################################################

list_dir <- c('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/SpaceRanger/Poholek_A1/outs',
              '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/SpaceRanger/Poholek_B1/outs',
              '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/SpaceRanger/Poholek_C1/outs')


LoadData <- function(list_dir){
  list_seurat <- NULL
  for (i in 1:length(list_dir)){
    data <- Seurat::Load10X_Spatial(
      # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
      data.dir = list_dir[i], 
      filename = 'filtered_feature_bc_matrix.h5',
      assay = "Spatial", # specify name of the initial assay
      slice = "slice1", # specify name of the stored image
      filter.matrix = TRUE, 
      to.upper = FALSE)
    list_seurat <- append(list_seurat, data)
  }
  return(list_seurat)
}

list_seurat <- LoadData(list_dir)


#################################################################
##           Go Through Each Objs to Define Clusters           ##
#################################################################
#map the name of the seurat into index in the list of seurats
name2obj <- list()
keys <- c("A1", "B1", "C1")
vals <- c(1, 2, 3)
name2obj[keys] <- vals

name <- 'A1'
name2obj[[name]]
data <- list_seurat[[name2obj[[name]]]]


data$orig.ident <- name
data$project.name <- name


data <- PercentageFeatureSet(data, "^mt-", col.name = "percent.mito")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
#does the hemo globin gene need to be removed?
#brain <- PercentageFeatureSet(brain, "^Hb.*-", col.name = "percent_hb")

VlnPlot(data, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mito"), pt.size = 0.1, ncol = 2) + NoLegend()

data <- subset(data, subset = nFeature_Spatial < 6000 & nFeature_Spatial > 3000 & 
    nCount_Spatial < 50000 & percent.mt < 10)

SpatialFeaturePlot(
  data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "bottom")

data <- SCTransform(data, assay = "Spatial", verbose = FALSE)

data <- RunPCA(data, assay = "SCT", verbose = FALSE)
data <- FindNeighbors(data, reduction = "pca", dims = 1:30, k.param = 20) # SNN graph
data <- FindClusters(data, verbose = FALSE) # 
data <- RunUMAP(data, reduction = "pca", dims = 1:30)

plot3 <- DimPlot(data, reduction = "umap", label = TRUE) + NoLegend()
plot4 <- SpatialDimPlot(data, label = TRUE, label.size = 3) + NoLegend()
plot3 
plot4

SpatialDimPlot(data, cells.highlight = CellsByIdentities(object = data, idents = c(6)), facet.highlight = TRUE)

data <- FindSpatiallyVariableFeatures(data, assay = "SCT", features = VariableFeatures(data)[1:1000],
                                       selection.method = "markvariogram")

top.features <- head(SpatiallyVariableFeatures(data, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(data, features = top.features, ncol = 3, alpha = c(0.1, 1))

#saveRDS(data, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/Seurat_Analysis_Results&Plots/PostQC_A1.RDS')


##################################################################
##                  Pool Regions From Clusters                  ##
##################################################################
path =c('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/Spatial_ER_050422/Seurat_Analysis_Results&Plots/PostQC_A1.RDS',
        '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/Spatial_ER_050422/Seurat_Analysis_Results&Plots/PostQC_B1.RDS',
        '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/Spatial_ER_050422/Seurat_Analysis_Results&Plots/PostQC_C1.RDS')

list_seurat <- NULL
for (i in 1:length(path)){
  list_seurat <- append(list_seurat, readRDS(path[i]))
}

GetSameCells <- function(list_seurat, pool_idents){
  all_idents <- vector("list", length = length(pool_idents))
  for (i in 1:length(list_seurat)){
    data <- list_seurat[[i]]
    cell_ident <- row.names(data@meta.data)
    cell_ident <- cell_ident[which(data@meta.data$seurat_clusters %in% pool_idents[[i]])]
    all_idents[[i]] <- cell_ident
  }
  return(all_idents)
}

GetIntersectGenes <- function(list_seurat){
  common_genes <- intersect(intersect(list_seurat[[1]]@assays$SCT@data@Dimnames[[1]],
                                      list_seurat[[2]]@assays$SCT@data@Dimnames[[1]]),
                            list_seurat[[3]]@assays$SCT@data@Dimnames[[1]])
  return(common_genes)
}

GetCountMatrices <- function(list_seurat, common_genes, all_idents){
  counts <- vector('list', length = 3)
  for (i in 1:length(list_seurat)){
    data <- list_seurat[[i]]
    mat <- subset(x = data, features = common_genes, cells = all_idents[[i]])
    counts[[i]] <- data.frame(mat@assays$SCT@data)
  }
  return(counts)
}



#each element of the vector is cluster identity of a seurat object. 
#cluster 3 and 4 are B cells in A1.
BC_idents <- list(c(3, 4), c(1, 2, 3, 5), c(5,6))
DC_idents <- list(c(0, 1, 2), c(0), c(0, 1, 2))
#int_idents <- list(c(5, 6), c(4, 6), c(4))

# get the UMI for the same cell type from different objects
all_idents <- GetSameCells(list_seurat, BC_idents)
# get the intersected gene names
common_genes <- GetIntersectGenes(list_seurat)
# get the count matrices 
counts <- GetCountMatrices(list_seurat, common_genes, all_idents)
colnames(counts[[1]]) <- paste('A1', colnames(counts[[1]]), sep = "-")
colnames(counts[[2]]) <- paste('B1', colnames(counts[[2]]), sep = "-")
colnames(counts[[3]]) <- paste('C1', colnames(counts[[3]]), sep = "-")
final_BC <- rbind(t(counts[[1]]),
               t(counts[[2]]),
               t(counts[[3]]))



all_idents <- GetSameCells(list_seurat, DC_idents)
common_genes <- GetIntersectGenes(list_seurat)
counts <- GetCountMatrices(list_seurat, common_genes, all_idents)
colnames(counts[[1]]) <- paste('A1', colnames(counts[[1]]), sep = "-")
colnames(counts[[2]]) <- paste('B1', colnames(counts[[2]]), sep = "-")
colnames(counts[[3]]) <- paste('C1', colnames(counts[[3]]), sep = "-")
final_DC <- rbind(t(counts[[1]]),
                       t(counts[[2]]),
                       t(counts[[3]]))

# all_idents <- GetSameCells(list_seurat, int_idents)
# common_genes <- GetIntersectGenes(list_seurat)
# counts <- GetCountMatrices(list_seurat, common_genes, all_idents)
# colnames(counts[[1]]) <- paste('A1', colnames(counts[[1]]), sep = "-")
# colnames(counts[[2]]) <- paste('B1', colnames(counts[[2]]), sep = "-")
# colnames(counts[[3]]) <- paste('C1', colnames(counts[[3]]), sep = "-")
# final_int <- rbind(t(counts[[1]]),
#                   t(counts[[2]]),
#                   t(counts[[3]]))



#################################################################
##           concatenate all Cell type and produce Y           ##
#################################################################
rownames(final_BC) <- paste('BC', rownames(final_BC), sep = "-")
rownames(final_DC) <- paste('DC', rownames(final_DC), sep = "-")
#rownames(final_int) <- paste('INT', rownames(final_int), sep = "-")

#ER_X <- rbind(final_BC, final_DC, final_int)
ER_X <- rbind(final_BC, final_DC)
print(dim(ER_X))
print(dim(final_BC))
print(dim(final_DC))
print(dim(final_int))
Y<-rep(c(0,1),times=c(dim(final_BC)[[1]], dim(final_DC)[[1]]))

write.csv(ER_X, '/Users/xiaoh/Desktop/Research/Project_Scripts/Poholek/Poholek_Spatial/Spatial_ER_061722/Scripts/Concat.csv')
write.csv(Y, '/Users/xiaoh/Desktop/Research/Project_Scripts/Poholek/Poholek_Spatial/Spatial_ER_061722/Scripts/Y.csv')


