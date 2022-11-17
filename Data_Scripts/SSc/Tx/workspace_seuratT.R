library(Seurat)
library(dplyr)



## LOAD IN DATA

setwd("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Lafyatis/Pre_Post_Treatment/RawData/")
load("postcompV1.RData")

sample_id <- skinmetadata[skinmetadata$treatment == 'TOFA-W6', ]
sample_id <- sample_id[-c(8, 11), ]
sample_id <- as.character(sample_id$sample)
data <- subset(x = postcomp, subset = orig.ident == sample_id)



## GET TOP 20 HIGH VARIANCE GENES PER CLUSTER


# define function for row variance
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

rowMedian2 <- function(z){
  if (length(median(z[z != 0])) == 0){
    return (0)
  }
  else{
    return(median(z[z!=0]))
  }
}

rowMedian3 <- function(z){
  return(median(z[z!=0]))
}


rowMedian<- function(x){
  #apply(x, 1, median, na.rm=T)
  apply(x, 1, rowMedian3)
}


#rowMedian(tmp_df)
clusters <- c(0:33)
remove <- c( 5, 7, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33 )
clusters <- clusters[!clusters %in% remove]

varList<-c()
dfVar<-NULL


# for (i in 1:length(clusters)){
#   cat("this is cluster ", i-1, "\n")
#   tmp <- subset(x = data, subset = `seurat_clusters` == i-1)
#   print(dim(tmp)[[2]])
#   }
  
# counts_per_patient<- c()
# for (i in 1:length(clusters)){
#   cat("this is cluster ", i-1, "\n")
#   tmp <- subset(x = data, subset = `seurat_clusters` == i-1)
#   print(dim(tmp)[[2]])
#   for (j in 1:length(sample_id)){
#     check <- subset(x = tmp, subset = orig.ident == sample_id[[j]])
#     #print(dim(check))
#     counts_per_patient <- append(counts_per_patient, dim(check)[[2]])
#   }
#   if (min(counts_per_patient) > 1){
#     print(i-1)
#   }
#   #print(dim(tmp))
# }

# for each cluster
for (i in 1:length(clusters)){
  print(i-1)
  tmp <- subset(x = data, subset = `seurat_clusters` == clusters[i])
  print(clusters[i])
  
  # drop clusters that have less than 20 cells in total
  if (dim(tmp)[[2]] > 20){
    
    matS<-matrix(nrow = 23751,ncol=length(sample_id)) # place holder
    
    # for each patient
    for (j in 1:length(sample_id)){
      
      vec<-c()
      tryCatch(
        {
          tmp<-subset(data, idents = clusters[i] ,subset = orig.ident == sample_id[j])
          tmp_df <- as.data.frame(tmp@assays[["SCT"]]@data)
          
          #tmp_df$avg<-rowMeans(tmp_df)
          #tmp_df$avg <-rowMedian(tmp_df)
          #tmp_df$avg[is.na(tmp_df$avg)] <- 0
          #matS[,j]<-tmp_df$avg
          
          ah <- AverageExpression(tmp, return.seurat = TRUE)
          tmp_df<-as.data.frame(ah@assays[["SCT"]]@data)
          matS[,j]<-tmp_df[, 1]
        
          },
        error= function(cond){
          matS[,j]<-c(rep.int(0,nrow(tmp_df)))
        },
        finally = print("done")
      )
    }
    dfBig_H<-as.data.frame(matS)
    dfBig_H[is.na(dfBig_H)] <- 0
    
    # repeat gene :/
    dfBig_H<-dfBig_H[-22664,]
    
    names<-c()
    clstrs<-substr(clusters[i], 1,2)
    names<-paste(clstrs, rownames(tmp_df),sep=".")
    names<-gsub(" ", "", names)
    names<-gsub("-", "", names)
    names<-gsub("/", "", names)
    names<- paste("c", names, sep=".")
    names<-names[-22664]
    rownames(dfBig_H)<-names
    
    # get variation
    dfBig_H$rowVar <- RowVar(dfBig_H)
    dfBig_H <- dfBig_H[order(dfBig_H$rowVar,decreasing=T),]
    
    varOrder <- head(dfBig_H,50)
    if (i==1){
      dfVar<-varOrder
    } else {
      dfVar<- rbind(dfVar, varOrder)
    }
    varList<-c(varList, row.names(varOrder))
    
  }

}

Var.df<-as.data.frame(t(dfVar[,-length(dfVar)]))
rownames(Var.df) <- sample_id
write.csv(Var.df, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Lafyatis/Pre_Post_Treatment/HER_080422/Data/avg_SCT_X.csv")

Y <- c(1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1)
Y <- as.data.frame(Y)
rownames(Y) <- sample_id
#write.csv(Y, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Lafyatis/Pre_Post_Treatment/HER_080322/Data/Y.csv")
