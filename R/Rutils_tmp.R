## filter and match function ###################################################


filter_match <- function(a,b,a_column=NA,b_column=NA,verbose=F,order_by_a=T){
  
  message("The first argument is a")
  
  if(is.null(a_column) || is.null(b_column)){
    
    stop("The column number cannot be empty")
    
    }
  
  
  
  
  ##  Column Slicing
  
  a_i <- a[,a_column,drop=F]
  clean_a <- get_rid_of_NA(a,a_column)
  b_j <- b[,b_column,drop=F]
  clean_b <- get_rid_of_NA(b,b_column)
  
  ## Subseting columns
  
  a_subset <- clean_a[clean_a[,a_column]  %in% clean_b[,b_column],]
  b_subset <- clean_b[clean_b[,b_column] %in% clean_a[,a_column],]
  
  ## subset column
  a_subset_i <- a_subset[,a_column,drop=T]
  b_subset_i <- b_subset[,b_column,drop=T] 
  
  if(verbose==T){
    
   message("The data set cleaning is done! \n") 
    
  }
  
  ##  matching
  
  if(isTRUE(order_by_a)){                                    # Order by a
  
  b_matched_by_a  <- b_subset[match(a_subset_i,b_subset_i),]
  return(list(b=b_matched_by_a,a=a_subset))
  
  }else{ 
  
    
    
  a_matched_by_b  <- a_subset[match(b_subset_i,a_subset_i),] # Order by b
  return(list(a=a_matched_by_b,b=b_subset))
  
  }
  
  if(verbose==T){
    
    message("The matching is done! \n") 
    
  }
  
  }
  

## Cleaning data ###############################################################



get_rid_of_NA <- function(a,a_column){
  
  
  ## Sanity check
  if(is.null(a)){       stop("The vector can not be null"       )}
  if(is.null(a_column)){stop("The column number can not be null")}
  a_i <- a[,a_column,drop=F]
  if(sum(is.null(a_i))==dim(a)[1]){stop("All values of vector can not be null")}
  
  ## Removing missing value
  ii <- is.null(a_i)
  a_clean <- a[!ii,]

  return(data.frame(a_clean))  
  
}


## Merge Two data sets based on columns ########################################

merge_two_data_sets <- function(a,b,a_column=NA,b_column=NA){

if(is.null(a_column)||is.null(b_column)){stop("The merge columns cannot be empty")}  
    
# a_column to match
# b_column to match
  
matched_res  <- filter_match(a,b,a_column=a_column,b_column=b_column,verbose=F,order_by_a=T)
  
a_matched  <- matched_res[[1]]  
b_matched  <- matched_res[[2]]
merge_df  <- cbind(a_matched,b_matched)

return(merge_df)

}

################################################################################

##  Make an example for data set cleaning

# a <- data.frame(a1=c(1,NA,NA,2,4), a2=c(1,5,NA,6,2))
# b <- data.frame(b1=c(NA,NA,2,5,1),b2=c(1,5,6,2,4))
# 
# merge_two_data_sets(a,b,1,2)

################################################################################
## Preparing results for final


## Add column to results


add_column <- function(res,columns){
  
if(is.null(res)||is.null(columns)){
  
  stop("The arguments can not be empty!")
  
  
}  
  

colnames(res) <- columns

return(res)
  
}

## Get zscore value from estimate and standard deviation

get_Zscore <- function(result,est_col_num,std_col_num){
  

if(is.na(result)||is.na(est_col_num)||is.na(std_col_num)){
    
    stop("The arguments can not be empty!")

  }    
  
  
est_col     <- result[,est_col_num]
std_col_num <- result[,std_col_num] 
Zscore      <- est_col/std_col
  
return(Zscore)

}



## Get coverage from Zscore and annotion coverage value


get_pvalue <- function(Z_score=NULL){
  
if(is.null(Z_score)){
  
  message("The Z_score can not be null!")
} 
  

pval <-   pnorm(Z_score,lower.tail = F)



return(pval)

}


# Get bh adjust

get_bh_adjust <- function(pval){
  
  
  if(is.null(pval)){
    
    message("The pval can not be null!")
  }   
  
  
  Adjusted_pval <- p.adjust(pval,"BH")
  
  return(Adjusted_pval)
  
}


# Adjusted pval coverage

get_adjusted_pval_coverage <- function(Adjusted_pval,coverage){
  
  if(is.null(Adjusted_pval)||is.null(coverage)){
    
    stop("The arguments can not be empty!")
    
  }   
  
  
  coverage_adj_pval     <- Adjusted_pval/coverage
  
  return(coverage_adj_pval)  
}



get_log_pval <- function(pval){
  
  if(is.null(pval)){
    
    stop("The arguments can not be empty!")
    
  } 
  
  log_pval <- -log10(pval)
   
  return(log_pval)
  
}


# write talbe

write_table <- function(obj,path){

write.table(obj,file=path,
            col.names = T,
            row.names = F,
            quote = F)
}


############################################ Graphics #############################################

# http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/#modifying-the-text-of-legend-titles-and-labels

ggplot_group <- function(bar_data,x_column_name,x_column_order_name,y_column_name,fill_column_name,xlab,ylab,title){
  
  bp <- ggplot(bar_data, aes(x = reorder(bar_data[,x_column_name],bar_data[,x_column_order_name]),
                             y = as.numeric(bar_data[,y_column_name]),
                             fill = bar_data[,fill_column_name] ))
  
  
  bp <-  bp + geom_bar(stat="identity",position="dodge")+ ylab(ylab)+ xlab(xlab)
  bp <- bp + theme(axis.text.x=element_text(color = "black", size=7, angle=90, vjust=.8, hjust=0.8))
  bp <- bp + ggtitle(title)+theme(plot.title = element_text(size = 15)) +scale_fill_discrete(name="Ethnic")
  
  return(bp)
  
}


ggplot_group <- function(bar_data,x_column_name,x_column_order_name,y_column_name,fill_column_name,xlab,ylab,title){
  
  bp <- ggplot(bar_data, aes(x = reorder(bar_data[,x_column_name],bar_data[,x_column_order_name]),
                             y = as.numeric(bar_data[,y_column_name]),
                             fill = bar_data[,fill_column_name] ))
  
  
  bp <-  bp + geom_bar(stat="identity",position="dodge")+ ylab(ylab)+ xlab(xlab)
  bp <- bp + theme(axis.text.x=element_text(color = "black", size=7, angle=90, vjust=.8, hjust=0.8))
  bp <- bp + ggtitle(title)+theme(plot.title = element_text(size = 15)) +scale_fill_discrete(name="Ethnic")
  
  return(bp)
  
}

#########################################################################################
### Here I am trying to get the features, the length and cluster number for each X
### The point is do a pc on each correlated features and summarize each cluster into a pc
#########################################################################################

getcorCluster <- function(X,ld_thre){
  clustered_X   <- list()
  cluster_length<- list()
  clust_ii_list <- list()
  ################################################################################################
  distMat       <- as.dist(1-cor(X))
  res           <- hclust(distMat)
  clusters      <- cutree(res,h=1-ld_thre)
  num           <- max(clusters)
  ################################################################################################
  
  for(i in 1:num){
    clust_ii        <- which(clusters==i)
    X_clust_ii      <- X[,clust_ii]  
    clust_ii_list   <- c(clust_ii_list,list(X_clust_ii)) 
    clustered_X     <- c(clustered_X,list(as.matrix(X_clust_ii)))  
    cluster_length  <- c(cluster_length,dim(clustered_X[[i]])[2])
  }
  
  return(list(clustered_X=clustered_X,cluster_length=cluster_length,number_cluster=num))
}

######################## Summarizing each feature and get the X values ##############################





getLDPC <- function(res,thre){
  prX <- NULL
  for(i in 1:length(res$cluster_length)){
    
    
    if(res[["cluster_length"]][[i]]>1){
      
      index <- pc_index(res[["clustered_X"]][[i]],thre = thre)
     # cat(paste0('high_',i,"index_",index,"\n"))
      prRes <- prcomp(res$clustered_X[[i]])
      prx   <- prRes$x[,1:index]
      prX <- cbind(prX,prx)
    }else{
      #cat(paste0('low_',i,"\n"))
      prX <- cbind(prX,X=matrix(res[["clustered_X"]][[i]],ncol=1)) 
      
    }
  }
  return(prX)
}



########################################################################################
### Here I am asking this question on what index the variance is above 0.95
########################################################################################
pc_index <- function(X,thre=0.5){
  # inputs : X
  pc_res <- prcomp(X)
  cum_sum_vec <- cumsum(pc_res$sdev)
  pc_precent<- cum_sum_vec/sum(pc_res$sdev)
  pc_0.05<- which(pc_precent>thre)
  index<- min(pc_0.05)
  return(index)
}


# 
# barplot(height,
#         width = 1,
#         names.arg = NULL,
#         legend.text = NULL,
#         beside = FALSE,
#         horiz = FALSE,
#         angle = 45,
#         col = NULL,
#         border = par("fg"),
#         main = NULL, 
#         sub = NULL, 
#         xlab = NULL, 
#         ylab = NULL,
#         xlim = NULL,
#         ylim = NULL,
#         xpd = TRUE, 
#         axes = TRUE,
#         axisnames = TRUE,
#         inside = TRUE, 
#         plot = TRUE,
#         axis.lty = 0, 
#         offset = 0,
#         add = FALSE,
#         args.legend = NULL, ...){
#   
#   
#   
#   require(reshape2)
#   require(ggplot2)
#   
#   data_1_df  <- read.table(data_set1,header=T,stringsAsFactors = F)
#   data_2_df  <- read.table(data_set2,header=T,stringsAsFactors = F)
#   
#   
#   
#   bar_plot_m <- cbind(data1= data_1_df[,col_number],data2=data_2_df[,col_number])
#   row.names(bar_plot_m) <- data_1_df[,1]
#   bar_plot <- bar_plot_m[order(bar_plot_m[,1],(bar_plot_m[,1]-bar_plot_m[,2]),decreasing=F),]
#   
#   
#   bar_plot <- as.data.frame(cbind(row.names(bar_plot),bar_plot),stringAsfactor=F)
#   colnames(bar_plot) <- data_set_names
#   #ratio <- get_ratio(bar_plot);
#   df.long<-melt(bar_plot,id.vars =col_names[1])
#   df.long$value <- as.numeric(df.long$value)
#   
#   
#   #aes(x=reorder(Cell,-ratio$ratio)  Fix this later
#   
#   plotF <- ggplot(df.long,aes(x=Var1,value,fill=Var2))+ geom_bar(stat="identity",position="dodge") +
#     ylab(ylab) + xlab(xlab) + theme(axis.text.x=element_text(color = "black", size=4, angle=20, vjust=.8, hjust=0.8))+ labs(fill ="Population")+ theme_bw()+ theme(axis.text.x=element_text(angle=90,hjust=1,size=7)) + ggtitle(title)+theme(plot.title = element_text(size = 5)) 
#   
#   
#   return(plotF)
#   
#   
#   
#   
#   
# }

