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


# get_pvalue <- function(Z_score=NULL){
#   
# if(is.null(Z_score)){
#   
#   message("The Z_score can not be null!")
# } 
#   
# 
# pval <-   pnorm(Z_score,lower.tail = F)
# 
# 
# 
# return(pval)
# 
# }


get_pvalue <- function(z){
  
  pvalue=2 * pnorm( abs(z), lower.tail=FALSE)
  
  return(pvalue)
}


log_get_pvalue <- function(pval){
  
  
  
  return(-log10(pval))
}



# Get bh adjust

get_bh_adjust <- function(pval){
  
  
  if(is.null(pval)){
    
    message("The pval can not be null!")
  }   
  
  
  Adjusted_pval <- p.adjust(pval,"BH")
  
  return(Adjusted_pval)
  
}

##################  Deleting strings from columns##################
delete_string_from_column <- function(strings,data){
  
  col_names <- colnames(data) 
  
  for(i in 1:length(strings)){
    
    col_names <-  sub(strings[i],"",col_names)    
    
  }
  
  colnames(data) <- col_names
  return(data)
}
###################################################################
normalize_mat <- function(data,row_wise=T){
  
  if(row_wise){normalized_data <-  diag(1/rowSums(data))%*%data }
  else{     normalized_data  <- data %*% diag(1/colSums(data))}
  
  row.names(normalized_data) <- row.names(data)
  colnames(normalized_data) <- colnames(data)
  
  return(normalized_data)
}
###################################################################

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

get_similarity <- function(a,row_wise=T,row_names){
  if(row_wise){
    M <- matrix(0,nrow = dim(a)[1],ncol=dim(a)[1])    
    diag(M) <- apply(a,1,function(x){return(1/sqrt(crossprod(x)))})
    sim_mat<- (M%*%tcrossprod(a))%*%M
    
    print(paste0("The rownames are length",length(row_names)))
    row.names(sim_mat) <- row_names
    print(paste0("The row dimensions are",dim(sim_mat)[2]))
    #colnames(sim_mat)  <- row_names
  }
  else{
    a <- t(a)
    M <- matrix(0,nrow = dim(a)[1],ncol=dim(a)[1])    
    diag(M) <- apply(a,1,function(x){return(1/sqrt(crossprod(x)))})
    sim_mat<- (M%*%tcrossprod(a))%*%M
    row.names(sim_mat) <- row_names
    colnames(sim_mat)  <- row_names
  }
  return(sim_mat)
}



get_data_from_dir <-  function(results_path,results_pattern,which_row,which_column,row_names,row_wise,log_scale){
  
  setwd(results_path)
  files      <- list.files()
  ii         <- grep(results_pattern,list.files())
  results    <- files[ii]
  
  
  
  
  res_list            <- sapply(results,function(x){return(as.matrix(read.delim(x,header = T)))},simplify = F)
  res_row_filter_list <- lapply(res_list,function(x){x[1:which_row,]}) ## filter row
  res_col_filter_df   <- sapply(res_row_filter_list,function(x){as.numeric(x[,which_column])},simplify = T) ## filter column
  
  
  
  if(log_scale){
  row.names(res_col_filter_df) <- res_row_filter_list[[1]][,1]
  res_pavl_filter_df <- apply(res_col_filter_df,c(1,2),log_get_pvalue)
  row.names(res_pavl_filter_df) <- row_names
  return(res_pavl_filter_df)}
  else{
    row.names(res_col_filter_df) <- row_names
    return(abs(res_col_filter_df))
      }
}

simulateX <- function(n=24,p=100,k=3,amplitude=0.25){
  
  #### Generating data
  # n = 24          # number of observations
  # p = 100         # number of variables
  # k = 3           # number of variables with nonzero coefficients
  # amplitude = 0.25   # signal amplitude (for noise level = 1)
  # Generate the variables from a multivariate normal distribution
  mu = rep(0,p)
  rho = 0.25
  Sigma = toeplitz(rho^(0:(p-1)))
  Sigma <- (t(Sigma)+Sigma)/2
  X = matrix(rnorm(n*p),n) %*% chol(Sigma)
  
  #### Generate the response from a linear model
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero) / sqrt(n)
  y.sample = function(X) X %*% beta + rnorm(n)
  y = y.sample(X)
  
  return(list(y=y,X=X,original=nonzero))
  
  
  
}


make_positive_definite <- function(Sigma){
  
  eigen_res <- eigen(Sigma)
  new_eigen_value <- ifelse(eigen_res$values<10e-6,10e-4,eigen_res$values)
  Sigma2 <- eigen_res$vectors%*%diag(new_eigen_value)%*%t(eigen_res$vectors)
  Sigma2 <- (Sigma2 +t(Sigma2))/2
  Sigma <- Sigma2 
  return(Sigma2)
  
}

knockoffsG  <- function(x) create.gaussian(X,mu=mu,Sigma = Sigma)
knockoffsF  <- function(x) create.fixed(X,sigma = Sigma)


## Gaussian Knockoffs

get_gaussian_knockoff <- function(X,y,statistic=stat.glmnet_lambdasmax,fdr=0.05){
  
  for(i  in 1:iter){
    
    result         <-     knockoff.filter(X, y, knockoffs = knockoffsG, statistic = statistic,offset=0,fdr=fdr)
    selected_list  <-     c(selected_list,list(c(result$selected)))
    knockoff_list  <-     c(knockoff_list,list(as.matrix(result$Xk)))
    
  }
  
  lenList        <- lapply(selected_list,function(x){length(x)})
  mm             <- which.max(unlist(lenList))
  selected_vars  <- selected_list[[mm]]
  selected_Xk    <- knockoff_list[[mm]]
  
  return(list(selected_vars=selected_vars))
  
}
## Fixed Knockoffs

get_fixed_knok_off <- function(X,statistic=stat.glmnet_lambdasmax,fdr){
  
  result         <-     knockoff.filter(X, y, knockoffs = knockoffsF, statistic = statistic,offset=0,fdr=fdr)
  
  return(list(selected_vars=result$selected))
  
}



#get_snp <- function()



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





































