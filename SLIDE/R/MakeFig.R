# ER_res
# Threshold
# Index

SigGenes_A <-function(ER_Res,k,x,y,thre=0.1,negcol="blue",posCol="red",orderbycol=T){
  
  ## Sanity check
  if(is.null(ER_Res$A)){stop("ES_Res and assign matrix can not be empty")}
  if(is.null(k)){stop("K can not be empty")}
  if(is.null(colnames(x))){stop("The row names can not be empty we need it for correlation")}
  if(is.null(y)){stop("The y cannot be empty")}
  if(length(y)!=dim(x)[1]){stop("y and x should be in same sample size! and y should not have more than one dimension")}
  
 
  ## Getting the significant genes
  ER_Res_abs    <- abs(ER_Res$A)
  ii            <- order(abs(ER_Res_abs[,k]),decreasing = T)
  sigGenes      <- which(ER_Res_abs[ii,k]>thre)
  sigGenesload  <- ER_Res_abs[names(sigGenes),k,drop=F]
  colnames(sigGenesload) <- paste0("A",k)
  
  
  ## Getting the correlation
  report <- NULL
  for(i in names(sigGenes)){
  
  report <- rbind(report,data.frame(correlation=cor(x[,i],y,method = "spearman"),sign=sign(cor(x[,i],y,method = "spearman"))))
    
  }  

  report$color <- ifelse(report$sign<0,negcol,posCol)
  corReport <- data.frame(load=sigGenesload,report)
  
  
  if(orderbycol){
    
    corReport <- corReport[order(corReport$color),]
    
  }
  
 return(corReport)
  
}

