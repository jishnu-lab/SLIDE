# ER_res
# Threshold
# Index

SigGenes_A <-function(ER_Res,
                      k,
                      x,
                      y,
                      thre=0.1,
                      negcol="blue",
                      posCol="red",
                      orderbycol=T,
                      filterTop=NULL,
                      baseonCorrlation=F,
                      basedonLoad=F){
  
    ## Sanity check
    if(is.null(ER_Res$A)){stop("ES_Res and assign matrix can not be empty")}
    if(is.null(k)){stop("K can not be empty")}
    if(is.null(colnames(x))){stop("The col names can not be empty we need it for correlation")}
    if(is.null(y)){stop("The y cannot be empty")}
    if(length(y)!=dim(x)[1]){stop("y and x should be in same sample size! and y should not have more than one dimension")}
    
    
    ## Getting the significant genes
    ER_Res_abs    <- abs(ER_Res$A)
    ii            <- order(abs(ER_Res_abs[,k]),decreasing = T)
    sigGenes      <- which(ER_Res_abs[ii,k]>thre)
    sigGenesload  <- ER_Res_abs[names(sigGenes),k,drop=F]
    colnames(sigGenesload) <-"A"
    
    
    ## Getting the correlation
    report <- NULL
    for(i in names(sigGenes)){
      
      report <- rbind(report,
                      data.frame(correlation=unname(cor(x[,i],y,method = "spearman")),
                                 sign=unname(sign(cor(x[,i],y,method = "spearman")))))
      
    }  
    
    report$color <- ifelse(report$sign<0,negcol,posCol)
    corReport <- data.frame(load=sigGenesload,report)
    
    
    if(orderbycol){
      
      corReport <- corReport[order(corReport$color),]
      
    }
    
    corReport$gene <- row.names(corReport)
    
    
    
    

# Addin filteration to the genes ------------------------------------------

    
    
    if(!is.null(filterTop)){
      
      if(!baseonCorrlation || !basedonLoad){
        
        stop("When you want to filter the gene both of the correlation and load cannot be False")
      }
      
      if(baseonCorrlation & !basedonLoad){
        
        corReport <- corReport[order(corReport$correlation,decreasing = T),][1:filterTop,]     ## Get top correlation
        
        
      }
      
      
      if(!baseonCorrlation & basedonLoad){
        
        corReport <- corReport[order(abs(corReport$A),decreasing = T),][1:filterTop,]               ## Get top loading
        
        
      }
      
      
      if(baseonCorrlation & basedonLoad){                                        ## Get the top correlation
        ## and top loading
        corReport1 <- corReport[order(corReport$A,decreasing = T),][1:filterTop,]   
        corReport2 <- corReport[order(corReport$correlation,decreasing = T),][1:filterTop,]
        corReport <- rbind(corReport1,corReport2)
        corReport <- corReport[!duplicated(corReport$gene),]
        
      }
      
      
      
      
      
      
    }
    
    
    
    
    
    
    
    
    return(corReport)
  
}

