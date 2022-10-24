### Running v2 PLS regression on data

library(pls)

fam_data <- read.table("/bgfs/djishnu/Javad/racer/racer-genetics/racer-imputed/racer-mega-imputed-f-matched-v3-QC.fam")
PCs      <- read.table("/bgfs/djishnu/Javad/racer/racer-genetics/racer-imputed/racer-mega-imputed-f-matched-v2-PC_5.PC") 
PCs      <- as.matrix(PCs[,-c(1,2)]) 
y <- fam_data$V6-1
plsFit <- plsr(fam_data$V6~gene,validation="CV")
pls.RMSEP = RMSEP(plsFit, estimate="CV")
min_comp = which.min(pls.RMSEP$val)
lmod1 <- lm(y ~ plsFit$scores[,1:700]+PCs)
lmod2 <- lm(y ~ plsFit$scores[,1:700])

sum               <- anova(lmod2,lmod1,test = "LRT")
#write.table(cbind(i,sum$`Pr(>Chi)`[2]),quote = F,append = T,row.names=F,col.names=F,file=paste0("/bgfs/djishnu/Javad/racer/racer-genetics/racer-imputed/output/geneScore_Sig",assoc_thre,"pcrThre",pcThre,".csv"))



################ 
PtPtP <- (crossprod(PCs)%*% t(PCs))
M     <- diag(x=1,nrow = length(y))-(PCs %*% (PtPtP))
K     <- tcrossprod(PCs)
Kg    <- tcrossprod(gene)
X     <- Kg[lower.tri(Kg,diag = F)]
yM    <- M%*%y

summary(lm(Y~X))


##### Functions List



## CorrectionMat 

getCorrectionMap <- function(X){

PtPtP <- (crossprod(X)%*% t(X))
M     <- diag(x=1,nrow = length(y))-(X %*% (PtPtP))
return(M)

}


## Get features

getvectorFeature<- function(xM){
K_x      <- tcrossprod(xM)
vecX     <- K_x[lower.tri(K_x,diag=F)]
return(vecX)

}

## Correct Feature
correctFeature <- function(M,X){
  
  XM <- M%*%X
  return(XM)
  
  
}

doRegression<- function(vecX,VecY){
  
  Sum <- summary(lm(VecY~vecX))  
  Pval <- Sum$coefficients[,'Pr(>|t|)']
  
  Sum$coefficients
  
}

################################################################################
PCs <- read.table()

M  <- getCorrectionMap(PCs)
yM <- correctFeature(M,y)
xM <- correctFeature(M,gene)
vecX <- getvectorFeature(xM)
VecY <- getvectorFeature(yM)
res  <- doRegression(vecX,VecY)


#summary(lm(fitG$g~as.factor(y)))
