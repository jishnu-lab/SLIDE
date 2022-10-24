
############################################################################
############################################################################
###                                                                      ###
###                                MIC,R2                                ###
###                                                                      ###
############################################################################
############################################################################

library(ggplot2)
library(minerva)
library(reshape2)
library(RColorBrewer)
library(viridis)

source("/ix/djishnu/Javad/Hierarchical_ER/R/MakeFig.R", echo=TRUE)
ER_final <- readRDS("/ix/djishnu/Javad/Lafayatis/Final_SLIDE_081222/final_er_output.rds")
Y        <- read.delim("/ix/djishnu/Javad/Lafayatis/Final_SLIDE_081222/SkinScore_MRSS.csv",sep=",")
X        <- read.csv("/ix/djishnu/Javad/Lafayatis/Final_SLIDE_081222/Var50_mtrp.csv",row.names = 1)
breaksList = seq(0, 1, by = 0.1)



for(k in c(10,12,99,6,47,85,56)){
  
resA <- SigGenes_A(ER_final,k=k,X,Y[,2],thre=0.1,negcol="blue",posCol="red",orderbycol=T)

geneNames  <- row.names(resA)
report <- data.frame(R2=numeric(),MIC_R2=numeric())

for(i in 1:length(geneNames)){

res <- mine(X[,geneNames[i],drop=F],Y[,2])
res$R2 <- (cor(X[,geneNames[i],drop=F],Y[,2]))^2
report <- rbind(report,data.frame(R2=res$R2,MIC_R2=res$`MIC-R2`))

}

row.names(report) <- geneNames
colnames(report)  <- c("R2","MIC-R2")
mat_breaks <- seq(100, 0, length.out = 50)
breaksList <- seq(0,0.5,length.out = 50)

png(sprintf("/ix/djishnu/Javad/Lafayatis/Plots/LV%s.png",k))
p<- pheatmap::pheatmap(report[order(report$R2,decreasing = T),],
                   color= colorRampPalette(c("#FFE900", "black"))(100)[mat_breaks],
                   cluster_rows = F,
                   cluster_cols = F,
                   legend = T,
                   main=paste0("Z",k),
                   breaks = breaksList)


#print(p)
dev.off()
}


