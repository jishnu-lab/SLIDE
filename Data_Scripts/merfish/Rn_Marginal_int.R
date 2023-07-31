library(SLIDE)
## a) Real, b) Random size-matched marginal + interactors,
################################################################################



Z <- read.table("/ix/djishnu/Javad/SLIDEbench_merfish/z_matrix.csv", sep = ",", header = T, row.names = 1)
y <- read.table("/ix/djishnu/Javad/SLIDEbench_merfish/y.csv", row.names = 1, sep = ",", header = T)
colnames(y) <- "y"
source("/ix/djishnu/Javad/Hierarchical_ER_v2/R/pairwise_interactions.R")


# Real ##########################################################################

# put your marginal and interaction variables here!
sigK <- c(10,12)
sigIn <- c("Z10.Z12","Z10.Z18","Z10.Z21","Z10.Z27","Z12.Z1","Z12.Z8","Z12.Z10","Z12.Z18","Z12.Z19","Z12.Z21","Z12.Z25")


IntData <- pairwise_interactions(sigK,Z)
Dataint <- IntData$interaction[, sigIn]
Data_real <- data.frame(y = y, Z[, sigK], Dataint)

SumReal <- summary(lm(y ~ ., data = Data_real))
SumReal$r.squared


# All random ####################################################################
Fullreport <- NULL

for (i in 1:100) {
  sigKRandom <- sample(ncol(Z), size = length(sigK)) ## Random marginal
  IntDataRandom <- pairwise_interactions(sigKRandom, Z)
  sigInRandom <- sample(ncol(IntDataRandom$interaction), size = length(sigIn)) ## Ranodm interaction
  IntDataRandom <- IntDataRandom$interaction[, sigInRandom]
  Data_fullRandom <- data.frame(y = y, Z[, sigKRandom], IntDataRandom)
  SumInt <- summary(lm(y ~ ., data = Data_fullRandom))
  SumInt$r.squared
  Fullreport <- rbind(Fullreport, sqrt(SumInt$r.squared))
}


Partialreport <- NULL
for (i in 1:100) {
  #IntData <- interUnion(sigK, Z)
  IntData  <- pairwise_interactions(sigK, Z)
  sigInRandom <- sample(ncol(IntData$interaction), size = length(sigIn)) ## Ranodm interaction
  IntDataRandom <- IntData$interaction[, sigInRandom]
  Data_partialRandom <- data.frame(y = y, Z[, sigK], IntDataRandom)
  SumPInt <- summary(lm(y ~ ., data = Data_partialRandom))
  SumPInt$r.squared
  Partialreport <- rbind(Partialreport, sqrt(SumPInt$r.squared))
}
################################################################################
## Report
library(ggplot2)
library(reshape2)
rawdf <- data.frame(FullRandom = Fullreport, PartialRandom = Partialreport)
df <- melt(rawdf)
colnames(df) <- c("group", "value")

cols <- c("#0000FF", "#00FF00")

# Basic density plot in ggplot2

P2 <- ggplot(df, aes(x           = value, fill = group)) +
  geom_density(alpha            = 0.7) +
  scale_fill_manual(values      = cols) +
  theme_light() +
  geom_vline(xintercept         = c(sqrt(SumReal$r.squared)), 
             linetype = "longdash",
             colour = "red",size=2) +
  ylab("Density") +
  xlab("Pearson Correlation")
P2 <- P2 + annotate("text", x     = sqrt(SumReal$r.squared) + 0.01, 
                    y = 55, 
                    label = " ", 
                    angle = "90") +
  xlim(0.25, 1)
P2 <- P2 + theme(panel.border = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 axis.text = element_text(size = 20),
                 axis.title.x =element_text(size = 20),
                 axis.title.y = element_text(size = 20))




P2

saveRDS(df, file                = "/ix/djishnu/Javad/SLIDEbench_merfish/ggplotobject_randompar_data2.rds")
saveRDS(P2, file                 = "/ix/djishnu/Javad/SLIDEbench_merfish/ggplotobject_randompar2.rds")
ggsave(P2, filename              = "/ix/djishnu/Javad/SLIDEbench_merfish/randomvsPartialPlot_limited2.pdf")
ggsave(P2, filename              = "/ix/djishnu/Javad/SLIDEbench_merfish/randomvsPartialPlot_limited2.png")
# 
