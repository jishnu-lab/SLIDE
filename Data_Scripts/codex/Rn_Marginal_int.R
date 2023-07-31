library(SLIDE)
library(pROC)
## a) Real, b) Random size-matched marginal + interactors,
################################################################################
Z <- read.table("/ix/djishnu/Javad/SLIDEbench_Codex/z_matrix.csv", sep = ",", header = T, row.names = 1)
y <- read.table("/ix/djishnu/Javad/SLIDEbench_Codex/y.csv", row.names = 1, sep = ",", header = T)
colnames(y) <- "y"
zz <- row.names(Z)
y <- y[zz,,drop=F]
source("/ix/djishnu/Javad/Hierarchical_ER_v2/R/pairwise_interactions.R")


# Real ##########################################################################

# put your marginal and interaction variables here!
sigK <- c(1,5)
sigIn <- c("Z5.Z2","Z5.Z3","Z5.Z4","Z5.Z6","Z5.Z7","Z5.Z8","Z5.Z9")



IntData <- pairwise_interactions(sigK,Z)
Dataint <- IntData$interaction[, sigIn]
Data_real <- data.frame(y = y, Z[, sigK], Dataint)
lmod  <- lm(y~.,data=Data_real)
yhat <- predict(lmod,Data_real[,-1],type = 'response')
aucreal <- auc(response=as.matrix(y), predictor=as.matrix(yhat))


# All random ####################################################################
Fullreport <- NULL

for (i in 1:1000) {
  sigKRandom <- sample(ncol(Z), size = length(sigK)) ## Random marginal
  IntDataRandom <- pairwise_interactions(sigKRandom, Z)
  sigInRandom <- sample(ncol(IntDataRandom$interaction), size = length(sigIn)) ## Ranodm interaction
  IntDataRandom <- IntDataRandom$interaction[, sigInRandom]
  Data_fullRandom <- data.frame(y = y, Z[, sigKRandom], IntDataRandom)
  SumInt <- summary(lm(y ~ ., data = Data_fullRandom))
  lmod  <- lm(y~.,data=Data_fullRandom)
  yhat <- predict(lmod,Data_fullRandom[,-1],type = 'response')
  aucrandom <- auc(response=as.matrix(y), predictor=as.matrix(yhat))
  Fullreport <- rbind(Fullreport, aucrandom)
}


Partialreport <- NULL
for (i in 1:1000) {
 
  IntData  <- pairwise_interactions(sigK, Z)
  sigInRandom <- sample(ncol(IntData$interaction), size = length(sigIn)) ## Ranodm interaction
  IntDataRandom <- IntData$interaction[, sigInRandom]
  Data_partialRandom <- data.frame(y = y, Z[, sigK], IntDataRandom)
  SumPInt <- summary(lm(y ~ ., data = Data_partialRandom))
  lmod  <- lm(y~.,data=Data_partialRandom)
  yhat <- predict(lmod,Data_partialRandom[,-1],type = 'response')
  aucPrandom <- auc(response=as.matrix(y), predictor=as.matrix(yhat))
  Partialreport <- rbind(Partialreport, aucPrandom)
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
  geom_vline(xintercept         =aucreal , 
             linetype = "longdash",
             colour = "red",size=2) +
  ylab("Density") +
  xlab("auc")
P2 <- P2 + annotate("text", x     = aucreal + 0.01, 
                    y = 55, 
                    label = " ", 
                    angle = "90") +
  xlim(0.25, 0.85)
P2 <- P2 + theme(panel.border = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 axis.text = element_text(size = 20),
                 axis.title.x =element_text(size = 20),
                 axis.title.y = element_text(size = 20))




P2

saveRDS(df, file                = "/ix/djishnu/Javad/SLIDEbench_Codex//ggplotobject_randompar_data2.rds")
saveRDS(P2, file                 = "/ix/djishnu/Javad/SLIDEbench_Codex/ggplotobject_randompar2.rds")
ggsave(P2, filename              = "/ix/djishnu/Javad/SLIDEbench_Codex/randomvsPartialPlot_limited2.pdf")
ggsave(P2, filename              = "/ix/djishnu/Javad/SLIDEbench_Codex/randomvsPartialPlot_limited2.png")
# 
