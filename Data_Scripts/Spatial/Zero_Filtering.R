rm(list = ls())
cat("\014")
library(ggplot2)

x <- read.csv('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/Spatial_ER_061722/Data/Concat.csv',
              row.names = 1)


hist(colSums(x == 0, na.rm=TRUE))

i <- colSums(x == 0, na.rm=TRUE) < 900
myData <- x[i]
hist(colSums(myData == 0, na.rm=TRUE))

write.csv(myData, '/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Poholek/Spatial/Spatial_ER_061722/Data/ER_X.csv')
