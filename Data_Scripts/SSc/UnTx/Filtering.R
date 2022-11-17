data <- read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Lafyatis/All_Cell_Type/HER_081222/Data/Var50_thresh5.csv",
                 row.names = 1)

genes <- colnames(data)



for (i in 0:17){
  genes <- gsub(paste0("c.", i, "."), "", fixed = TRUE, genes)
}


idx <- grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-|MTRNR|MT4|MT3|MT2A|MT1E|MT1M|MT1A|MT1B|MT1F|MT1G|MT1H|^MTND", genes)

for (i in 1:length(idx)){
  print(genes[idx[i]])
}

data <- data[, -idx]

write.csv(data, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Lafyatis/All_Cell_Type/HER_081222/Data/Var50_mtrp.csv")
