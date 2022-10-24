library(tidyverse)
## Read Immport data

immport_cel_seq_meta_data <- read_tsv('/bgfs/djishnu/Priyamvada/Immport/celseq_meta.tsv.725591')
## Modify sample name to include disease type 
#immport_cel_seq_meta_data$sample_diseasetype <- paste(immport_cel_seq_meta_data$sample,'-', immport_cel_seq_meta_data$disease)
## Map sample by patient
sample_plate<- immport_cel_seq_meta_data %>% select(sample,plate) %>% table()
sample_plate_sum <- apply(sample_plate,1,sum)


sample_plate %>% group_by(sample)


sample_cell_map <- immport_cel_seq_meta_data %>% group_by(sample)%>% select(cell_name)

patient_sample <- list()
for (i in unique(sample_cell_map$sample)){
  
patient_sample <-c(patient_sample,list(sample_cell_map[sample_cell_map$sample==i,"cell_name"]))
  
  
  
}
#sample_plat_column<- apply(sample_plate,2,function(x){if(sum(x)>0) return(1)})

names(patient_sample) <- unique(sample_cell_map$sample)
 
## filter data for each patient And save

cellSeq <- read.delim('/bgfs/djishnu/Priyamvada/Immport/celseq_matrix_ru10_molecules.tsv.725585.gz')

for(i in  unique(sample_cell_map$sample)){
  
  write.table(cellSeq[,unlist(patient_sample[[i]])],file=paste0(i,".csv"),col.names = T,row.names = T,quote = F)
  
}


## save
print(sample_plate_sum)
summary_cel_seq_meta_data_filtered <- immport_cel_seq_meta_data %>% group_by(disease, sample) %>% summarise(n = n())
