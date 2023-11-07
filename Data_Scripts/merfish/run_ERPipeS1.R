library(EssReg)
library(doParallel)


cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# cores <- 6
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

yaml_path = '/ix/djishnu/Hanxi/SLIDE_merFISH/ordinal_analysis/pipeline1.yaml'
pipelineER1(yaml_path, 'all')
