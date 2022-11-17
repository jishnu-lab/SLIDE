library(EssReg)
library(doParallel)


cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# cores <- 6
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

yaml_path = '/ix/djishnu/Hanxi/Joglekar/CD4_Expansion/HER_091822/pipeline1_no_scale.yaml'
pipelineER1(yaml_path, 'all')
