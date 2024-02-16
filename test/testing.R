library(doParallel)
cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# cores <- 6
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

checkDataParams("/ix/djishnu/Hanxi/SLIDE/test/test.yaml")
#zeroFiltering(yaml_path = "/ix/djishnu/Hanxi/SLIDE/test/test.yaml", 24, 804)

yaml_path = "/ix/djishnu/Hanxi/test/SLIDE_main_s3Vs2Day28.yaml"
input_params <- yaml::yaml.load_file(yaml_path)
Optimize_SLIDE(input_params, sink_file = FALSE)
