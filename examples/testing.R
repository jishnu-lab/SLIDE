yaml_path = "/ix/djishnu/Hanxi/SLIDE/examples/example_binary.yaml"
input_params <- yaml::yaml.load_file(yaml_path)

SLIDE::checkDataParams(input_params)


##################################################################
##                            STEP 1                            ##
##################################################################

##################################################################
##                           option 1                           ##
##################################################################
## Approximate the SLIDE model's performance for different delta and lambda. 
SLIDE::optimizeSLIDE(input_params, sink_file = FALSE)


##################################################################
##                           option 2                           ##
##################################################################

## Finding the  SLIDE model's performance  using cross-validation
SLIDE::SLIDEcv(yaml_path, nrep = 20, k = 5)

## Depending on the limitations you can run optimizeSLIDE or slidecv

##################################################################
##                           option 3                           ##
##################################################################

## You can also run the ```optimizeSLIDE``` and to find an approximate optimal model performance 
## Then you can evaluate the performance of the optimum combination of delta and lambda using ```SLIDEcv```


##################################################################
##                            STEP 2                            ##
##################################################################

## After finding the optimum delta and lambda one you can get a correlation network using the below command:
SLIDE::plotCorrelationNetworks(input_params)

