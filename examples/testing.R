yaml_path = "/ix/djishnu/Hanxi/SLIDE/examples/example_binary.yaml"
input_params <- yaml::yaml.load_file(yaml_path)

SLIDE::checkDataParams(input_params)

## Approximate the SLIDE model's performance for different delta and lambda. 
SLIDE::optimizeSLIDE(input_params, sink_file = FALSE)

## Finding the  SLIDE model's performance  using cross-validation
SLIDE::SLIDEcv(yaml_path, nrep = 20, k = 5)

## Depending on the limitations you can run optimizeSLIDE or slidecv


## You can also run the ```optimizeSLIDE``` and for a single combination of delta and lambda you can run ```SLIDEcv```





SLIDE::plotCorrelationNetworks(input_params)

