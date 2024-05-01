yaml_path = "/ix/djishnu/Hanxi/SLIDE/examples/example_binary.yaml"
input_params <- yaml::yaml.load_file(yaml_path)

SLIDE::checkDataParams(input_params)
SLIDE::optimizeSLIDE(input_params, sink_file = FALSE)
SLIDE::plotCorrelationNetworks(input_params)
SLIDE::SLIDEcv(yaml_path, nrep = 20, k = 5)
