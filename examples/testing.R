yaml_path = "/ix/djishnu/Hanxi/SLIDE/examples/test.yaml"
checkDataParams(yaml_path)
input_params <- yaml::yaml.load_file(yaml_path)
optimizeSLIDE(input_params, sink_file = FALSE)
