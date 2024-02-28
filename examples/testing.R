yaml_path = "/ix/djishnu/Hanxi/SLIDE/examples/example_continuous.yaml"
input_params <- yaml::yaml.load_file(yaml_path)

checkDataParams(input_params)
optimizeSLIDE(input_params, sink_file = FALSE)
