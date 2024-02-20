checkDataParams("/ix/djishnu/Hanxi/SLIDE/test/test.yaml")
#zeroFiltering(yaml_path = "/ix/djishnu/Hanxi/SLIDE/test/test.yaml", 24, 804)

#yaml_path = "/ix/djishnu/Hanxi/test/SLIDE_main_s3Vs2Day28.yaml"
yaml_path = "/ix/djishnu/Hanxi/SLIDE/test/test.yaml"
input_params <- yaml::yaml.load_file(yaml_path)
optimizeSLIDE(input_params, sink_file = FALSE)
