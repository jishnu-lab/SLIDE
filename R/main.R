##################################### set up the parameters #####################################
input_params <- yaml::yaml.load_file("/ix/djishnu/Hanxi/SLIDE/test/test.yaml")
sink_file = FALSE

# check if output path exists
if (dir.exists(input_params$out_path)){
     cat("Populating all outputs to ", input_params$out_path, '.\n')
   } else{
    cat ("Output path not found, creating an output folder at ", input_params$out_path, ".\n")
    dir.create(file.path(input_params$out_path), showWarnings = F, recursive = T)
   }

if(sink_file == TRUE){
  sink(paste0(input_params$out_path, "/standarad_put.txt"))
}

cat("Input Delta:\n")
for (item in input_params$delta) {
  cat(item, "\n")
  }
cat("Input Lambda:\n")
for (item in input_params$lambda) {
  cat(item, "\n")
  }

if (is.null(input_params$alpha) == TRUE){alpha_level = 0.05} else {alpha_level = input_params$alpha}
if (is.null(input_params$thresh_fdr) == TRUE){thresh_fdr = 0.2} else {thresh_fdr = input_params$thresh_fdr}
if (is.null(input_params$rep_cv) == TRUE){rep_cv = 50} else {rep_cv = input_params$rep_cv}
if (is.null(input_params$sigma) == TRUE){
  sigma = NULL
  cat("Setting sigma as Null.\n")
  } else {sigma = input_params$rep_cv}
if (is.null(input_params$spec) == TRUE){spec = 0.1} else {spec = input_params$spec}

cat("Setting alpha_level at ", alpha_level, ".\n")
cat("Setting thresh_fdr at ", thresh_fdr, ".\n")
cat("Setting rep_cv at ", rep_cv, ".\n")
##################################### Code #####################################
x <- as.matrix(utils::read.csv(input_params$x_path, row.names = 1))
y <- as.matrix(utils::read.csv(input_params$y_path, row.names = 1))
x_std <- scale(x, T, T)
y_std <- scale(y, T, T)


for (d in input_params$delta){
  for (l in  input_params$lambda){
    ##################################### Get LF.#####################################
    cat("Getting latent factors for delta, ", d, ", and lambda, ", l, ". \n")
    
    dir.create(file.path(paste0(input_params$out_path, '/', d, '_', l, '_', 'out' )), showWarnings = F, recursive = T)
    if (input_params$y_factor) {
      y <- toCont(y, input_params$y_order)
      saveRDS(y, file = paste0(input_params$out_path, "plainER_y_mapping.rds"))
      orig_y <- y$cat_y
      y <- y$cont_y
    }
    
    ## final output
    # er_output <- plainER(x = x,
    #                      x_std = x_std,
    #                      y = y,
    #                      std_y = y_std,
    #                      sigma = NULL,
    #                      delta = d,
    #                      lambda = l,
    #                      rep_cv = input_params$rep_cv,
    #                      alpha_level = input_params$alpha_level,
    #                      thresh_fdr = input_params$thresh_fdr,
    #                      out_path = input_params$out_path)
    # 
    # saveRDS(er_output, paste0(input_params$out_path, 'final_er_output.rds'))
    # 
  }
}
