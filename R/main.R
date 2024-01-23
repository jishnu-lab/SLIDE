##################################### set up the parameters #####################################
input_params <- yaml::yaml.load_file("/ix/djishnu/Hanxi/SLIDE/test/test.yaml")
sink_file = FALSE

##################################### check and print key parameters #####################################
# check if output path exists
if (dir.exists(input_params$out_path)){
     cat("Populating all outputs to ", input_params$out_path, '.\n')
   } else{
    cat ("Output path not found, creating an output folder at ", input_params$out_path, ".\n")
    dir.create(file.path(input_params$out_path), showWarnings = F, recursive = T)
   }

if(sink_file == TRUE){
  sink(paste0(input_params$out_path, "/standarad_out.txt"))
}


if (is.null(input_params$delta) == TRUE){delta = c(0.01, 0.1)}else{delta = input_params$delta}
if (is.null(input_params$lambda) == TRUE){lambda = c(0.5, 1.0)}else{lambda = input_params$lambda}

cat("Input Delta:\n")
for (item in delta) {
  cat(item, "\n")
  }
cat("Input Lambda:\n")
for (item in lambda) {
  cat(item, "\n")
  }

# take care of all input params 
if (is.null(input_params$alpha)){alpha_level = 0.05} else {alpha_level = input_params$alpha}
if (is.null(input_params$thresh_fdr)){thresh_fdr = 0.2} else {thresh_fdr = input_params$thresh_fdr}
if (is.null(input_params$rep_cv)){rep_cv = 50} else {rep_cv = input_params$rep_cv}
if(is.null(input_params$spec) ){spec = 0.1} else(spec = input_params$spec)
if (is.null(input_params$sigma)){
  sigma = NULL
  cat("Setting sigma as Null.\n")
  } else {sigma = input_params$rep_cv}
if (is.null(input_params$spec)){spec = 0.1} else {spec = input_params$spec}
if (is.null(input_params$SLIDE_iter)){SLIDE_iter = 500} else{SLIDE_iter = input_params$SLIDE_iter}

cat("Setting alpha_level at ", alpha_level, ".\n")
cat("Setting thresh_fdr at ", thresh_fdr, ".\n")
cat("Setting rep_cv at ", rep_cv, ".\n")
cat("Setting spec at ", spec, ".\n")
cat("Setting SLIDE_iter at ", SLIDE_iter, ".\n")



##################################### Code #####################################
x <- as.matrix(utils::read.csv(input_params$x_path, row.names = 1))
y <- as.matrix(utils::read.csv(input_params$y_path, row.names = 1))
x_std <- scale(x, T, T)


for (d in delta){
  for (l in  lambda){
    ##################################### Get LF.#####################################
    cat("Getting latent factors for delta, ", d, ", and lambda, ", l, ". \n")
    
    loop_outpath = paste0(input_params$out_path, '/', d, '_', l, '_', 'out/' )
    dir.create(file.path(loop_outpath), showWarnings = F, recursive = T)
    
    if (input_params$y_factor) {
      y <- toCont(y, input_params$y_order)
      
      saveRDS(y, file = paste0(input_params$out_path, "plainER_y_mapping.rds"))
      orig_y <- y$cat_y
      y <- y$cont_y
    }
    
    #final output
    all_latent_factors <- getLatentFactors(x = x,
                         x_std = x_std,
                         y = y,
                         sigma = NULL,
                         delta = d,
                         lambda = l,
                         rep_cv = rep_cv,
                         alpha_level = alpha_level,
                         thresh_fdr = thresh_fdr,
                         out_path = loop_outpath)

    saveRDS(er_output, paste0(loop_outpath, 'AllLatentFactors.rds'))

    z_matrix <- calcZMatrix(x_std, all_latent_factors, x_path = NULL, lf_path = NULL, loop_outpath)
    SLIDE_res <- runSLIDE(y, y_path = NULL, z_path = NULL, z_matrix, all_latent_factors, lf_path = NULL, niter = SLIDE_iter)
    saveRDS(SLIDE_res, paste0(loop_outpath, 'SLIDE_LFs.rds'))

      }
}
