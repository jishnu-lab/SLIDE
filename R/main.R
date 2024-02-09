#' Main Pipeline of Running SLIDE (get LF and run SLIDE without CV)
#'
#' @param yaml_path the string to a yaml path
#' @param sink_file boolean flag of saving to a sink file or not
#' @export

##################################### set up the parameters #####################################

main <- function(yaml_path, sink_file){
  
  input_params <- yaml::yaml.load_file(yaml_path)
  ##################################### check and print key parameters #####################################
  # check if output path exists
  if (dir.exists(input_params$out_path)){
    cat("Populating all outputs to ", input_params$out_path, '.\n')
  } else{
    cat ("Output folder not found, creating at ", input_params$out_path, ".\n")
    dir.create(file.path(input_params$out_path), showWarnings = F, recursive = T)
  }
  
  
  if (is.null(input_params$delta) == TRUE){delta = c(0.01, 0.1)}else{delta = input_params$delta}
  if (is.null(input_params$lambda) == TRUE){lambda = c(0.5, 1.0)}else{lambda = input_params$lambda}
  
  # cat("Input Delta:\n")
  # for (item in delta) {
  #   cat(item, "\n")
  # }
  # cat("Input Lambda:\n")
  # for (item in lambda) {
  #   cat(item, "\n")
  # }
  
  # take care of all input params 
  if (is.null(input_params$alpha)){alpha_level = 0.05} else {alpha_level = input_params$alpha}
  if (is.null(input_params$thresh_fdr)){thresh_fdr = 0.2} else {thresh_fdr = input_params$thresh_fdr}
  if (is.null(input_params$rep_cv)){rep_cv = 50} else {rep_cv = input_params$rep_cv}
  if (is.null(input_params$spec)){spec = 0.1} else(spec = input_params$spec)
  if (is.null(input_params$do_interacts)){do_interacts = TRUE} else(do_interacts = do_interacts)
  if (is.null(input_params$sigma)){
    sigma = NULL
    cat("Setting sigma as Null.\n")
  } else {sigma = input_params$rep_cv}
  
  if (is.null(input_params$SLIDE_iter)){SLIDE_iter = 500} else{SLIDE_iter = input_params$SLIDE_iter}
  if (is.null(input_params$eval_type)){
    if (length(unique(y))<=2){eval_type = "auc"}
    else{eval_type = "corr"}
  }else{eval_type = input_params$eval_type}
  if (is.null(input_params$SLIDE_top_feats)){SLIDE_top_feats = 10} else {SLIDE_top_feats = input_params$SLIDE_top_feats}
  

  ##################################### Code #####################################
  x <- as.matrix(utils::read.csv(input_params$x_path, row.names = 1))
  y <- as.matrix(utils::read.csv(input_params$y_path, row.names = 1))
  x_std <- scale(x, T, T)
  
  #initiate the summary table
  summary_table <- as.data.frame(matrix(NA, nrow = length(delta) * length(lambda), ncol = 6))
  colnames(summary_table) <- c('delta', 'lambda', '# of LFs', '# of Sig LFs', '# of Interactors', 'sampleCV Performance')
  
  i = 1
  for (d in delta){
    for (l in  lambda){
      ##################################### Get LF.#####################################
      loop_outpath = paste0(input_params$out_path, '/', d, '_', l, '_', 'out/' )
      dir.create(file.path(loop_outpath), showWarnings = F, recursive = T)
      
      if(sink_file == TRUE){
        sink(paste0(loop_outpath, "/standard_out.txt"))
      }
      cat("Getting latent factors for delta, ", d, ", and lambda, ", l, ". \n")
      cat("Setting alpha_level at ", alpha_level, ".\n")
      cat("Setting thresh_fdr at ", thresh_fdr, ".\n")
      cat("Setting rep_cv at ", rep_cv, ".\n")
      cat("Setting spec at ", spec, ".\n")
      cat("Setting eval_type as ", eval_type, ".\n")
      cat("Setting SLIDE_iter at ", SLIDE_iter, ".\n")
      cat("Setting SLIDE_top_feats as ", SLIDE_top_feats, ".\n")
      
      
      if (input_params$y_factor) {
        y_temp <- toCont(y, input_params$y_order)
        
        saveRDS(y_temp, file = paste0(input_params$out_path, "plainER_y_mapping.rds"))
        orig_y <- as.matrix(y_temp$cat_y)
        y <- as.matrix(y_temp$cont_y)
        row.names(y) <- row.names(y_temp$cat_y)
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
      
      saveRDS(all_latent_factors, paste0(loop_outpath, 'AllLatentFactors.rds'))
      
      # get Z matrix
      z_matrix <- calcZMatrix(x_std, all_latent_factors, x_path = NULL, lf_path = NULL, loop_outpath)
      
      # run SLIDE
      SLIDE_res <- runSLIDE(y, y_path = NULL, z_path = NULL, z_matrix, all_latent_factors, lf_path = NULL, niter = SLIDE_iter, spec = spec, do_interacts=do_interacts)
      saveRDS(SLIDE_res, paste0(loop_outpath, 'SLIDE_LFs.rds'))
      
      # get top features txt files and latent factor plots
      SLIDE_res <- getTopFeatures(x, y, all_latent_factors, loop_outpath, SLIDE_res, num_top_feats = SLIDE_top_feats, condition = eval_type)
      plotSigGenes(SLIDE_res, plot_interaction = TRUE, out_path = loop_outpath)
      
      #the SLIDE_res has to be the output from getTopFeatures
      #calculate the control performance plot
      calcControlPerformance(z_matrix = z_matrix, y, SLIDE_res, niter = SLIDE_iter, condition = eval_type, loop_outpath)
      
      # calculate the sampleCV performance 
      performance = sampleCV(y, z_matrix, SLIDE_res, fraction = 2/3, condition = eval_type, sampleCV_iter = 20, logistic = FALSE, out_path = loop_outpath)
      
      # fill in the summary table
      interactors = c(SLIDE_res$interaction$p1, SLIDE_res$interaction$p2)[which(!(c(SLIDE_res$interaction$p1, SLIDE_res$interaction$p2) %in% SLIDE_res$marginal_vals))]
      if (sum(interactors %in% SLIDE_res$marginal_vals) != 0) {stop("getting interactor code is wrong.")}
    
      loop_summary = c(d, l, all_latent_factors$K, length(SLIDE_res$marginal_vals), length(interactors), performance)
      summary_table[i, ] = loop_summary
      i = i+1
    }
  }
  write.csv(summary_table, paste0(input_params$out_path, "/summary_table.csv"))
  return(summary_table)
}

