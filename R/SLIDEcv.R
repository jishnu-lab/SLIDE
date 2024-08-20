#' Cross validation for SLIDE latent factors
#'
#' @param yaml_path the string to a yaml path where optimizeSLIDE has been run; can either be the original yaml file (to run CV on all runs),
#'  or the yaml for a specific run
#' @param nrep number of replicates for cross validation (default 20)
#' @param k number of folds for cross validation (default 10)
#' @export

SLIDEcv <- function(yaml_path=NULL, nrep=20, k=10){

  if(is.null(yaml_path)){stop("\nyaml path can't be empty! \n")}

  slide_input <- yaml::yaml.load_file(as.character(yaml_path))

  # check if we are running SLIDE on all the runs from optimizeSLIDE or just one run
  run_dirs = list.files(slide_input$out_path, pattern = "[0-9]+(\\.)?[0-9]+?_[0-9]+(\\.)?([0-9]+)?_out|*_out$", include.dirs = TRUE,
                        full.names = TRUE, recursive = FALSE)

  # make sure they are directories
  run_dirs = base::intersect(run_dirs, list.dirs(slide_input$out_path, recursive = FALSE, full.names = TRUE))

  if (length(run_dirs) == 0) {
    # if we don't find any output folders from optimizeSLIDE, just use the outpath in the yaml (this is case when
    # running SLIDEcv using a standard yaml file/from the yaml file in a single run from optimizeSLIDE)
    run_dirs = slide_input$out_path
  }
  for (r in run_dirs) {
    slide_input$out_path = paste0(r, "/")

    cat("\nRunning SLIDE cross validation for ", r, "\n")
    y <- as.matrix(utils::read.csv(slide_input$y_path, row.names = 1)) ## not standardized

    x <- as.matrix(utils::read.csv(slide_input$x_path, row.names = 1)) ## not standardized

    cat("\nRunning SLIDE cross validation for ", r, "\n")

    # check for SLIDE results
    all_latent_factors_path = list.files(r, pattern = "AllLatentFactors.rds", full.names = T, recursive = F)
    z_matrix_path = list.files(r, pattern = "z_matrix.csv", full.names = T, recursive = F)

    if (length(all_latent_factors_path) == 0) {
      stop("\nNo file for AllLatentFactors.rds found. Run SLIDE first or check path \n")
    } else if (length(z_matrix_path) == 0) {
      stop("\nNo file for z_matrix.csv found. Run SLIDE first or check path \n")
    }

    all_latent_factors = readRDS(all_latent_factors_path)
    z = as.matrix(read.csv(z_matrix_path, row.names = 1))

#
#     all_latent_factors <- getLatentFactors(x=x,
#                      y=y,
#                      std_y = T,
#                      x_std = scale(x,T,T),
#                      delta =slide_input$delta,
#                      lambda = slide_input$lambda,
#                      sigma = NULL)
#
#
#     z <- calcZMatrix(x_std = scale(x,T,T),all_latent_factors = all_latent_factors,out_path = slide_input$out_path)

    z <- scale(z,T,T)

    f_size <- calcDefaultFsize(y=y,all_latent_factors = all_latent_factors )

    slide_input$std_cv <- T
    slide_input$std_y <- T
    slide_input$permute <- T
    slide_input$parallel <- T
    slide_input$fdr <- 0.1
    slide_input$do_interacts = ifelse(is.null(slide_input$do_interacts), TRUE, slide_input$do_interacts)


    # make a new directory to store the replicate folders for SLIDECV

    folder_for_slide_reps = paste0(r, "/SLIDE_CV_replicate_data/")
    if (!dir.exists(folder_for_slide_reps)) {
      dir.create(folder_for_slide_reps)
    }

    for(j in 1:nrep){

    benchCV2(k =k,
             x = x,
             y = y,
             z=z,
             std_y = slide_input$std_y,
             std_cv = slide_input$std_cv,
             delta = slide_input$delta,
             permute = slide_input$permute,
             eval_type = slide_input$eval_type,
             y_levels = slide_input$y_levels,
             lambda = slide_input$lambda,
             out_path = folder_for_slide_reps,
             rep_cv = slide_input$rep_cv,
             alpha_level = slide_input$alpha_level,
             thresh_fdr = slide_input$thresh_fdr,
             rep = j,
             benchmark = F,
             niter=slide_input$SLIDE_iter,
             spec=slide_input$spec,
             fdr=slide_input$fdr,
             f_size=f_size,
             parallel = slide_input$parallel,
             ncore=20,
            do_interacts = slide_input$do_interacts)

    }


    pathLists <- list.files(folder_for_slide_reps,recursive = T,pattern = "results", full.names = T)
    perfList <- lapply(pathLists, readRDS)
    perRes <- do.call(rbind,lapply(perfList,function(x){x$final_corr}))


    if (slide_input$eval_type == "corr") {
      lambda_boxplot = ggpubr::ggboxplot(data = perRes, x = "method", y = "corr", palette = "aaas",
                                         fill = "method" ) +
        ggpubr::stat_compare_means(label = "p.signif")

    } else {
      lambda_boxplot = ggpubr::ggboxplot(data = perRes, x = "method", y = "auc", palette = "aaas",
                                         fill = "method" ) +
        ggpubr::stat_compare_means(label = "p.signif")
    }
    ggplot2::ggsave(plot = lambda_boxplot, filename = paste0(slide_input$out_path, "SLIDECV_boxplot.pdf"), height = 6, width = 6)
    saveRDS(perRes,file=paste0(slide_input$out_path,"SLIDECV_boxplot_data.rds"))
  }
}
