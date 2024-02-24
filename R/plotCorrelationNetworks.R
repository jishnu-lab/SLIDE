#' Plot correlation network for latent factors. Saves in a new directory "correlation_networks" in the output directory
#'
#' @param input_params loaded yaml file that has already been run through optimizeSLIDE
#' @return none
#' @export
#'
plotCorrelationNetworks = function(input_params) {

  original_wd = getwd()
  # iterate through each run (different values of delta/lambda)
  run_dirs = list.files(input_params$out_path, pattern = "[0-9]+(\\.)?[0-9]+?_[0-9]+(\\.)?([0-9]+)?_out", include.dirs = TRUE,
                        full.names = TRUE, recursive = FALSE)

  # make sure they are directories
  run_dirs = base::intersect(run_dirs, list.dirs(input_params$out_path, recursive = FALSE, full.names = TRUE))

  x = as.matrix(read.csv(input_params$x_path, row.names = 1))
  y = as.matrix(read.csv(input_params$y_path, row.names = 1))

  for (r in run_dirs) {
    # going to store correlation plots in a new folder
    dir_name = paste0(r, "/correlation_networks")

    if ( !dir.exists(dir_name) ) {
      dir.create(dir_name, recursive = T)
    }

    gene_list_files = list.files(r, pattern = "gene_list_Z", full.names = TRUE)

    if (length(gene_list_files) == 0) {
      # no gene lists found
      cat("\n No feature lists found. Run optimizeSLIDE first \n")
    }

<<<<<<< HEAD
=======
    setwd(dir_name)
>>>>>>> main

    for (f in gene_list_files) {

      LF_num = unlist(stringr::str_match(f, pattern = "Z[0-9]+"))
      temp_list = read.table(f, header = TRUE) %>% tidyr::drop_na()

<<<<<<< HEAD
      setwd(dir_name)

=======
>>>>>>> main
      # subset genes from expression data
      x_gene = as.matrix(x[, temp_list$names])

      # color code nodes as red = associated with higher Y/Y=1; blue = associated with lower Y/Y=0
      if (length(unique(y)) == 2) {
        col_auc = round(apply(x_gene, 2, function(xs) glmnet:::auc(as.matrix(y), as.matrix(xs))), 2)

<<<<<<< HEAD
        temp_cols = ifelse(col_auc > 0.5, "salmon", ifelse(col_auc < 0.5, "skyblue", "lightgray"))
      } else {
        col_cor = round(apply(x_gene, 2, function(xs) cor(as.matrix(y), as.matrix(xs), method = "spearman")), 2)

        temp_cols = ifelse(col_cor > 0, "salmon", ifelse(col_cor < 0, "skyblue", "lightgray"))
=======
        temp_cols = ifelse(col_auc > 0.55, "salmon", ifelse(col_auc < 0.45, "skyblue", "lightgray"))
      } else {
        col_cor = round(apply(x_gene, 2, function(xs) glmnet:::auc(as.matrix(y), as.matrix(xs))), 2)

        temp_cols = ifelse(col_cor > 0.1, "salmon", ifelse(col_cor < -0.1, "skyblue", "lightgray"))
>>>>>>> main
      }

      x_temp = cor(x_gene)

      pl = qgraph::qgraph(x_temp, filename=LF_num,
                          layout = "spring", minimum=0.25, repulsion=0.1,
                          labels = colnames(x_temp), color = temp_cols,
                          title = LF_num,
                          label.scale.equal=FALSE,label.prop=0.95,shape="ellipse",
                          posCol="#40006D", negCol="#59A14F",filetype='pdf',height=5,width=7)
<<<<<<< HEAD
      setwd(original_wd)

    }
  }
=======
    }
  }
  setwd(original_wd)
>>>>>>> main
}
