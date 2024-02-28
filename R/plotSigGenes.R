#' Plot significant genes for significant (marginal) latent factors
#'
#' @importFrom magrittr '%>%'
#' @param slide_results list - SLIDE results output from runSLIDE function
#' @param out_path string - Output path to save results
#' @param plot_interactions logical - whether to plot interaction variables
#' @export



plotSigGenes = function(slide_results, plot_interactions = F, out_path = NULL) {

  slide_vars = slide_results$feature_res

  # this is for the upper y limit in the plot
  max_num_genes_in_any_lf = 0

  sg_df = data.frame()

  # make results into plottable dataframe
  for (lf in 1:length(slide_vars)) {
    lf_df = list2DF(slide_vars[[lf]])
    lf_df[,"lf_num"] = as.numeric(stringr::str_replace(names(slide_vars)[lf],
                                                   pattern = "Z", replacement = ""))
    sg_df = rbind(sg_df, lf_df)

    max_num_genes_in_any_lf = ifelse(nrow(lf_df) > max_num_genes_in_any_lf,
                                     nrow(lf_df),
                                     max_num_genes_in_any_lf)
  }

  sg_plot_df = data.frame()
  for (lf in unique(sg_df$lf_num)) {
    # get genes for this latent factor
    lf_temp = sg_df[sg_df$lf_num == lf, ]

    # add the heights to each gene
    lf_temp$plot_height = seq(nrow(lf_temp), 1)

    sg_plot_df = rbind(sg_plot_df, lf_temp)
  }

  text_color_vals = sort(unique(stringr::str_to_lower(sg_plot_df$color)))

  plot_list = list()

if(!is.null(slide_results$SLIDE_res$marginal_vars)){

  slide_results$SLIDE_res$marginal_vars <- gsub("z", "", slide_results$SLIDE_res$marginal_vars) ## With interact equal to true we miss this
}
  # plot marginals
  marg_plot = sg_plot_df %>% dplyr::filter(sg_plot_df$lf_num %in% slide_results$SLIDE_res$marginal_vars) %>%
    ggplot2::ggplot(., ggplot2::aes(x = factor(lf_num), y = plot_height, label = names)) +
    ggplot2::geom_text(ggplot2::aes(color = factor(color))) +
    ggplot2::scale_color_manual(values = text_color_vals, guide = "none") + ggplot2::theme_void() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(), axis.title.x = ggplot2::element_text(),
                   axis.title.y = ggplot2::element_text(angle = 90)) +ggplot2::xlab("Significant Latent Factor") +
    ggplot2::ylab("Genes Associated with Significant Latent Factors") +
    ggplot2::ylim(0, max_num_genes_in_any_lf) +
    ggplot2::ggtitle("SLIDE Marginal Variables")

  plot_list[[1]] = marg_plot

  if (plot_interactions) {
    # plot all
    # add font face labels to data
    sg_plot_df$is_marginal = ifelse(sg_plot_df$lf_num %in% slide_results$SLIDE_res$marginal_vars,
                                    list(c("bold", "italic")), c("plain"))
    plt = sg_plot_df %>% ggplot2::ggplot(., ggplot2::aes(x = factor(lf_num), y = plot_height, label = names)) +
      ggplot2::geom_text(ggplot2::aes(color = factor(color),
                             fontface = ifelse(lf_num %in% slide_results$SLIDE_res$marginal_vars,
                                                                    "bold.italic", "plain"))) +
      ggplot2::scale_color_manual(values = text_color_vals, guide = "none") + ggplot2::theme_void() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(), axis.title.x = ggplot2::element_text(),
                     axis.title.y = ggplot2::element_text(angle = 90)) +
      ggplot2::xlab("Significant Latent Factor") +
      ggplot2::ylab("Genes Associated with Significant Latent Factors") +
      ggplot2::ylim(0, max_num_genes_in_any_lf) +
      ggplot2::ggtitle("Significant Latent Factors - Marginals (bold/italic) and Interactions")
    
    make_interaction_adj = function(slide_results) {
      # make edgelist
      edges = data.frame()
    
      # add marginals first
      for (e in slide_results$marginal_vars) {
        mvar = paste0("Z", e)
        elist = list(A = mvar, B = mvar, C = 0)
        edges = rbind.data.frame(edges, elist)
      }
    
      for (e in slide_results$interaction_vars) {
        elist = stringr::str_split(e, pattern = "\\.")[[1]]
        elist = list(A = elist[1], B = elist[2], C = 1)
        edges = rbind.data.frame(edges, elist)
      }
      return(edges)
    }
    
    edges = make_interaction_adj(slide_results$SLIDE_res)
    
    unique_nodes = unique(unlist(edges[, 1:2]))
    
    node_colors = ifelse(unique_nodes %in% paste0("Z", slide_results$SLIDE_res$marginal_vars),
                         "salmon", "lightgray")
    
    group_list = list(
      marginal = which(unique_nodes %in% paste0("Z", slide_results$SLIDE_res$marginal_vars)),
      interaction = which(unique_nodes %in% setdiff(unique_nodes, paste0("Z", slide_results$SLIDE_res$marginal_vars)))
    )
    
    lf_graph = qgraph::qgraph(edges, layout = "spring", directed = F, edge.color = "black", esize = 4,
               groups = group_list, color = c("salmon", "lightgray"), DoNotPlot = TRUE)

    plot_list[[2]] = plt
    plot_list[[3]] = lf_graph
  }


  if ( !is.null(out_path) ) {

    saveRDS(sg_plot_df, paste0(out_path, '/plotSigGenes_data.RDS'))

    ggplot2::ggsave(plot = marg_plot, filename = paste0(out_path, '/plotSigGenes_marginals.png'),

                    device = "png",
                    width = 1.5 * length(slide_results$SLIDE_res$marginal_vars), height = 7,
                    limitsize = FALSE)

    if (plot_interactions) {

      ggplot2::ggsave(plot = plt, filename = paste0(out_path, '/plotSigGenes.png'),
                      device = "png",
                      width = 1.5 * length(unique(sg_plot_df$lf_num)), height = 7,
                      limitsize = FALSE)
      
      # need to change wd to save correctly
      orig_wd = getwd()
      new_wd = out_path

      setwd(new_wd)
      qgraph::qgraph(edges, layout = "spring", directed = F, edge.color = "black", esize = 4,
           groups = group_list, color = c("salmon", "lightgray"), filename = "plotInteractions",
                    filetype = "png", height = 5, width = 7)
      
      setwd(orig_wd)
    }
  }
  return(plot_list)
}
