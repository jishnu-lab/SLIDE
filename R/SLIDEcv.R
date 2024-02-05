require(doParallel)
require(dplyr)
library(yaml)


SLIDEcv <- function(yaml_path){


  sl_input <- yaml::yaml.load_file(as.character(yaml_path))




for(j in 1:sl_input$nreps){

  paperBench(yaml_path, replicate=j)

}

  sl_input <- yaml::yaml.load_file(yaml_path)

pathLists <- list.files(sl_input$out_path,recursive = T,pattern = "results")
perfList <- lapply(paste0(sl_input$out_path,pathLists), readRDS)
perRes <- do.call(rbind,lapply(perfList,function(x){x$final_corr}))


if (sl_input$eval_type == "corr") {
  lambda_boxplot <- ggplot2::ggplot(data = perRes,
                                    ggplot2::aes(x = method,
                                                 y = corr,
                                                 fill = method)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(fill = "Method") +
    ggplot2::scale_alpha(guide = 'none')
} else {
  lambda_boxplot <- ggplot2::ggplot(data = perRes,
                                    ggplot2::aes(x = method,
                                                 y = auc,
                                                 fill = method)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(fill = "Method") +
    ggplot2::scale_alpha(guide = 'none')
}
library(ggplot2)
ggsave(paste0(sl_input$out_path,"/delta",sl_input$delta,"lambda",sl_input$lambda,"_boxplot.pdf"),lambda_boxplot)
saveRDS(perRes,file=paste0(sl_input$out_path,"/delta",sl_input$delta,"lambda",sl_input$lambda,"_boxplot_data.rds"))
}
