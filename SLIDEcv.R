require(doParallel)
require(dplyr)
library(SLIDE)
library(EssReg)
library(yaml)


args <- commandArgs(trailingOnly = T)
yaml_path  <- args[1]
#replicate  <- args[2]
sprintf("The ymal path is:  %s",yaml_path)
er_input <- yaml::yaml.load_file(as.character(yaml_path))




for(j in 1:er_input$nreps){

  paperBench(yaml_path, replicate=j)

}

er_input <- yaml::yaml.load_file(yaml_path)

pathLists <- list.files(er_input$out_path,recursive = T,pattern = "results")
perfList <- lapply(paste0(er_input$out_path,pathLists), readRDS)
perRes <- do.call(rbind,lapply(perfList,function(x){x$final_corr}))


if (er_input$eval_type == "corr") {
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
ggsave(paste0(er_input$out_path,"/delta",er_input$delta,"lambda",er_input$lambda,"_boxplot.pdf"),lambda_boxplot)
