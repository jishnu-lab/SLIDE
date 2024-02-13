#' calcControlPerformance
#'
#' Calculate control model performances including the real model, the model with random interaction terms and full random model.
#'
#' @param z_matrix the z_matrix, the output of CalcZMatrix function.
#' @param y the response matrix.
#' @param SLIDE_res the output of GetTopFeatures function.
#' @param niter a integer for number of iterations.
#' @param condition a string on whether to run auc or corr.
#' @param out_path the output path.
#' @return The z matrix calculated by the CalcZMatrix function
#' @export


calcControlPerformance <- function(z_matrix, y, SLIDE_res, niter, condition, out_path){

  colnames(y) <- "y"
  zz <- row.names(z_matrix)
  if (sum(row.names(y) != row.names(z_matrix)) != 0){
    cat("The row names of y does not match the row names of the z matrix. Matching the names now...")
    row.names(y) <- row.names(z_matrix)
  }
  y <- y[zz,,drop=F]


  # Real ##########################################################################

  sigK <- SLIDE_res$SLIDE_res$marginal_vars
  sigK<- toupper(sigK)
  sigK <- as.numeric(gsub("Z","",sigK))

  sigIn <- as.vector(SLIDE_res$SLIDE_res$interaction_vars)



  if(!is.null(sigIn)){
  IntData <- pairwiseInteractions(as.numeric(sigK),z_matrix)
  Dataint <- IntData$interaction[, sigIn]
  Data_real <- data.frame(y = y, z_matrix[, sigK], Dataint)

  }else{

  Data_real <- data.frame(y = y, z_matrix[, sigK])

  }

  if (condition == 'auc'){
    lmod  <- lm(y~.,data=Data_real)
    yhat <- predict(lmod,Data_real[,-1],type = 'response')
    aucreal <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat))
  }else if (condition == 'corr'){
    SumReal <- summary(lm(y ~ ., data = Data_real))
    aucreal <- sqrt(SumReal$r.squared)
  }


  # All random ####################################################################
  Fullreport <- NULL
  Partialreport <- NULL

  for (i in 1:niter) {
    sigKRandom <- sample(ncol(z_matrix), size = length(sigK)) ## Random marginal

    if(!is.null(sigIn)){
    IntDataRandom <- pairwiseInteractions(sigKRandom, z_matrix)

    sigInRandom <- sample(ncol(IntDataRandom$interaction), size = length(sigIn)) ## Random interaction
    IntDataRandom <- IntDataRandom$interaction[, sigInRandom]

    Data_fullRandom <- data.frame(y = y, z_matrix[, sigKRandom], IntDataRandom) ## With interaction

    }else{

      Data_fullRandom <- data.frame(y = y, z_matrix[, sigKRandom]) # No interaction

    }


    SumInt <- summary(lm(y ~ ., data = Data_fullRandom))
    if (condition == 'auc'){
      lmod  <- lm(y~.,data=Data_fullRandom)
      yhat <- predict(lmod,Data_fullRandom[,-1],type = 'response')
      aucrandom <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat))
      Fullreport <- rbind(Fullreport, aucrandom)
    }else if (condition == "corr"){
      SumInt$r.squared
      Fullreport <- rbind(Fullreport, sqrt(SumInt$r.squared))

    }

  }
  Partialreport <- NULL
  if(!is.null(sigIn)){

      for (i in 1:niter) {

        IntData  <- pairwiseInteractions(sigK, z_matrix)
        sigInRandom <- sample(ncol(IntData$interaction), size = length(sigIn)) ## Ranodm interaction

        IntDataRandom <- IntData$interaction[, sigInRandom]
        Data_partialRandom <- data.frame(y = y, z_matrix[, sigK], IntDataRandom)}
        SumPInt <- summary(lm(y ~ ., data = Data_partialRandom))

      if (condition == 'auc'){
        lmod  <- lm(y~.,data=Data_partialRandom)
        yhat <- predict(lmod,Data_partialRandom[,-1],type = 'response')
        aucPrandom <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat))
        Partialreport <- rbind(Partialreport, aucPrandom)
      }else if(condition == 'corr'){
        SumPInt$r.squared
        Partialreport <- rbind(Partialreport, sqrt(SumPInt$r.squared))
      }


        rawdf <- data.frame(FullRandom = Fullreport,Partialreport=Partialreport)
        df <- reshape2::melt(rawdf)
        colnames(df) <- c("group", "value")

        cols <- c("#0000FF", "#00FF00")


        P2 <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) + ggplot2::geom_density(alpha = 0.7) + ggplot2::scale_fill_manual(values = cols) +
          ggplot2::theme_light()  +
          ggplot2::ylab("Density") +
          ggplot2::xlab(condition)

        P2 <- P2 + ggplot2::annotate("text", x = aucreal + 0.01, y = 55, label = " ", angle = "90") + ggplot2::xlim(0.25, max(df$value) + 0.1)+
          ggplot2::geom_vline(xintercept = aucreal, linetype = "longdash", colour = "red",size=2)

        P2 <- P2 + ggplot2::theme(panel.border = ggplot2::element_blank(),
                                  panel.grid.major = ggplot2::element_blank(),
                                  panel.grid.minor = ggplot2::element_blank(),
                                  panel.background = ggplot2::element_blank(),
                                  axis.line = ggplot2::element_line(colour = "black"),
                                  axis.text = ggplot2::element_text(size = 20),
                                  axis.title.x = ggplot2::element_text(size = 20),
                                  axis.title.y = ggplot2::element_text(size = 20))
        P2 <- P2 + ggplot2::annotate("text", x = aucreal + 0.01, y = 55, label = " ", angle = "90") + ggplot2::xlim(0.25, max(df$value) + 0.1)+
          ggplot2::geom_vline(xintercept = aucreal, linetype = "longdash", colour = "red",size=2)

       P2
        }else{







      rawdf <- data.frame(FullRandom = Fullreport)
      df <- reshape2::melt(rawdf)
      colnames(df) <- c("group", "value")

      cols <- c("#0000FF")

      # Basic density plot in ggplot2

      P2 <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) + ggplot2::geom_density(alpha = 0.7) + ggplot2::scale_fill_manual(values = cols) +
        ggplot2::theme_light() + ggplot2::geom_vline(xintercept = aucreal, linetype = "longdash", colour = "red",size=2) +
        ggplot2::geom_vline(xintercept = aucreal, linetype = "longdash", colour = "red",size=2)+
        ggplot2::ylab("Density") +
        ggplot2::xlab(condition)

      P2 <- P2  + ggplot2::xlim(0.25, max(df$value) + 0.1)

      P2 <- P2 + ggplot2::theme(panel.border = ggplot2::element_blank(family='Arial'),
                                panel.grid.major = ggplot2::element_blank(family='Arial'),
                                panel.grid.minor = ggplot2::element_blank(family='Arial'),
                                panel.background = ggplot2::element_blank(family='Arial'),
                                axis.line = ggplot2::element_line(colour = "black"),
                                axis.text = ggplot2::element_text(size = 20,family = 'Arial'),
                                axis.title.x = ggplot2::element_text(size = 20,family = 'Arial'),
                                axis.title.y = ggplot2::element_text(size = 20,family = 'Arial'))
      P2





    }

    saveRDS(df, file = paste0(out_path, "ControlPerformance.rds"))
    saveRDS(P2, file = paste0(out_path, "ControlPerformancePlot.rds"))
    ggplot2::ggsave(P2, filename = paste0(out_path, "ControlPerfomancePlot.png"))


}




