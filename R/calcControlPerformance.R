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


calcControlPerformance <- function(z_matrix, y, do_interacts, SLIDE_res, condition, out_path){

  colnames(y) <- "y"
  zz <- row.names(z_matrix)
  if (sum(row.names(y) != row.names(z_matrix)) != 0){
    cat("The row names of y does not match the row names of the z matrix. Matching the names now...")
    row.names(y) <- row.names(z_matrix)
  }
  y <- y[zz,,drop=F]


  # Real ##########################################################################

  cat("Getting performance for real model... \n")

  sigK <- SLIDE_res$SLIDE_res$marginal_vars
  sigK<- toupper(sigK)
  sigK <- as.numeric(gsub("Z","",sigK))
  sigIn <- as.vector(SLIDE_res$SLIDE_res$interaction_vars)


  ## If we have interaction
  if(do_interacts){

    Fullreport    <- NULL
    Partialreport <- NULL

## Make the real model and get the performance

    if(length(sigIn)!=0){

    IntData   <- pairwiseInteractions(as.numeric(sigK),z_matrix)

    Dataint   <- IntData$interaction[, sigIn,drop=F]

    Data_real <- data.frame(y = y, z_matrix[, sigK,drop=F], Dataint)

    }else{


    Data_real <- data.frame(y = y, z_matrix[, sigK,drop=F])


    }

    if (condition == 'auc'){
      lmod    <- lm(y~.,data=Data_real)
      yhat    <- predict(lmod,Data_real[,-1,drop=F],type = 'response')
      perreal <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat), levels = sort(unique(y[, 1])), quite=T)
    }else if (condition == 'corr'){
      SumReal <- summary(lm(y ~ ., data = Data_real))
      perreal <- sqrt(SumReal$r.squared)
    }



    ## Make Full Random
    cat("Getting performance for fully random model... \n")

    if(length(sigIn)!=0){

    for (i in 1:1000) {


      sigKRandom      <- sample(ncol(z_matrix), size = length(sigK))
      IntDataRandom   <- pairwiseInteractions(sigKRandom, z_matrix)
      sigInRandom     <- sample(ncol(IntDataRandom$interaction), size = length(sigIn)) ## Random interaction
      IntDataRandom   <- IntDataRandom$interaction[, sigInRandom,drop=F]


      ## Make the final simulated  data
      Data_fullRandom <- data.frame(y = y, z_matrix[, sigKRandom,drop=F], IntDataRandom)

      ## Do performance check
      SumInt <- summary(lm(y ~ ., data = Data_fullRandom))

       if (condition == 'auc'){
        lmod  <- lm(y~.,data=Data_fullRandom)
        yhat <- predict(lmod,Data_fullRandom[,-1,drop=F],type = 'response')
        prerandom <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat), levels = sort(unique(y[, 1])), quite=T)
        Fullreport <- rbind(Fullreport, prerandom)

        }else if (condition == "corr"){

        SumInt$r.squared
        Fullreport <- rbind(Fullreport, sqrt(SumInt$r.squared))

      }}

  cat("Getting performance for partially random model... \n")

    for (i in 1:1000) {

      IntData  <- pairwiseInteractions(sigK, z_matrix)
      sigInRandom <- sample(ncol(IntData$interaction), size = length(sigIn)) ## Ranodm interaction
      IntDataRandom <- IntData$interaction[, sigInRandom,drop=F]
      Data_partialRandom <- data.frame(y = y, z_matrix[, sigK,drop=F], IntDataRandom)

      SumPInt <- summary(lm(y ~ ., data = Data_partialRandom))

      if (condition == 'auc'){
        lmod  <- lm(y~.,data=Data_partialRandom)
        yhat <- predict(lmod,Data_partialRandom[,-1,drop=F],type = 'response')
        perPrandom <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat),quite=T)
        Partialreport <- rbind(Partialreport, perPrandom)
      }else if(condition == 'corr'){

        SumPInt$r.squared
        Partialreport <- rbind(Partialreport, sqrt(SumPInt$r.squared))

      }

      }



      rawdf <- data.frame(FullRandom = Fullreport,PartialRandom=Partialreport)
      df <- reshape2::melt(rawdf)
      colnames(df) <- c("group", "value")

      cols <- c("#0000FF", "#00FF00")


      P2 <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) + ggplot2::geom_density(alpha = 0.7) + ggplot2::scale_fill_manual(values = cols) +
        ggplot2::theme_light()  +
        ggplot2::ylab("Density") +
        ggplot2::xlab(condition)

      P2 <- P2 + ggplot2::annotate("text", x = perreal + 0.01, y = 55, label = " ", angle = "90") + ggplot2::xlim(0.25, max(df$value) + 0.2)+
        ggplot2::geom_vline(xintercept = perreal, linetype = "longdash", colour = "red",size=2)

      P2 <- P2 + ggplot2::theme(panel.border = ggplot2::element_blank(),
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank(),
                                panel.background = ggplot2::element_blank(),
                                axis.line = ggplot2::element_line(colour = "black"),
                                axis.text = ggplot2::element_text(size = 20),
                                axis.title.x = ggplot2::element_text(size = 20),
                                axis.title.y = ggplot2::element_text(size = 20))
      P2 <- P2 + ggplot2::annotate("text", x = perreal + 0.01, y = 55, label = " ", angle = "90") + ggplot2::xlim(0.25, max(df$value) + 0.2)+
        ggplot2::geom_vline(xintercept = perreal, linetype = "longdash", colour = "red",size=2)


      P2
      saveRDS(df, file = paste0(out_path, "ControlPerformance.rds"))
      saveRDS(P2, file = paste0(out_path, "ControlPerformancePlot.rds"))
      ggplot2::ggsave(P2, filename = paste0(out_path, "ControlPerformancePlot.png"))

      }else{

        Data_real <- data.frame(y = y, z_matrix[, sigK,drop=F])

        if (condition == 'auc'){
          lmod  <- lm(y~.,data=Data_real)
          yhat <- predict(lmod,Data_real[,-1,drop=F],type = 'response')
          perreal <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat),quite=T)
        }else if (condition == 'corr'){
          SumReal <- summary(lm(y ~ ., data = Data_real))
          perreal <- sqrt(SumReal$r.squared)
        }


        ## Generating the random
        for (i in 1:1000) {

          sigKRandom      <- sample(ncol(z_matrix), size = length(sigK))

          ## Make the final simulated  data

          Data_fullRandom <- data.frame(y = y, z_matrix[, sigKRandom,drop=F])

          ## Do performance check

          SumInt <- summary(lm(y ~ ., data = Data_fullRandom))

          if (condition == 'auc'){

            lmod  <- lm(y~.,data=Data_fullRandom)
            yhat <- predict(lmod,Data_fullRandom[,-1,drop=F],type = 'response')
            prerandom <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat),quite=T)
            Fullreport <- rbind(Fullreport, prerandom)

          }else if (condition == "corr"){

            SumInt$r.squared
            Fullreport <- rbind(Fullreport, sqrt(SumInt$r.squared))

          }
        }

        rawdf <- data.frame(FullRandom = Fullreport)
        df <- reshape2::melt(rawdf)
        colnames(df) <- c("group", "value")

        cols <- c("#0000FF")


        P2 <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) + ggplot2::geom_density(alpha = 0.7) + ggplot2::scale_fill_manual(values = cols) +
          ggplot2::theme_light()  +
          ggplot2::ylab("Density") +
          ggplot2::xlab(condition)

        P2 <- P2 + ggplot2::annotate("text", x = perreal + 0.01, y = 55, label = " ", angle = "90") + ggplot2::xlim(0.25, max(df$value) + 0.2)+
          ggplot2::geom_vline(xintercept = perreal, linetype = "longdash", colour = "red",size=2)

        P2 <- P2 + ggplot2::theme(panel.border = ggplot2::element_blank(),
                                  panel.grid.major = ggplot2::element_blank(),
                                  panel.grid.minor = ggplot2::element_blank(),
                                  panel.background = ggplot2::element_blank(),
                                  axis.line = ggplot2::element_line(colour = "black"),
                                  axis.text = ggplot2::element_text(size = 20),
                                  axis.title.x = ggplot2::element_text(size = 20),
                                  axis.title.y = ggplot2::element_text(size = 20))
        P2 <- P2 + ggplot2::annotate("text", x = perreal + 0.01, y = 55, label = " ", angle = "90") + ggplot2::xlim(0.25, max(df$value) + 0.2)+
          ggplot2::geom_vline(xintercept = perreal, linetype = "longdash", colour = "red",size=2)


        P2
        saveRDS(df, file = paste0(out_path, "ControlPerformance.rds"))
        saveRDS(P2, file = paste0(out_path, "ControlPerformancePlot.rds"))
        ggplot2::ggsave(P2, filename = paste0(out_path, "ControlPerfomancePlot.png"))

        }


   }else{

      ## No interaction flag

      ## Getting real performance
      Data_real <- data.frame(y = y, z_matrix[, sigK,drop=F])

      if (condition == 'auc'){
        lmod  <- lm(y~.,data=Data_real)
        yhat <- predict(lmod,Data_real[,-1,drop=F],type = 'response')
        perreal <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat),quite=T)
      }else if (condition == 'corr'){
        SumReal <- summary(lm(y ~ ., data = Data_real))
        perreal <- sqrt(SumReal$r.squared)
      }

      Fullreport <- NULL
      ## Generating the random
      for (i in 1:1000) {

        sigKRandom      <- sample(ncol(z_matrix), size = length(sigK))

        ## Make the final simulated  data

        Data_fullRandom <- data.frame(y = y, z_matrix[, sigKRandom,drop=F])

        ## Do performance check

        SumInt <- summary(lm(y ~ ., data = Data_fullRandom))

        if (condition == 'auc'){

          lmod  <- lm(y~.,data=Data_fullRandom)
          yhat <- predict(lmod,Data_fullRandom[,-1],type = 'response')
          prerandom <- pROC::auc(response=as.matrix(y), predictor=as.matrix(yhat),quite=T)
          Fullreport <- rbind(Fullreport, prerandom)

        }else if (condition == "corr"){

          SumInt$r.squared
          Fullreport <- rbind(Fullreport, sqrt(SumInt$r.squared))

        }
        }

      rawdf <- data.frame(FullRandom = Fullreport)
      df <- reshape2::melt(rawdf)
      colnames(df) <- c("group", "value")

      cols <- c("#0000FF")


      P2 <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) + ggplot2::geom_density(alpha = 0.7) + ggplot2::scale_fill_manual(values = cols) +
        ggplot2::theme_light()  +
        ggplot2::ylab("Density") +
        ggplot2::xlab(condition)

      P2 <- P2 + ggplot2::annotate("text", x = perreal + 0.01, y = 55, label = " ", angle = "90") + ggplot2::xlim(0, max(df$value) + 0.2)+
        ggplot2::geom_vline(xintercept = perreal, linetype = "longdash", colour = "red",size=2)

      P2 <- P2 + ggplot2::theme(panel.border = ggplot2::element_blank(),
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank(),
                                panel.background = ggplot2::element_blank(),
                                axis.line = ggplot2::element_line(colour = "black"),
                                axis.text = ggplot2::element_text(size = 20),
                                axis.title.x = ggplot2::element_text(size = 20),
                                axis.title.y = ggplot2::element_text(size = 20))
      P2 <- P2 + ggplot2::annotate("text", x = perreal + 0.01, y = 55, label = " ", angle = "90") + ggplot2::xlim(0, max(df$value) + 0.2)+
        ggplot2::geom_vline(xintercept = perreal, linetype = "longdash", colour = "red",size=2)







    saveRDS(df, file = paste0(out_path, "ControlPerformance.rds"))
    saveRDS(P2, file = paste0(out_path, "ControlPerformancePlot.rds"))
    ggplot2::ggsave(P2, filename = paste0(out_path, "ControlPerfomancePlot.pdf"))
    P2

   }
}




