#' Get Latent Factors from Input Data.
#'
#' Perform Essential Regression without prior knowledge.
#'
#' @importFrom foreach '%dopar%'
#' @param y a response vector of dimension \eqn{n}, must be continous
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param std_y a boolean flag of wether z score Y or not
#' @param x_std a data matrix of dimensions \eqn{n \times p}, a scaled and mean centered x matrix
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}, default is \code{NULL}
#' @param delta \eqn{\delta}, a numerical value used for thresholding, or a vector of values
#' @param thresh_fdr a numerical constant used for thresholding the correlation matrix to
#' control the false discovery rate, default is 0.2
#' @param lambda \eqn{\lambda}, a numerical constant used in thresholding
#' @param rep_cv number of replicates for cross-validation
#' @param alpha_level \eqn{\alpha}, a numerical constant used in confidence interval calculation
#' @param out_path a string path to where to save output
#' @return a list of results from the Essential Regression framework including: \eqn{K} = number of clusters,
#' \eqn{\hat{A}}, \eqn{\hat{C}}, \eqn{\hat{I}}, the indices of the pure variables, \eqn{\hat{\Gamma}},
#' \eqn{\hat{\beta}}, \eqn{\alpha}-level confidence intervals (if requested), prediction results (if requested),
#' the optimal value of \eqn{\lambda} determined by cross-validation, the optimal value of \eqn{\delta}
#' determined by cross-validation, \eqn{Q}, and the variances of \eqn{\hat{\beta}}
#' @export

getLatentFactors <- function(y, x, x_std, std_y = TRUE, sigma = NULL, delta, thresh_fdr = 0.2, lambda = 0.1,
                    rep_cv = 50, alpha_level = 0.05, out_path = NULL) {
  ## Data Housekeeping #########################################################
  raw_y <- y
  raw_x <- x  # input raw x is now raw_x

  ## standardization but only x
  x <- x_std # scaled x is now x

  if (std_y) {
    y <- scale(y, T, T)
  }


  opt_lambda <- lambda ## no longer CV lambda
  ## ER Housekeeping ###########################################################
  n <- nrow(x);  p <- ncol(x) #### feature matrix dimensions
  se_est <- apply(x, 2, stats::sd) #### get sd of columns for feature matrix

  if (is.null(sigma)) {
    sigma <- cor(x)
  }

  #### scale delta by this value to satisfy some requirements so that the
  #### statistical guarantees in the paper hold
  delta_scaled <- delta * sqrt(log(max(p, n)) / n)

  ## Sigma Thresholding ########################################################
  #### save correlation matrix heatmap
  if (!is.null(out_path)) {
    pdf_file <- paste0(out_path, "delta_", delta[1], "/corr_mat_heatmap.pdf")
    dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
    grDevices::pdf(file = pdf_file)
    makeHeatmap(sigma, "Correlation Matrix Heatmap", T, T)
    grDevices::dev.off()
  }

  #### threshold sigma to control for FDR
  if (!is.null(thresh_fdr)) {
    control_fdr <- threshSigma(x = x,
                               sigma = sigma,
                               thresh = thresh_fdr)
    sigma <- control_fdr$thresh_sigma
    kept_entries <- control_fdr$kept_entries
  } else {
    kept_entries <- matrix(1, nrow = nrow(sigma), ncol = ncol(sigma))
  }

  #### save threshold correlation matrix heatmap
  if (!is.null(out_path)) {
    pdf_file <- paste0(out_path, "delta_", delta[1], "/thresh_corr_mat_heatmap.pdf")
    dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
    grDevices::pdf(file = pdf_file)
    makeHeatmap(sigma, "FDR Thresholded Correlation Matrix Heatmap", T, T)
    grDevices::dev.off()
  }

  ## Delta Cross-Validation ####################################################
  #### if delta has more than 1 element, then do rep_CV # of replicates
  #### of CV_Delta and select median of replicates
  #### use the unstandardized version of x (if available) to avoid signal leakage in CV
  if (length(delta_scaled) > 1) {
    foreach::foreach(i = 1:rep_cv, .combine = c) %dopar% {
      cvDelta(raw_x = raw_x,
              fdr_entries = kept_entries,
              deltas_scaled = delta_scaled)
    } -> cv_delta_reps
    opt_delta <- stats::median(cv_delta_reps)
  } else {
    opt_delta <- delta_scaled
  }

  #### estimate membership matrix Ai
  #### also returns a vector of the indices of estimated pure variables
  #### and a list of the indices of estimated pure variables
  result_AI <- estAI(sigma = sigma, delta = opt_delta, se_est = se_est)
  pure_numb <- sapply(result_AI$pure_list, FUN = function(x) {length(c(x$pos, x$neg))})

  A_hat <- result_AI$AI
  I_hat <- result_AI$pure_vec
  I_hat_list <- result_AI$pure_list

  if (is.null(I_hat)) {
    cat("Algorithm fails due to non-existence of pure variable.\n")
    stop()
  }

  C_hat <- estC(sigma = sigma, AI = A_hat)

  #### ER Supplement (2.1)
  #### Gamma_hat_{ii} = Sigma_hat_{ii} - A_hat_{i.}^T Sigma_hat_{Z} A_hat_{i.} for i in I_hat
  #### Gamma_hat_{ji} = 0 for all j ≠ i
  
  Gamma_hat <- rep(0, p)
  Gamma_hat[I_hat] <- diag(sigma[I_hat, I_hat]) - diag(A_hat[I_hat, ] %*% C_hat %*% t(A_hat[I_hat, ]))
  Gamma_hat[Gamma_hat < 0] <- 1e-2 #### replace negative values with something very close to 0

  #### Prediction --- what is Q??
  #### use standardized x and y
  pred_result <- prediction(y = y, x = x, sigma = sigma, A_hat = A_hat,
                            Gamma_hat = Gamma_hat, I_hat = I_hat)

  #### theta_hat (supplement 2.2)
  theta_hat <- pred_result$theta_hat
  #### Inference In Latent Factor Regression With Clusterable Features
  #### Z_tilde = Q*X

  Q <- try(theta_hat %*% solve(crossprod(x %*% theta_hat) / n, crossprod(theta_hat)), silent = T)
  if (class(Q)[1] == "try-error") {
    Q <- theta_hat %*% MASS::ginv(crossprod(x %*% theta_hat) / n) %*% crossprod(theta_hat)
  }

  #### Beta Estimation
  if (length(result_AI$pure_vec) != nrow(sigma)) { ## check if all vars are pure vars?
    sigma_TJ <- estSigmaTJ(sigma = sigma, AI = A_hat, pure_vec = result_AI$pure_vec)
    # opt_lambda <- ifelse(length(lambda) > 1,
    #                      stats::median(unlist(replicate(rep_cv,
    #                                                     cvLambda(x = x,
    #                                                              fdr_entries = kept_entries,
    #                                                              lambdas = lambda,
    #                                                              AI = result_AI$AI,
    #                                                              pure_vec = result_AI$pure_ec),
    #                                                     simplify = F))),
    #                      lambda)
    opt_lambda <- lambda
    if (opt_lambda > 0) {
      AI_hat <- abs(A_hat[I_hat, ]) ## just rows of pure variables
      sigma_bar_sup <- max(solve(crossprod(AI_hat), t(AI_hat)) %*% se_est[I_hat]) ## not sure what this does
      AJ <- estAJDant(C_hat = C_hat, sigma_TJ = sigma_TJ,
                      lambda = opt_lambda * opt_delta * sigma_bar_sup,
                      se_est_J = sigma_bar_sup + se_est[-I_hat])
      if (is.null(AJ)) {
        return (NULL)
      }
    } else {
      AJ <- t(solve(C_hat, sigma_TJ))
    }
    A_hat[-result_AI$pure_vec, ] <- AJ
  }

  Gamma_hat[-I_hat] <- diag(sigma[-I_hat, -I_hat]) - diag(A_hat[-I_hat,] %*% C_hat %*% t(A_hat[-I_hat, ]))
  Gamma_hat[Gamma_hat < 0] <- 1e2 #### replace negative values with 100
  #### use standardized x and y
  
  res_beta <- estBeta(y = y, x = x, sigma = sigma, A_hat = A_hat,
                      C_hat = C_hat, Gamma_hat = Gamma_hat, I_hat = I_hat,
                      I_hat_list = I_hat_list, alpha_level = alpha_level)
  beta_hat <- res_beta$beta_hat
  beta_conf_int <- res_beta$conf_int
  beta_var <- res_beta$beta_var

  #### organization/renaming
  rownames(A_hat) <- colnames(x)
  colnames(A_hat) <- paste0("Z", 1:ncol(A_hat))

  rownames(C_hat) <- colnames(C_hat) <- paste0("Z", 1:ncol(C_hat))

  I_clust <- NULL
  for (i in 1:length(I_hat_list)) {
    clust_name <- paste0("Z", i)
    cluster <- I_hat_list[[i]]
    pos <- cluster$pos
    neg <- cluster$neg
    if (length(pos) > 0) {
      names(pos) <- colnames(x)[pos]
    } else {
      pos <- NULL
    }
    if (length(neg) > 0) {
      names(neg) <- colnames(x)[neg]
    } else {
      neg <- NULL
    }
    I_clust[[clust_name]] <- list("pos" = pos,
                                  "neg" = neg)
  }

  names(I_hat) <- colnames(x)[I_hat]

  return(list(K = ncol(A_hat),
              A = A_hat,
              C = C_hat,
              I = I_hat,
              I_clust = I_clust, ## original is I_ind
              Gamma = Gamma_hat,
              beta = beta_hat,
              beta_conf_int = beta_conf_int, ## original is beta_CIs
              beta_var = beta_var,
              pred = pred_result,
              opt_lambda = opt_lambda,
              opt_delta = opt_delta / sqrt(log(max(p, n)) / n),
              Q = Q,
              thresh_sigma = sigma))
}
