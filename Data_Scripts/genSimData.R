genSimData <- function(n = 100, p = 300, k = 10, sig_k = 5, sig_inx = 0) {
  #### generate data such that there are k clusters ####
  remainder <- p %% k
  vars_per_cluster <- p %/% k
  cluster_sizes <- rep(vars_per_cluster, k)
  ## if there is a remainder, then add until the sum of the cluster sizes = p
  if (remainder > 0) {
    cluster_sizes[1:remainder] <- cluster_sizes[1:remainder] + 1
  }

  ## initialization
  cluster_indices <- list()
  curr_ind <- 1
  x <- matrix(nrow = n, ncol = 0)
  ## loop through the cluster sizes (basically looping through each cluster)
  for (i in 1:length(cluster_sizes)) {
    ## get cluster size
    cluster_size <- cluster_sizes[i]
    ## make the covariance matrix: 1s on diagonal, 0.8 everywhere else
    mvtnorm_sigma <- diag(x = 1, nrow = cluster_size, ncol = cluster_size)
    mvtnorm_sigma[mvtnorm_sigma == 0] <- 0.8
    ## generate cluster data (all Xs in this cluster have high correlation)
    cluster_data <- mvtnorm::rmvnorm(n = n, mean = rep(0, cluster_size), sigma = mvtnorm_sigma)
    ## append the generated data to the X matrix
    x <- cbind(x, cluster_data)
    ## keep track of which X indices correspond to this cluster
    cluster_indices[[i]] <- seq(from = curr_ind, to = (curr_ind + cluster_size - 1), by = 1)
    curr_ind <- curr_ind + cluster_size
  }
  ## rename columns of X matrix
  colnames(x) <- paste0("X", 1:ncol(x))

  #### generate coefficients ####
  ## sample sig_k clusters (these will be used to generate Y)
  nonzero_k <- sample(1:k, sig_k)
  ## randomly sample coefficient values for these sig_k clusters
  nonzero_beta_k <- rnorm(sig_k, mean = 5, sd = 0.5)
  ## initialization of all X coefficients as 0
  nonzero_beta_x <- rep(0, p)
  ## loop through the nonzero clusters
  for (i in 1:length(nonzero_k)) {
    ## extract the corresponding X indices
    x_indices <- cluster_indices[[nonzero_k[i]]]
    ## give these X's nonzero coefficients
    nonzero_beta_x[x_indices] <- nonzero_beta_k[i]
  }



  # If there is no interaction -----------------------------------------------

  if (sig_inx == 0) {
    marg_inds <- nonzero_k

    marginal_vals <- as.matrix(x[, marg_inds])


    beta_marginal <- sapply(1:ncol(marginal_vals),
      function(x) {
        rnorm(1, 1, 0.5)
      },
      simplify = T
    )


    full_data <- data.frame(marginal_vals)

    beta_marginal_corss <- apply(beta_marginal * marginal_vals, 1, sum)


    full_betas <- c(beta_marginal)
    names(full_betas) <- colnames(full_data)

    y <- 1 / sd(c(beta_marginal_corss)) * c(beta_marginal_corss) +
      stats::rnorm(n, sd = 1) ## add a little bit of noise

    er_res <- getLatentFactors(
      y = y,
      x = x,
      x_std = scale(x, T, T),
      std_y = T,
      sigma = NULL,
      delta = 0.1,
      lambda = 1,
      thresh_fdr = 0.2,
      rep_cv = 50,
      alpha_level = 0.05,
      out_path = NULL
    )



    z <- EssReg::predZ(x = scale(x, T, T), er_res = er_res)
    colnames(z) <- paste0("Z", 1:ncol(z))


    return(list(
      "x" = x,
      "y" = y,
      "z" = z,
      "marginals" = marginal_vals,
      "interactions" = NULL,
      "full_data" = full_data,
      "betas" = full_betas
    ))
  }


  # If there are interactions -------------------------------------------------

  else if (sig_inx > 0) {
    #### generate interaction terms ####
    ## sample sig_inx interactions (these will be used to generate)
    marg_inds <- nonzero_k
    ## initialization of interaction terms and their coefficients
    # inx <- matrix(nrow = n, ncol = 0)
    inx_betas <- c()
    ## single interaction term per marginal
    interaction_terms <- interUnion(marg_inds, x)
    inx_x_ind <- sample(size = sig_inx, 1:ncol(interaction_terms$interaction))
    interaction_vals <- interaction_terms$interaction[, inx_x_ind] ## sampling from interaction
    marginal_vals <- interaction_terms$marginal



    ## Beta for marginals
    beta_marginal <- sapply(1:ncol(marginal_vals),
      function(x) {
        rnorm(1, 1, 0.5)
      },
      simplify = T
    )



    ## Beta For interaction
    beta_interact <- sapply(1:ncol(interaction_vals),
      function(x) {
        rnorm(1, 1, 0.5)
      },
      simplify = T
    )




    beta_marginal_cross <- apply(beta_marginal * marginal_vals, 1, sum)

    beta_interact_cross <- apply(beta_interact * interaction_vals, 1, sum)

    ## generate Y from clusters and interaction terms

    full_data <- data.frame(marginal_vals, interaction_vals)

    full_betas <- c(beta_marginal, beta_interact)
    names(full_betas) <- colnames(full_data)

    y <- 1 / sd(c(beta_marginal_cross + beta_interact_cross)) * c(beta_marginal_cross + beta_interact_cross) +
      stats::rnorm(n, sd = 1) ## add a little bit of noise

    ## run ER to get clusters
    er_res <- EssReg::plainER(
      y = y,
      x = x,
      x_std = scale(x, T, T),
      std_y = T,
      sigma = NULL,
      delta = 0.1,
      lambda = 1,
      thresh_fdr = 0.2,
      rep_cv = 50,
      alpha_level = 0.05,
      out_path = NULL
    )
    ## generate Zs
    z <- EssReg::predZ(x = scale(x, T, T), er_res = er_res)
    colnames(z) <- paste0("Z", 1:ncol(z))


    return(list(
      "x" = x,
      "y" = y,
      "z" = z,
      "marginals" = marginal_vals,
      "interactions" = interaction_vals,
      "full_data" = full_data,
      "betas" = full_betas
    ))
  }
}
