makePosDef <- function(sigma) {
  eigen_res <- eigen(sigma)
  new_eigen_value <- ifelse(eigen_res$values < 10e-6, 10e-4, eigen_res$values)
  sigma2 <- eigen_res$vectors %*% diag(x = new_eigen_value, nrow = length(new_eigen_value), ncol = length(new_eigen_value)) %*% t(eigen_res$vectors)
  sigma2 <- (sigma2 + t(sigma2)) / 2
  return(sigma2)
}
