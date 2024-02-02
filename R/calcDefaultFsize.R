#' calculate the default f_size
#' @param y a matrix of the response vector
#' @param all_latent_factors a object contain all latent factors from getLatentFactors.
#' @return a integer representing the default f_size.
#' @export


calcDefaultFsize <- function(y, all_latent_factors){
  
  if ((dim(y)[[1]] <= all_latent_factors$K) && (all_latent_factors$K < 100)){
    if (abs(dim(y)[[1]]-all_latent_factors$K) <= 2){
      f_size = dim(y)[[1]] - 2
    } else{
      f_size = dim(y)[[1]]
    }
  }
  
  if ((dim(y)[[1]] > all_latent_factors$K) && (all_latent_factors$K < 100)){
    f_size = all_latent_factors$K
  }

  if (all_latent_factors$K >= 100){
    f_size = 100
  }
return(f_size)
}