
#' @title RRR via stability approach.
#' @description The model is described in the introductioon.
#' @param responce the responce matrix.
#' @param lambda the path of lambda.
#' @param subset the proportion of subset from all samples.
#' @param ntimes the number of subsamples.
#' @return a list of coefficient matrix and stabilities and lambda.
#' @examples
#' \dontrun{
#' library(SC19074)
#' n = 76 # samples
#' q = 88 # number of response variables
#' rho = 0.5 # correlation between Xs
#' rX = 10 # rank of X
#' r = 5 # rank of C
#' p = 46 # number of explanatory variables
#' sim21 <- rrr.sim6(n, p, q, r, b = 0.5, rho, rX)
#' sim_ann21 <- rrr2(sim21$Y, sim21$X, penaltySVD = "ann",modstr = list(gamma = 2))
#' sim_res21 <- stability_ann_sim(sim21$Y, sim21$X, sim_ann21$lambda, subset = 0.7, ntimes = 100)
#' }
#' @export
stability_ann_sim <- function(response, X, lambda, subset = 0.7, ntimes = 50){
  # Function to run stability analysis for ANN reduced rank regression on simulated data
  # Stability of rank for each lambda
  # Rows are subsets, columns are lambdas
  st_res <- matrix(NA, nrow = ntimes, ncol = length(lambda))
  # List: first index is over subsamples, second index over lambdas
  coef_res <- list()
  # Number to sample each time
  nsample <- ceiling(subset*nrow(X))
  for(i in 1:ntimes){
    # Rows to sample
    sample_i <- sample.int(nrow(X), size = nsample, replace = FALSE)
    y_sub <- response[sample_i,]
    x_sub <- X[sample_i,]
    sub_res <- list()
    for(j in 1:length(lambda)){
      mod <- rrr2(as.matrix(y_sub), as.matrix(x_sub), penaltySVD = "ann",
                  modstr = list(gamma = 2, lambda = lambda[j]))
      st_res[i,j] <- mod$rank  
      sub_res[[j]] <- mod$coef
    }
    coef_res[[i]] <- sub_res
  }
  # Stability of coefficients for each lambda
  st_coef <-data.frame(lambda = lambda, var.subspace.coef = get_coef_subspace_sim(coef_res))
  return(list(st_res = st_res, st_coef = st_coef, lambda = lambda))
}