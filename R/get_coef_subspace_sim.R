
#' @title Obtain the variance under each lambda by the cofficient matrices.
#' @description Defination is shown on the introduction.
#' @param coef_res a dataframe of dimension(times,lambdas) and each element is a coefficient matrix. 
#' @return the defined variance under each lambda.
#' @examples
#' \dontrun{
#' coef <- list()
#' for(t.j in 1:10){
#'   coef[[t.j]] <- list()
#'   for(l.i in 1:5){
#'     coef[[t.j]][[l.i]] <- matrix(rnorm(30),nrow=6)
#'   }
#' }
#' res <- get_coef_subspace_sim(coef)
#' }
#' @export
get_coef_subspace_sim <- function(coef_res){
  nlambdas <- length(coef_res[[1]])
  ntimes <- length(coef_res)
  var.subspace.coef <- rep(0, nlambdas)
  n <- min(dim(coef_res[[1]][[1]]))
  for(lambda.i in 1:nlambdas){
    Bj <- list()
    Bjv <- list()
    Bjb <- list()
    for (time.j in 1:ntimes) {
      Bj[[time.j]] <- coef_res[[time.j]][[lambda.i]]
      r <- qr(Bj[[time.j]])$rank
      Bjv[[time.j]] <- qr.Q(qr(Bj[[time.j]]))[,(r+1):n]
      Bjb[[time.j]] <- qr.Q(qr(Bj[[time.j]]))[,1:r]
    }
    for (time.j1 in 2:ntimes) {
      for (time.j2 in 1:(time.j1-1)) {
        s1 <- svd(t(Bjv[[time.j1]])%*%Bjb[[time.j2]])$d[1]
        s2 <- svd(t(Bjv[[time.j2]])%*%Bjb[[time.j1]])$d[1]
        var.subspace.coef[lambda.i] <- var.subspace.coef[lambda.i]+max(s1,s2)
      }
    }
  }
  var.subspace.coef <- var.subspace.coef*2/(ntimes*(ntimes-1))
  return(var.subspace.coef)
}