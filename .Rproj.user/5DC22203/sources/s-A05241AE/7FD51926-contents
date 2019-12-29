
#' @title Auto generated samples.
#' @description Simulate according to Chen, Dong, Chan for case where p,q >n.
#' @param n samples
#' @param p number of explanatory variables
#' @param q number of response variables
#' @param rho correlation between Xs
#' @param rX rank of X
#' @param r rank of C
#' @param b signal to noise ratio
#' @return a list of Y and X and true signal to noise ratio
#' @examples
#' \dontrun{
#' n = 76 # samples
#' q = 88 # number of response variables
#' rho = 0.5 # correlation between Xs
#' rX = 10 # rank of X
#' r = 5 # rank of C
#' p = 46 # number of explanatory variables
#' b = 1 #signal to noise ratio
#' res <- rrr.sim6(n,p,q,r,b,rho,rX)
#' }
#' @export
rrr.sim6 <- function(n, p, q, r, b, rho, rX){
  library("MASS")
  # Generate C
  C1 <- matrix(rnorm(p*r), nrow = p, ncol = r)
  C2 <- matrix(rnorm(q*r), nrow = q, ncol = r)
  C <- b*C1 %*% t(C2)
  # Generate X
  gamma <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p){
    for (j in 1:p){
      gamma[i,j] <- rho^(abs(i-j))
    }
  }
  gamma.eig <- eigen(gamma)
  gamma_sqrt <- gamma.eig$vectors %*% diag(sqrt(gamma.eig$values)) %*% solve(gamma.eig$vectors)
  X1 <- matrix(rnorm(n*rX), nrow = n, ncol = rX)
  X2 <- matrix(rnorm(rX*p), nrow = rX, ncol = p)
  X0 <- X1 %*% X2
  X <- X0 %*% gamma_sqrt
  # Generate Y
  E <- matrix(rnorm(n*q), nrow = n, ncol = q)
  Y <- X %*% C + E
  # Signal to noise ratio
  svdxc <- svd(X%*%C)
  P <- X %*% ginv(t(X) %*% X) %*% t(X)
  svdpe <- svd(P %*% E)
  s2n <- svdxc$d[r]/svdpe$d[1]
  return(list(Y = Y, X = X, C = C, s2n = s2n))
}