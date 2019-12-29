#' @title inner function 1.
#' @description null.
#' @param sv.tol singular value tolerence
#' @param qr.tol qr decomposition tolerence
#' @return a list.
#' @examples
#' \dontrun{
#' null
#' }
#' @export
rrr.control <- function(sv.tol = 1e-7, qr.tol = 1e-7)
{
  list(sv.tol = sv.tol,               # singular value tolerence
       qr.tol = qr.tol)               # QR decomposition tolerence
}

#' @title inner function 2 to give The default value.
#' @description null.
#' @param gamma inner para with default value 2
#' @param nlambda number of lambda
#' @param lambda inner para
#' @return a list.
#' @examples
#' \dontrun{
#' null
#' }
#' @export
rrr.modstr <- function(gamma = 2,
                       nlambda = 100,
                       lambda = NULL)
{
  list(gamma = gamma,
       nlambda = nlambda,
       lambda = lambda)
}



#' @title Soft thresholding
#' @description lambda can be a matrix with same dim as C.
#' @param c a inner para
#' @param Lam a inner para
#' @param tol tolerence
#' @return Soft thresholding
#' @examples
#' \dontrun{
#' null
#' }
#' @export
softTH <- function(C, Lam, tol = 1e-4) {
  Ct <- abs(C) - Lam
  Ct[Ct < tol] <- 0
  sign(C) * Ct
}


#' @title RRR model.
#' @description The model is developed by Kun Chen.
#' @param Y the resopnce matrix
#' @param X the design matrix
#' @param penaltySVD the penalty term
#' @param ic.type "GIC", "AIC", "BIC", "BICP" or "GCV"
#' @param df.type "exact" or "naive") 
#' @param maxrank the max rank to be expected
#' @param modstr the mod
#' @param control the control of SVD and QR
#' @return the solution of rrr.
#' @examples
#' \dontrun{
#' null
#' }
#' @export
rrr2 <- function (Y, X, penaltySVD = c("rank", "ann"), ic.type = c("GIC", 
                                                           "AIC", "BIC", "BICP", "GCV"), df.type = c("exact", "naive"), 
          maxrank = min(dim(Y), dim(X)), modstr = list(), control = list()) 
{
  Call <- match.call()
  penaltySVD <- match.arg(penaltySVD)
  ic.type <- match.arg(ic.type)
  df.type <- match.arg(df.type)
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  if (n != nrow(X)) 
    stop("'nrow(X)' has to be equal to 'nrow(Y)'.")
  control <- do.call("rrr.control", control)
  modstr <- do.call("rrr.modstr", modstr)
  qrX <- qr(X, tol = control$qr.tol)
  C_ls <- qr.coef(qrX, Y)
  C_ls <- ifelse(is.na(C_ls), 0, C_ls)
  rX <- qrX$rank
  nrank <- min(q, rX, maxrank)
  rmax <- min(rX, q)
  XC <- qr.fitted(qrX, Y)
  svdXC <- svd(XC, nu = rmax, nv = rmax)
  A <- svdXC$v[, 1:rmax]
  Ad <- (svdXC$d[1:rmax])^2
  Ad <- Ad[Ad > control$sv.tol]
  rmax <- length(Ad)
  if (identical(penaltySVD, "ann")) {
    gamma <- modstr$gamma
    nlambda <- modstr$nlambda
    lambda <- modstr$lambda
    Adx <- Ad^((1 + gamma)/2)
    if (is.null(lambda)) {
      lambda <- exp(seq(log(max(Adx)), log(Adx[min(nrank + 
                                                     1, rmax)]), length = nlambda))
    }
    f.d <- function(lambda, Ad) {
      Adx <- (Ad)^((1 + gamma)/2)
      softTH(Adx, lambda)/Adx
    }
    f.d.derv <- function(lambda, Ad) {
      lambda * (gamma + 1) * Ad^((-gamma - 2)/2)
    }
  }
  else if (identical(penaltySVD, "rank")) {
    lambda <- modstr$lambda
    if (is.null(lambda)) 
      lambda <- seq_len(nrank)
    f.d <- function(lambda, Ad) rep(1, lambda)
    f.d.derv <- function(lambda, Ad) rep(0, lambda)
  }
  nlam <- length(lambda)
  Spath <- matrix(0, nrank, nlam)
  df.exact <- df.naive <- rep(0, nlam)
  for (i in seq_len(nlam)) {
    f <- f.d(lambda[i], Ad[seq_len(nrank)])
    Spath[seq_along(f), i] <- f
    r <- sum(f > control$sv.tol)
    if (r >= 1) {
      if (identical(df.type, "exact")) {
        f.derv <- f.d.derv(lambda[i], Ad[seq_len(r)])
        term1 <- max(rX, q) * sum(f)
        a <- vector()
        count = 1
        for (k in seq_len(r)) {
          for (s in (r + 1):rmax) {
            a[count] <- (Ad[k] + Ad[s]) * f[k]/(Ad[k] - 
                                                  Ad[s])
            count <- count + 1
          }
        }
        term2 <- sum(a)
        if (r == rmax) 
          term2 <- 0
        if (r == rmax & r == min(p, q)) 
          term2 <- 0
        b <- vector()
        count = 1
        for (k in seq_len(r)) {
          for (s in seq_len(r)) {
            if (s == k) {
              b[count] <- 0
            }
            else {
              b[count] <- Ad[k] * (f[k] - f[s])/(Ad[k] - 
                                                   Ad[s])
            }
            count <- count + 1
          }
        }
        term3 <- sum(b)
        term4 <- sum(sqrt(Ad[seq_len(r)]) * f.derv)
        df.exact[i] <- term1 + term2 + term3 + term4
      }
      df.naive[i] <- r * (rX + q - r)
    }
  }
  tempFit <- X %*% C_ls
  rankall <- sse <- rep(0, nlam)
  for (i in seq_len(nlam)) {
    dsth <- Spath[, i]
    rank <- sum(dsth != 0)
    rankall[i] <- rank
    if (rank != 0) {
      tempC <- A[, seq_len(rank)] %*% (dsth[seq_len(rank)] * 
                                         t(A[, seq_len(rank)]))
      tempYhat <- tempFit %*% tempC
      sse[i] <- sum((Y - tempYhat)^2)
    }
    else {
      sse[i] <- sum(Y^2)
    }
  }
  logsse <- log(sse)
  df <- switch(df.type, exact = df.exact, naive = df.naive)
  ic <- switch(ic.type, GCV = n * q * sse/(n * q - df)^2, AIC = n * 
                 q * log(sse/n/q) + 2 * df, BIC = n * q * log(sse/n/q) + 
                 log(q * n) * df, BICP = n * q * log(sse/n/q) + 2 * log(p * 
                                                                          q) * df, GIC = n * q * log(sse/n/q) + log(log(n * q)) * 
                 log(p * q) * df)
  min.id <- which.min(ic)
  rankest <- rankall[min.id]
  dsth <- Spath[, min.id]
  if (rankest != 0) {
    U <- C_ls %*% A[, seq_len(rankest)] %*% diag(1/svdXC$d[seq_len(rankest)], 
                                                 nrow = rankest, ncol = rankest) * sqrt(n)
    D <- diag(svdXC$d[seq_len(rankest)] * dsth[seq_len(rankest)], 
              nrow = rankest, ncol = rankest)/sqrt(n)
    V <- A[, seq_len(rankest)]
    C <- U %*% D %*% t(V)
  }
  else {
    C <- matrix(nrow = p, ncol = q, 0)
    U <- matrix(nrow = p, ncol = 1, 0)
    V <- matrix(nrow = q, ncol = 1, 0)
    D <- 0
  }
  rval <- list(call = Call, Y = Y, X = X, A = A, Ad = Ad, coef.ls = C_ls, 
               Spath = Spath, df.exact = df.exact, df.naive = df.naive, 
               penaltySVD = penaltySVD, sse = sse, ic = ic, coef = C, 
               U = U, V = V, D = D, rank = rankest, lambda = lambda)
  class(rval) <- "rrr"
  rval
}