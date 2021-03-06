% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_coef_subspace_sim.R
\name{get_coef_subspace_sim}
\alias{get_coef_subspace_sim}
\title{Obtain the variance under each lambda by the cofficient matrices.}
\usage{
get_coef_subspace_sim(coef_res)
}
\arguments{
\item{coef_res}{a dataframe of dimension(times,lambdas) and each element is a coefficient matrix.}
}
\value{
the defined variance under each lambda.
}
\description{
Defination is shown on the introduction.
}
\examples{
\dontrun{
coef <- list()
for(t.j in 1:10){
  coef[[t.j]] <- list()
  for(l.i in 1:5){
    coef[[t.j]][[l.i]] <- matrix(rnorm(30),nrow=6)
  }
}
res <- get_coef_subspace_sim(coef)
}
}
