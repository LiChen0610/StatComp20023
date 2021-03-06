% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp20023R.R
\name{simu_sirs_generation}
\alias{simu_sirs_generation}
\title{Function generating data for simulation 1 in SICA paper}
\usage{
simu_sirs_generation(s, n, p, r, beta)
}
\arguments{
\item{s}{the number of support}

\item{n}{nrow in design matrix}

\item{p}{ncol in design matrix}

\item{r}{the correlation coefficients}

\item{beta}{the true coefficients of the model}
}
\value{
the design matrix as well as its response vector as a list
}
\description{
Function that generates data for simulation 1 in the paper of Lv and Fan
}
\examples{
\dontrun{
s <- 7; n<- 35; p <- 100; r <- 0.1; beta <- c(1,-0.5,0.7,-1.2,-0.9,3,0.55);
data <- simu_sirs_generation(s,n,p,r,beta)
}
}
\references{
Lv, Jinchi; Fan, Yingying. A unified approach to model selection and sparse recovery using regularized least squares. Ann. Statist. 37 (2009), no. 6A, 3498--3528. doi:10.1214/09-AOS683.
}
\seealso{
\code{\link{sirs}}
}
