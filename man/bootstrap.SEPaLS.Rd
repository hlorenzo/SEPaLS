% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEPaLS.R
\name{bootstrap.SEPaLS}
\alias{bootstrap.SEPaLS}
\title{Bootstrap function for SEPaLS estimator.}
\usage{
bootstrap.SEPaLS(
  X,
  Y,
  yn,
  type = c("vMF", "Laplace"),
  mu0 = NULL,
  kappa0 = NULL,
  lambda = NULL,
  B = 20
)
}
\arguments{
\item{X}{\eqn{(n\times p)}-dimensional matrix of the covariates.}

\item{Y}{\eqn{(n)}-dimensional vector of the response.}

\item{yn}{\eqn{y_n} the quantile corresponding to lowest values of \eqn{Y}s
to put in the tail.}

\item{type}{character, whether \code{vMF} for von Mises-Fisher prior or
\code{Laplace} for Laplace prior. See details.}

\item{mu0}{\eqn{\mu_0}, unitary \eqn{(p)}-dimensional vector. The direction
parameter for the \code{vMF} prior.}

\item{kappa0}{\eqn{\kappa_0}, positive. The concentration parameter for the
\code{vMF} prior.}

\item{lambda}{\eqn{\lambda}, positive. The concentration parameter for the
\code{Laplace} prior.}

\item{B}{positive integer. The number of bootstrap samples on which estimate
the SEPaLS directions. Default to 20.}
}
\value{
A list with two elements:
\itemize{
    \item \code{ws}: A \eqn{(B\times p)}-dimensional matrix  with each
    row corresponding to the \emph{SEPaLS} direction estimated on each
    bootstrap sample.
    \item \code{cor}: The correlation of each estimate direction on the
    Out-Of-Bag (OOB) sample with the response.
}
}
\description{
Bootstrap function for SEPaLS estimator.
}
\examples{
set.seed(5)
n <- 3000
p <- 10
X <- matrix(rnorm(n*p),n,p)
beta <- c(5:1,rep(0,p-5)) ; beta <- beta/sqrt(sum(beta^2))
Y <- (X\%*\%beta)^3 + rnorm(n)
boot.sepals_Laplace <- bootstrap.SEPaLS(X,Y,yn=1,type="Laplace",lambda=0.01,
B=100)
boxplot(boot.sepals_Laplace$ws);abline(h=0,col="red",lty=2)
}
\seealso{
\code{\link{SEPaLS}}
}
