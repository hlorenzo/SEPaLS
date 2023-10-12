#' Maximum Likelihood estimator
#'
#' @param X \eqn{(n\times p)}-dimensional matrix of the covariates.
#' @param Y \eqn{(n)}-dimensional vector of the response.
#' @param yn the quantile corresponding to the lowest values of \eqn{Y}s to put
#' in the tail.
#'
#' @return The maximum likelihood estimator.
#' @export
#'
#' @examples
#' n <- 3000
#' p <- 10
#' X <- matrix(rnorm(n*p),n,p)
#' beta <- c(5:1,rep(0,p-5)) ; beta <- beta/sqrt(sum(beta^2))
#' Y <- X%*%beta + rnorm(n,sd=1/3)
#' estimators <- do.call(rbind,lapply(seq(0,1,length.out=100),function(pp){
#'   yn <- quantile(Y,probs = pp)
#'   maximum_Likelihood_SEPaLS(X,Y,yn)
#' }))
#' matplot(estimators,type="l",lty=1,col=c(rep(2,5),rep(1,p-5)))
#' abline(h=beta/sqrt(sum(beta^2)),col=c(rep(2,5),rep(1,p-5)))
maximum_Likelihood_SEPaLS <- function(X,Y,yn){
  id <- which(Y>yn)
  n <- length(Y)
  n_tail <- length(id)
  ## Select data in tail
  X_here <- X[id,,drop=FALSE]
  Y_here <- Y[id]
  ## Build statistics
  m_Y <- sum(Y_here)/n
  m_X <- colSums(X_here)/n
  m_XY <- as.vector(t(X_here)%*%Y_here)/n
  F_survival <- n_tail/n
  v_n <- F_survival*m_XY-m_X*m_Y
  return(v_n)
}

#' Function to estimate SEPaLS estimators
#'
#' @param X \eqn{(n\times p)}-dimensional matrix of the covariates.
#' @param Y \eqn{(n)}-dimensional vector of the response.
#' @param yn \eqn{y_n} the quantile correponding to lowest values of \eqn{Y}s to put in th
#' e tail.
#' @param type character, wether \code{vMF} for von Mises-Fisher prior or
#' \code{Laplace} for Laplace prior. See details.
#' @param mu0  \eqn{\mu_0}, unitary \eqn{(p)}-dimensional vector. The direction
#'  parameter for the \code{vMF} prior.
#' @param kappa0 \eqn{\kappa_0}, positive. The concentration parameter for the
#' \code{vMF} prior.
#' @param lambda \eqn{\lambda}, positive. The concentration parameter for the
#' \code{Laplace} prior.
#'
#' @details The SEPaLS estimators are built depending on the value given to
#' \code{type}:
#' \itemize{
#'     \item \code{vMF}: then the estimator is proportional to
#'     \deqn{ \hat{\beta}_{ml}(y_n) + \kappa_0\mu_0, }
#'     where \eqn{\hat{\beta}_{ml}(y_n)} is the EPLS estimator, which coincides
#'     with the maximum-likelihood estimator of SEPaLS for a threshold \eqn{y_n}.
#'     \item \code{Laplace}: then the estimator is proportional to
#'     \deqn{ S_\lambda\left(\hat{\beta}_{ml}(y_n)\right), }
#'     where \eqn{S_\lambda} is the soft-thresholding operator of threshold
#'     \eqn{\lambda}.
#' }
#'
#' @return A SEPaLS estimator
#' @export
#'
#' @seealso \code{\link{bootstrap.SEPaLS}}
#'
#' @examples
#' set.seed(1)
#' n <- 3000
#' p <- 10
#' X <- matrix(rnorm(n*p),n,p)
#' beta <- c(5:1,rep(0,p-5)) ; beta <- beta/sqrt(sum(beta^2))
#' Y <- (X%*%beta)^3 + rnorm(n,sd=1/3)
#' mu0 <- rnorm(p) ; mu0 <- mu0/sqrt(sum(mu0^2))
#' sepals_vMF <- SEPaLS(X,Y,yn=1,type="vMF",mu0=mu0,kappa0=1)
#' sepals_Laplace <- SEPaLS(X,Y,yn=1,type="Laplace",lambda=0.01)
SEPaLS <- function(X,Y,yn,type=c("vMF","Laplace"),
                   mu0=NULL,kappa0=NULL,
                   lambda=NULL){
  v_n <- maximum_Likelihood_SEPaLS(X,Y,yn)
  out <- v_n
  problem <- FALSE
  if(type=="vMF")
  {
    if(is.null(kappa0)){
      problem <- TRUE
      message("Please give a value to kappa0")
    }
    if(is.null(mu0)){
      problem <- TRUE
      message("Please give a value to mu0")
    }
    if(!problem){
      out <- out + kappa0*mu0
    }
  }
  else if(type=="Laplace")
  {
    if(is.null(lambda)){
      problem <- TRUE
      message("Please give a value to lambda")
    }
    if(!problem){
      out_abs <- abs(out) - lambda
      out_abs[which(out_abs<0)] <- 0
      out <- sign(out)*out_abs
    }
  }
  if(!problem){
    norm_out <- sqrt(sum(out^2))
    if(norm_out<1e-9)
    {
      out <- out*0
    }
    else
    {
      out <- out/sqrt(sum(out^2))
    }
  }
  return(out)
}

#' Bootstrap function for SEPaLS estimator.
#'
#' @param X \eqn{(n\times p)}-dimensional matrix of the covariates.
#' @param Y \eqn{(n)}-dimensional vector of the response.
#' @param yn \eqn{y_n} the quantile corresponding to lowest values of \eqn{Y}s
#' to put in the tail.
#' @param type character, whether \code{vMF} for von Mises-Fisher prior or
#' \code{Laplace} for Laplace prior. See details.
#' @param mu0  \eqn{\mu_0}, unitary \eqn{(p)}-dimensional vector. The direction
#'  parameter for the \code{vMF} prior.
#' @param kappa0 \eqn{\kappa_0}, positive. The concentration parameter for the
#' \code{vMF} prior.
#' @param lambda \eqn{\lambda}, positive. The concentration parameter for the
#' \code{Laplace} prior.
#' @param B positive integer. The number of bootstrap samples on which estimate
#' the SEPaLS directions. Default to 20.
#'
#' @return A list with two elements:
#' \itemize{
#'     \item \code{ws}: A \eqn{(B\times p)}-dimensional matrix  with each
#'     row corresponding to the \emph{SEPaLS} direction estimated on each
#'     bootstrap sample.
#'     \item \code{cor}: The correlation of each estimate direction on the
#'     Out-Of-Bag (OOB) sample with the response.
#' }
#'
#' @export
#'
#' @seealso \code{\link{SEPaLS}}
#'
#' @examples
#' set.seed(5)
#' n <- 3000
#' p <- 10
#' X <- matrix(rnorm(n*p),n,p)
#' beta <- c(5:1,rep(0,p-5)) ; beta <- beta/sqrt(sum(beta^2))
#' Y <- (X%*%beta)^3 + rnorm(n)
#' boot.sepals_Laplace <- bootstrap.SEPaLS(X,Y,yn=1,type="Laplace",lambda=0.01,
#' B=100)
#' boxplot(boot.sepals_Laplace$ws);abline(h=0,col="red",lty=2)
bootstrap.SEPaLS <- function(X,Y,yn,type=c("vMF","Laplace"),
                             mu0=NULL,kappa0=NULL,
                             lambda=NULL,B=20){
  p <- ncol(X)
  n <- nrow(X)
  ws <- matrix(NA,B,p)
  cor_oob <- rep(NA,B)
  for(i in 1:B)
  {
    id <- sample(1:n,size = n,replace = TRUE)
    X_train <- X[id,,drop=FALSE]
    Y_train <- Y[id]
    w <- SEPaLS(X=X_train,Y=Y_train,yn=yn,type = type,mu0=mu0,kappa0=kappa0,
                lambda=lambda)
    ws[i,] <- w
    id_oob <- (1:n)[-unique(id)]
    id_oob_thresh <- id_oob[which(Y[id_oob]>=yn)]
    if(length(id_oob_thresh)>2){
      X_oob <- X[id_oob_thresh,,drop=FALSE]
      Y_oob <- Y[id_oob_thresh]
      t_oob <- X_oob%*%w
      mu_t_oob <- sum(t_oob)/n
      mu_Y_oob <- sum(Y_oob)/n
      t_oob_scaled <- t_oob-mu_t_oob
      t_oob_scaled <- t_oob_scaled/sqrt(sum(t_oob_scaled^2)/n)
      Y_oob_scaled <- Y_oob-mu_Y_oob
      Y_oob_scaled <- Y_oob_scaled/sqrt(sum(Y_oob_scaled^2)/n)
      cor_oob[i] <- sum(t_oob_scaled*Y_oob_scaled)/n
    }
  }
  return(list(ws=ws,cor=cor_oob) )
}
