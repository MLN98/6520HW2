#' Elastic Net
#'
#' Fitting elastic net to data using coordinate descent algorithm
#'
#'@param y The response variable
#'@param X The covariate
#'@param beta The initialization of the coeffient. The default value is 0.
#'@param alpha Elastic net parameter, 0<=alpha<=1
#'@param lambda Regularization parameter
#'@param maxit Maximum Number of Iteration. The default value is 10000.
#'@param tol stop criterion. If the difference between the current result
#' and result from the last iteration is less than the Tolerance, then we say the algorithm converges
#' The default value is 1e-6
#'@return The solution that fits elastic net to the given data using
#'coordinate descent algorithm.
#'
#'@export
#'
elnet_coord <- function(y,X, beta, alpha,lambda, maxit= 10000, tol = 1e-6){
  X = as.matrix(X)
  n = dim(X)[1]
  m = dim(X)[2]

  Y = (y-mean(y))/sd(y)
  Xmean = apply(X,2,mean)
  Xsd = apply(X,2,sd)
  X = (X - matrix(Xmean,nrow(X),ncol(X),byrow=TRUE))%*%diag(1/Xsd)

  if(any(is.na(beta))){
    beta = rep(0, m)
  }
  tol.met=F
  beta1=beta
  iter=0

  while(!tol.met){

    beta0=beta1
    iter=iter+1

    for (j in 1:m){
      ej=Y-X[,-j]%*% beta1[-j]
      xj=X[,j]
      z=mean(xj*ej)
      num=sign(z)*max(0,abs(z)-lambda*alpha)
      denom=mean(xj^2)+lambda*(1-alpha)
      beta.j=num/denom
      beta1[j]=beta.j
    }

    if (any(abs(beta1) == Inf))
      stop("The algorithm diverges")

    if (iter == maxit) {
      stop("Max iteration reached")
    }

    if (norm(beta0-beta1,'2')<tol){
      tol.met=T
    }
  }

  return(beta1)
}
