#' Leveraging Implementation
#'
#' Implements of algorithmic leveraging for linear regression using uniform and leverage score based subsampling of rows
#'
#'@param y The response in vector form
#'@param X The covariate in matrix form
#'@param sample_size sample size for leveraging
#'@param method The method to implement, 'weighted' or 'uniform'
#'@return  coefficients of the linear regression model
#'@export
algo_leverage <- function(y,X,sample_size,method ="uniform"){
  X = as.matrix(X)
  n = dim(X)[1]
  if(method == 'uniform'){
    pi = rep(1/n, n)
    n=dim(X)[1]
    m=dim(X)[2]
    index = seq(1,n,1)
    idx = sample(index, size = sample_size, prob = pi,replace = TRUE)
    sampleX=X[idx,]
    sampleY=y[idx]
    pi_=1/ (sample_size*pi[idx]^0.5)
    model <- lm(sampleY ~ sampleX, weights=pi_)
    beta = model$coefficients
    return(beta)
  }

  if(method == 'weighted'){
    H =X%*%solve(t(X)%*%X)%*%t(X)
    pi = diag(H)/sum(diag(H))
    n=dim(X)[1]
    p=dim(X)[2]
    idx=sample(n, size = sample_size, replace = TRUE, prob=pi)
    sampleX=X[idx,]
    sampleY=y[idx]
    pi_=1/ (sample_size*pi[idx]^0.5)
    model <- lm(sampleY ~ sampleX, weights=pi_)
    beta = model$coefficients
    return(beta)
  }
}
