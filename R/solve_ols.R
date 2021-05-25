#' Linear system solver
#'
#' Iterative Methods for Solving Linear System of Equations
#'
#'@param A The A part for linear system of Ax = b, a matrix object
#'@param b The b part for linear system of Ax = b, a matrix object
#'@param x0 Initial value of the solution
#'@param method The method used for solver: "Gauss" or "Jacobi"
#'@param parallel TRUE for applying parallel computing for Jacobi method
#'@param num_cores Number of cores for parallel Jacobi.
#'@param maxit Maximum Number of Iteration
#'
#'@return The solution to the linear system
#'
#'@export
solve_ols <- function(A,b, x0,num_cores = 1, parallel = F, method = "Gauss", maxit = 10000){

  A = as.matrix(A)
  n = dim(A)[1]
  p = dim(A)[2]
  U= getU(A);L=getL(A); D=diag(diag(A));

  if (any(is.na(x0))){
    x0 = rep(0, p)
  }

  ###Solve the problem using designated method.
  if (method == "Gauss"){

    if (norm(solve(L+D,U),'2')>1){
      warning("Don't converge")
    }
    else{
      res =  f.iter(A, b, x0, GS,maxit)
    }
  }

  if (method == "Jacobi"){

    if(norm(-solve(D,L+U),'2')>1){
      warning("Don't converge")
    }

    else{
      if(!parallel){
        res =  f.iter(A, b, x0, Jac,maxit)
      }

      else{
        registerDoParallel(makeCluster(num_cores))
        res = f.iter(A, b, x0, Jac.P,maxit)
      }
    }
  }

  return (res)
}


#' obtain the upper triangular matrix of A
#'@param A Linear System of interest
#'@return the  upper triangular matrix of A
#'@author Yixiao Lin
#'@export
getU <- function(A){
  A[lower.tri(A)] = 0
  diag(A)=0
  return(A)
}

#' obtain the lower triangular matrix of A
#'@param A Linear System of interest
#'@return the lower triangular matrix of A
#'@author Yixiao Lin
#'@export
getL <- function(A){
  A[upper.tri(A)] = 0
  diag(A)=0
  return(A)
}


#' Gauss-Seidel Method
#'@param A Linear System of interest
#'@param b Target Vector
#'@param x Initial guess of the solution
#'@export
GS = function(A, b, x) {
  D = diag(A)
  n=length(x)
  diag(A) = 0
  for (i in 1:n) {
    x[i] = (b[i] - sum(A[i, ]* x))/D[i]
  }
  return(x)
}



#' Jacobi Method
#'@param A Linear System of interest
#'@param b Target Vector
#'@param x Initial guess of the solution
#'@export
Jac <- function(A, b, x){
  D = diag(A)
  diag(A) = 0
  n=length(x)
  x1 = rep(0,n)
  for (i in 1:n) {
    x1[i] =(b[i] - sum(A[i, ]*x))/D[i]
  }
  return(x1)
}


#' Parallel Jacobi
#'@param A Linear System of interest
#'@param b Target Vector
#'@param x Initial guess of the solution
#'@export
Jac.P <- function(A, b, x){
  D = diag(A)
  n=length(x)
  diag(A) = 0
  outlist = foreach(i = 1:n,.multicombine = TRUE)  %dopar% {
    (b[i]- sum(A[i, ]*x))/D[i]
  }
  x1 = as.vector(unlist(outlist))
  return(x1)
}


#' #' iterative method for GS, Jacobi and Parallel Jacobi
#' #'@param A Linear System of interest
#' #'@param b Target Vector
#' #'@param x0 Initial guess of the solution
#' #'@param f function used for solving Ax=b. There are 3 options: "Gauss" "Jacobi"  and "Parallel Jacobi"
#' #'@param maxit Maximum Number of Iteration
#' #'@export
#' f.iter <- function(A, b, x0, f, maxit){
#'   for (i in 1:maxit){
#'     x <- f(A, b, x0)
#'     if (norm(x0-x,'2')<1e-5){
#'       return(x)
#'     }
#'     if (any(abs(x) == Inf))
#'       stop("The algorithm diverges")
#'     x0 = x
#'
#'   }
#'   return(x)
#' }
