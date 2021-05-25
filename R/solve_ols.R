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
solve_ols <- function(A,b, x0, num_cores = 1, parallel = F, method = "Gauss", maxit = 10000){
  A = as.matrix(A)
  n = dim(A)[1]
  m = dim(A)[2]
  U = A
  L = A
  L[upper.tri(L, diag = TRUE)] = 0
  U[lower.tri(U, diag = TRUE)] = 0
  D=diag(diag(A))

  if (method == "Gauss"){
    if (norm(solve(L+D,U),'2')>1){
      warning("The solution doesn't converge with Gauss-Seidel method.")
    }
    else{
      for (i in 1:maxit){
        # Apply Gauss-seidel
        x_prev = x0
        A_ = A
        D = diag(A_)
        diag(A_) = 0
        for (i in 1:length(x0)) {
          x0[i] = (b[i] - sum(A_[i,]* x0))/D[i]
        }
        # Return value when converge
        if (norm(x0-x_prev,'2')<1e-5){
          return(x0)
        }
      }
    }
  }

  if (method == "Jacobi"){
    if(norm(-solve(D,L+U),'2')>1){
      warning("The solution doesn't converge with Jacobi method.")
    }
    else{
      if(!parallel){
        for (i in 1:maxit){
          A_ = A
          D = diag(A_)
          diag(A_) = 0
          n=length(x0)
          x = rep(0,n)
          for (i in 1:n) {
            x[i] =(b[i] - sum(A_[i, ]*x0))/D[i]
          }
          if (norm(x0-x,'2')<1e-5){
            return(x)
          }
          x0 = x
        }
        res = x
      }

      if(parallel == TRUE){
        registerDoParallel(makeCluster(num_cores))
        for (i in 1:maxit){
          A_ = A
          D = diag(A_)
          diag(A_) = 0
          outlist = foreach(i = 1:length(x0),.multicombine = TRUE)  %dopar% {
            (b[i]- sum(A_[i, ]*x0))/D[i]
          }
          x = as.vector(unlist(outlist))
          if (norm(x0-x,'2')<1e-5){
            return(x)
          }
          x0 = x
        }
        res = x
      }
    }
  }
  return(res)
}
