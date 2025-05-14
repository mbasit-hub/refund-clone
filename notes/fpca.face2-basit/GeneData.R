################################################################################
# n: the number of subjects.
# J: the number of observations in each curve.
# design: generate regular or irregular spaced data, default as "regular". 
# level: the level of sparse
# sigma: the standard deviation of random errors
# equally: equally-spaced, default as TRUE. 
################################################################################
GeneData <- function(n = 100, J = 100, sigma = 1, equally = TRUE,
                     design = "regular", level = 0.5){
  
  if(equally){
    t <- seq(0, 1, length = J)
  } else {
    tJ = 5*J
    t0 <- seq(0, 1, length = tJ)
    idx = c(1, sort(sample(2:(tJ-1), J-2, replace=F)), tJ)
    t <- t0[idx]
  }
  
  
  # K <- 4
  # lambda <- c(2.5, 1.5, 1, 0.5)
  # ### True Eigenfunctions
  # case = 1
  # if(case==1) phi <- sqrt(2)*cbind(sin(2*pi*t),cos(2*pi*t),
  #                                  sin(4*pi*t),cos(4*pi*t))
  # if(case==2) phi <- cbind(rep(1,J),sqrt(3)*(2*t-1),
  #                          sqrt(5)*(6*t^2-6*t+1),
  #                          sqrt(7)*(20*t^3-30*t^2+12*t-1))
  # 
  
  K <- 3
  lambda <- c(1, 0.5, 0.25)
  ### True Eigenfunctions
  case = 2
  if(case==1) phi <- sqrt(2)*cbind(sin(2*pi*t),sin(4*pi*t),cos(4*pi*t))
  if(case==2) phi <- cbind(rep(1,J),sqrt(3)*(2*t-1),sqrt(5)*(6*t^2-6*t+1))
  
  # Generate scores
  xi <- matrix(0, nrow=n, ncol=K)
  for(k in 1:K) {
    xi[,k] <- rnorm(n, sd=sqrt(lambda[k]))
  }
  # Generate errors
  epsilon <- matrix(rnorm(n*J,sd=sigma), nc=J)
  
  # Generate dense data
  X <- xi %*% t(phi) 
  Y0 <- X +  epsilon 
  
  
  # Generate sparse data
  if (design == "regular") {
    Y <- Y0
  } else {
    nobs <- floor(level*J)
    Y <- matrix(NA,nrow=n,ncol=J)
    for (i in 1:n) {
      idx <- sample(1:J, nobs)
      Y[i,idx] <- Y0[i,idx]
    }
  }
  
  return(list(Y = Y, evalues = lambda, efunctions = phi, scores = xi, argvals = t))
}


  
  
  
  
  
  
  
  
  
  
                  
                  
                  
                  
                  
                  
                  
 
                  
                  
                  
                  
                  
                  
                  
                  
                                   
                  
                  