#' Michalewicz function
#'
#' @param X
#' @param noisy
#' @param sigma
#' @param m
#'
#' @return
#' @export
#'
#' @examples

michalewicz <- function(X, noisy = TRUE, sigma = 0.1, m=10){
  X <- X*pi
  if(is.null(nrow(X))) X <- matrix(X, nrow=1, byrow = TRUE )

  res <- numeric(nrow(X))
  for(i in 1:nrow(X)){
    xx <- X[i,]

    ii <- c(1:length(xx))
    sum <- sum(sin(xx) * (sin(ii*xx^2/pi))^(2*m))

    res[i] <- -sum
  }
  if(noisy) res <- res + rnorm(nrow(X), sd = sigma)
  return(res)
}

#' Shekel function
#'
#' @param X
#' @param noisy
#' @param sigma
#' @param m
#'
#' @return
#' @export
#'
#' @examples

shekel <- function(X, noisy = TRUE, sigma = 0.1, m=5){
  if(is.null(nrow(X))) X <- matrix(X, ncol = 4, byrow = TRUE)
  X <- X*10

  res <- numeric(nrow(X))
  b <- 0.1 * c(1, 2, 2, 4, 4, 6, 3, 7, 5, 5)[1:m]
  C <- c(4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0,
         4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6,
         4.0, 1.0, 8.0, 6.0, 3.0, 2.0, 5.0, 8.0, 6.0, 7.0,
         4.0, 1.0, 8.0, 6.0, 7.0, 9.0, 3.0, 1.0, 2.0, 3.6)
  C <- matrix(C, 4, 10, byrow=TRUE)
  Ct <- t(C)
  for(i in 1:nrow(X)){
    xx <- X[i,]
    xxmat <- matrix(rep(xx,times=m), m, 4, byrow=TRUE)
    inner <- rowSums((xxmat-Ct[1:m,1:4])^2)

    outer <- sum(1/(inner+b))

    res[i] <- -outer
  }

  if(noisy) res <- res + rnorm(nrow(X), sd = sigma)
  return(res)
}

#' Hartrmann 3 function
#'
#' @param X
#' @param noisy
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples

hart3 <- function(X, noisy = TRUE, sigma = 0.1){
  if(is.null(nrow(X))) X <- matrix(X, ncol = 3, byrow = TRUE)

  res <- numeric(nrow(X))
  alpha <- c(1.0, 1.2, 3.0, 3.2)

  A <- c(3.0, 10, 30,
         0.1, 10, 35,
         3.0, 10, 30,
         0.1, 10, 35)

  A <- matrix(A, 4, 3, byrow=TRUE)

  P <- 10^(-4) * c(3689, 1170, 2673,
                   4699, 4387, 7470,
                   1091, 8732, 5547,
                   381, 5743, 8828)

  P <- matrix(P, 4, 3, byrow=TRUE)
  for(i in 1:nrow(X)){
    xx <- X[i,]
    xxmat <- matrix(rep(xx,times=4), 4, 3, byrow=TRUE)
    inner <- rowSums(A[,1:3]*(xxmat-P[,1:3])^2)
    outer <- sum(alpha * exp(-inner))

    res[i] <- -outer
  }

  if(noisy) res <- res + rnorm(nrow(X), sd = sigma)
  return(res)
}

#' Hartmann 6 function
#'
#' @param X
#' @param noisy
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples

hart6 <- function(X, noisy = TRUE, sigma = 0.1){
  if(is.null(nrow(X))) X <- matrix(X, ncol = 6, byrow = TRUE)

  res <- numeric(nrow(X))
  alpha <- c(1.0, 1.2, 3.0, 3.2)
  A <- c(10, 3, 17, 3.5, 1.7, 8,
         0.05, 10, 17, 0.1, 8, 14,
         3, 3.5, 1.7, 10, 17, 8,
         17, 8, 0.05, 10, 0.1, 14)
  A <- matrix(A, 4, 6, byrow=TRUE)
  P <- 10^(-4) * c(1312, 1696, 5569, 124, 8283, 5886,
                   2329, 4135, 8307, 3736, 1004, 9991,
                   2348, 1451, 3522, 2883, 3047, 6650,
                   4047, 8828, 8732, 5743, 1091, 381)
  P <- matrix(P, 4, 6, byrow=TRUE)

  for(i in 1:nrow(X)){
    xx <- X[i,]
    xxmat <- matrix(rep(xx,times=4), 4, 6, byrow=TRUE)
    inner <- rowSums(A[,1:6]*(xxmat-P[,1:6])^2)
    outer <- sum(alpha * exp(-inner))

    res[i] <- -(2.58 + outer) / 1.94
  }

  if(noisy) res <- res + rnorm(nrow(X), sd = sigma)
  return(res)
}
