#' Multitask Expected Improvement
#'
#' @param Xinit
#' @param f
#' @param N
#' @param g
#' @param lower
#' @param upper
#' @param add
#' @param optim.multi.start
#'
#' @return
#' @export
#'
#' @examples

multitask_EI <- function(Xinit, f, N = 50, g = NULL, lower = 0, upper=1, add = "both", optim.multi.start = NULL){
  
  X <- Xinit

  y <- f(Xinit)
  y.mod <- (y-mean(y))/sd(y)
  m <- ncol(X)
  if(is.null(optim.multi.start)) optim.multi.start <- min(2*m, 20)
  
  da <- laGP::darg(d = NULL, X = lhs::randomLHS(3000, m)*(upper-lower)+lower) # Setting a prior for theta
  if(is.null(g)){
    ga <- laGP::garg(list(mle = TRUE, max = var(y.mod)), y.mod)
    gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = ga$start, dK = TRUE)
    laGP::mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))$theta
  } else{
    gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = g, dK = TRUE)
    laGP::mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$d
  }

  while (nrow(X) < N) {
    predcur <- laGP::predGPsep(gpi, X, lite = TRUE)
    curmax <- max(predcur$mean)
    curmin <- min(predcur$mean)

    start <-  lhs::maximinLHS(optim.multi.start, m)*(upper-lower)+lower

    optim.res.max <- cbind(start, rep(0, optim.multi.start))

    for(i in 1:optim.multi.start){
      tryCatch({
        out <- optim(start[i,],
                     ESMultiTask::EI,
                     method = "L-BFGS-B",
                     lower = lower,
                     upper = upper,
                     model = gpi,
                     fmax=curmax,
                     control = list(fnscale=-1)
        )
        optim.res.max[i,] <- c(out$par, out$value)
      },
      error = function(e) e
      )
    }
    optim.res.max <- optim.res.max[which.max(optim.res.max[,m+1]),]

    optim.res.min <- cbind(start, rep(0, optim.multi.start))

    for(i in 1:optim.multi.start){
      tryCatch(
        {
          out <- optim(start[i,],
                       ESMultiTask::EI,
                       method = "L-BFGS-B",
                       lower = lower,
                       upper = upper,
                       model = gpi,
                       fmin=curmin,
                       maximisation = FALSE,
                       control = list(fnscale=-1)
          )
          optim.res.min[i,] <- c(out$par, out$value)
        },
        error = function(e) e
      )

    }
    optim.res.min <- optim.res.min[which.max(optim.res.min[,m+1]),]

    if(add == "both"){
      xnew <- rbind(optim.res.max[1:m], optim.res.min[1:m])
    } else{
      if(optim.res.max[m+1] > optim.res.min[m+1]){
        xnew <- matrix(optim.res.max[1:m], ncol=m)
      } else{
        xnew <- matrix(optim.res.min[1:m], ncol=m)
      }
    }

    ynew <- f(xnew)

    X <- rbind(X, xnew)
    y <- c(y, ynew)
    y.mod <- (y-mean(y))/sd(y)
    if(is.null(g)){
      ga <- laGP::garg(list(mle = TRUE, max = var(y.mod)), y.mod)
      gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = ga$start, dK = TRUE)
      laGP::mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))$theta
    } else{
      gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = g, dK = TRUE)
      laGP::mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$d
    }


  }

  res <- list(X=X, y = y)

  return(res)

}


