#' EQI optimisation
#'
#' This algorithm uses Expected Quantile Improvement (EQI) with Gaussian processes
#' to find the minimum of a blackbox function
#'
#' @param Xinit
#' @param f
#' @param N
#' @param noisy
#' @param sigma
#' @param model.type
#' @param type
#' @param lower
#' @param upper
#' @param B
#' @param beta
#' @param predf
#' @param optim.multi.start
#'
#' @return
#' @export
#'
#' @examples

EQIoptim <- function(
  Xinit,
  f,
  N = 100,
  noisy = TRUE,
  sigma = 0,
  model.type = c("laGP", "km"),
  type = "UK",
  lower = 0,
  upper = 1,
  B = 2,
  beta = 0.9,
  predf = laGP::predGPsep,
  optim.multi.start = 5) {

  X <- Xinit

  y <- f(Xinit, noisy = noisy, sigma = sigma)
  m <- ncol(X)

  if(model.type == "laGP"){
    da <- laGP::darg(list(mle = TRUE, max = 0.5), lhs::randomLHS(1000, 2))
    if (noisy) {
      ga <- laGP::garg(list(mle = TRUE, max = var(y)), y)
      gpi <- laGP::newGPsep(X, y, d = da$start, g = ga$start, dK = TRUE)
      laGP::mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))$msg
    } else {
      gpi <- laGP::newGPsep(X, y, d = da$start, g = 1e-6, dK = TRUE)
      laGP::mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$msg
    }
  } else if(model.type == "km"){
    gpi <- DiceKriging::km(design = X, response = y, covtype="gauss", control=list(trace=FALSE), nugget.estim = noisy )
    colnames(X) <- colnames(gpi@X)
  } else {
    stop("model.type has to be either km or laGP")
  }




  while (nrow(X) < N) {

    Tn <- N-nrow(X)
    if(Tn < B) B <- Tn

    if(model.type == "laGP"){
      cur.pred <- predf(gpi, X, lite = TRUE)
      qmin <- min(cur.pred$mean + qnorm(beta) * sqrt(cur.pred$s2))

    } else {
      cur.pred <- DiceKriging::predict.km(object = gpi, newdata = X+sqrt(.Machine$double.eps), light.return = TRUE, cov.compute = FALSE, type = type)
      qmin <- min(cur.pred$mean + qnorm(beta) * cur.pred$sd)
    }
    cur.eqi <- ESMultiTask::EQI(X, model = gpi, model.from = model.type, comp.model.noise.var = sigma^2/Tn, x.design = X)

    start <- matrix(X[which.max(cur.eqi),], nrow=1)
    if(optim.multi.start > 1){
      start <- rbind(start, lhs::maximinLHS(optim.multi.start-1, m)*(upper-lower) + lower)
    }

    optim.res <- matrix(nrow = optim.multi.start, ncol = m+1)

    for(i in 1:optim.multi.start){
      out <- optim(start[i,],
                   ESMultiTask::EQI,
                   method = "L-BFGS-B",
                   lower = lower,
                   upper = upper,
                   model = gpi,
                   model.from = model.type,
                   comp.model.noise.var = sigma^2/Tn,
                   beta = beta,
                   predf = predf,
                   x.design = X,
                   qmin = qmin,
                   control = list(fnscale=-1)
      )
      optim.res[i,] <- c(out$par, out$value)
    }

    xnew <- optim.res[which.max(optim.res[,m+1]),1:m]
    if(all(xnew == start[1,])){
      xexists <- TRUE
    } else{
      xexists <- FALSE
    }
    xnew <- matrix(rep(xnew, B), ncol = m, byrow = TRUE)

    ynew <- f(xnew, noisy = noisy, sigma = sigma)

    if(model.type == "laGP"){
      laGP::updateGPsep(gpi, xnew, ynew)
      if (noisy) {
        mle <- laGP::mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))
      } else {
        mle <- laGP::mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)
      }
    } else{
      gpi <- DiceKriging::update(gpi, newX = xnew, newy = ynew, nugget.reestim = noisy, newX.alreadyExist = xexists)
    }

    X <- rbind(X, xnew)
    y <- c(y, ynew)

  }


  if(model.type == "laGP"){
    predq <- predf(gpi, X, lite = TRUE)
    m <- which.min(predq$mean + qnorm(beta) * sqrt(predq$s2))
    xmin <- X[m,]
    laGP::deleteGPsep(gpi)
  } else{
    predq <- DiceKriging::predict.km(gpi, X+sqrt(.Machine$double.eps), type = type, light.return = TRUE)
    m <- which.min(predq$mean + qnorm(beta) * predq$sd)
    xmin <- X[m,]
  }


  return(list(xmin = xmin, X = X, y = y, ytrue = f(matrix(xmin, nrow=1), noisy = FALSE)))

}
