

#' MQ algorithm
#'
#' MQ algorithm uses kriging quantiles with Gaussian processes to find the minimum of a blackbox function
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

MQ <- function(
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
  beta = 0.1,
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

      cur.pred <- ESMultiTasktim::KQ(X, gpi, model.from = model.type, beta = beta, type = type, predf = predf)
      cur.min <- X[which.min(cur.pred),]

    start <- matrix(cur.min, nrow=1)
    if(optim.multi.start > 1){
      start <- rbind(start, lhs::maximinLHS(optim.multi.start-1, m))
    }

    optim.res <- matrix(nrow = optim.multi.start, ncol = m+1)

    for(i in 1:optim.multi.start){
      out <- optim(start[i,],
                   ESMultiTasktim::KQ,
                   method = "L-BFGS-B",
                   lower = lower,
                   upper = upper,
                   model = gpi,
                   model.from = model.type,
                   beta = beta,
                   predf = predf,
                   type = type
                   )
      optim.res[i,] <- c(out$par, out$value)
    }

    xnew <- optim.res[which.min(optim.res[,m+1]),1:m]

    xnew <- matrix(rep(xnew, B), ncol = m, byrow = TRUE)

    ynew <- f(xnew, noisy = noisy, sigma = sigma)

    if(model.type == "laGP"){
      laGP::updateGPsep(gpi, xnew, ynew)
    if (noisy) {
      mle <- laGP::mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, 1), ab = c(da$ab, ga$ab))
    } else {
      mle <- laGP::mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)
    }
    } else{
      gpi <- DiceKriging::update(gpi, newX = xnew, newy = ynew, nugget.reestim = noisy)
    }

    X <- rbind(X, xnew)
    y <- c(y, ynew)

  }

  KQmin <- ESMultiTasktim::KQ(X, gpi, model.from = model.type, beta= beta, predf=predf, type=type)

  if(model.type == "laGP") laGP::deleteGPsep(gpi)
  xmin <- X[which.min(KQmin),]

    return(list(xmin = xmin, X = X, y = y, ytrue = f(xmin, noisy = FALSE)))

}

