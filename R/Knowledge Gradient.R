#' Knowledge Gradient
#'
#'This algorithm uses either Correlated KG or Approximate KG criteria with Gaussian processes
#'to find the minimum of a blackbox function
#'
#' @param Xinit
#' @param f
#' @param N
#' @param noisy
#' @param sigma
#' @param model.type
#' @param KG.type
#' @param type
#' @param lower
#' @param upper
#' @param B
#' @param optim.multi.start
#'
#' @return
#' @export
#'
#' @examples


KG <- function(
  Xinit,
  f,
  N = 100,
  noisy = TRUE,
  sigma = 0,
  model.type = c("laGP", "km"),
  KG.type = "Correlated",
  type = "UK",
  lower = 0,
  upper = 1,
  B = 2,
  optim.multi.start = 5) {

  if(KG.type == "Correlated"){
    KGfun <- ESMultiTask::CKG
  } else if(KG.type == "Approximate"){
    KGfun <- ESMultiTask::AKG
  } else{
    stop("KG.type has to be either Correlated or Approximate")
  }

  X <- Xinit

  y <- f(Xinit, noisy = noisy, sigma = sigma)
  m <- ncol(X)

  if(model.type == "laGP"){
    da <- darg(list(mle = TRUE, max = 0.5), randomLHS(1000, 2))
    if (noisy) {
      ga <- garg(list(mle = TRUE, max = var(y)), y)
      gpi <- newGPsep(X, y, d = da$start, g = ga$start, dK = TRUE)
      mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))$msg
    } else {
      gpi <- newGPsep(X, y, d = da$start, g = 1e-6, dK = TRUE)
      mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$msg
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

    cur.pred <- apply(X, 1, KGfun, model = gpi, new.noise.var = sigma^2, type = type)
    cur.min <- X[which.max(cur.pred),]



    start <- matrix(cur.min, nrow=1)
    if(optim.multi.start > 1){
      start <- rbind(start, lhs::maximinLHS(optim.multi.start-1, m)*(upper-lower)+lower)
    }

    optim.res <- matrix(nrow = optim.multi.start, ncol = m+1)

    for(i in 1:optim.multi.start){
      out <- optim(start[i,],
                   KGfun,
                   method = "L-BFGS-B",
                   lower = lower,
                   upper = upper,
                   model = gpi,
                   new.noise.var = sigma^2,
                   type = type,
                   control = list(fnscale=-1)
      )
      optim.res[i,] <- c(out$par, out$value)
    }

    xnew <- optim.res[which.max(optim.res[,m+1]),1:m]

    xnew <- matrix(rep(xnew, B), ncol = m, byrow = TRUE)

    ynew <- f(xnew, noisy = noisy, sigma = sigma)

    if(model.type == "laGP"){
      updateGPsep(gpi, xnew, ynew)
      if (noisy) {
        mle <- mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, 1), ab = c(da$ab, ga$ab))
      } else {
        mle <- mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)
      }
    } else{
      gpi <- DiceKriging::update(gpi, newX = xnew, newy = ynew, nugget.reestim = noisy)
    }

    X <- rbind(X, xnew)
    y <- c(y, ynew)

  }

  if(model.type == "laGP"){
    cur.pred <- predf(gpi, X, lite = TRUE)$mean
    xmin <- X[which.min(cur.pred),]
    deleteGPsep(gpi)
  } else {
    cur.pred <- DiceKriging::predict.km(object = gpi, newdata = X, light.return = TRUE, cov.compute = FALSE, type = type)$mean
    xmin <- X[which.min(cur.pred),]
  }



  return(list(xmin = xmin, X = X, y = y, ytrue = f(xmin, noisy = FALSE)))

}


