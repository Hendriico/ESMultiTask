#' Kriging Quantile
#'
#' @param x
#' @param model
#' @param model.from
#' @param comp.model.noise.var
#' @param beta
#' @param predf
#'
#' @return
#' @export
#'
#' @examples
KQ <- function(x,
                 model,
                 model.from = c("laGP", "km"),
                 beta = 0.1,
                 predf = laGP::predGPsep,
                 type = "UK"){
  if(beta<=0 | beta > 0.5) stop("beta has to be between (0, 0.5]")
  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)

  if(model.from == "laGP"){
    pred <- predf(model,  x, lite=TRUE, nonug = TRUE)
    kq <- pred$mean + qnorm(beta) * sqrt(pred$s2)

  } else if (model.from == "km"){
    colnames(x) <- colnames(model@X)
    pred <- DiceKriging::predict.km(model,  x, type = type, light.return = TRUE)
    kq <- pred$mean + qnorm(beta) * pred$sd
  } else{
    stop("model.from has to be one of laGP or km")
  }
  return(kq)
}


#' Expected improvement
#'
#' @param x
#' @param model
#' @param fmax
#' @param fmin
#' @param maximisation
#' @param model.type
#'
#' @return
#' @export
#'
#' @examples

EI <- function(x,
               model,
               fmax,
               fmin,
               maximisation = TRUE,
               model.type = "laGP"
               ){
  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)

  if(model.type == "laGP"){
    pred <- laGP::predGPsep(model,  x, lite=TRUE, nonug = TRUE)
    predmean <- pred$mean
    predsd <- sqrt(pred$s2)
  } else if(model.type == "robust"){
    pred <- RobustGaSP::predict.rgasp(model, x)
    predmean <- pred$mean
    predsd <- pred$sd
  } else if(model.type == "km"){
      pred <- DiceKriging::predict.km(model, x+sqrt(.Machine$double.eps), type = "SK", light.return = TRUE)
      predmean <- pred$mean
      predsd <- pred$sd
  } else{
      stop("model.type has to be either laGP, km or robust")
    }


  if(maximisation){
    meas <- (predmean-fmax)/predsd
  res <- (predmean - fmax)*pnorm(meas) + predsd *dnorm(meas)
  } else{
    meas <- (fmin-predmean)/predsd
    res <- (fmin-predmean)*pnorm(meas) + predsd *dnorm(meas)
  }

  return(res)
  }


#' Augmented Expected Improvement
#'
#' @param x
#' @param model
#' @param model.from
#' @param comp.model.noise.var
#' @param beta
#' @param fmin
#' @param x.design
#' @param predf
#'
#' @return
#' @export
#'
#' @examples
#'

AEI <- function(x,
                  model,
                  model.from = c("laGP", "km"),
                  comp.model.noise.var = 0,
                  beta = 0.9,
                  fmin = NULL,
                  x.design = NULL ,
                  predf = laGP::predGPsep){

  if(beta<0.5 | beta >= 1) stop("beta has to be between [0.5, 1)")
  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)

  if(model.from == "laGP"){

    if(is.null(fmin)){
      if(is.null(x.design)){
        return("x.design must be submitted when qmin is not given")
      } else{
        predq <- predf(model, x.design, lite = TRUE, nonug = TRUE)
        m <- which.min(predq$mean + qnorm(beta) * sqrt(predq$s2))
        fmin <- predq$mean[m]
      }}

    pred <- predf(model,  x, lite=TRUE, nonug = TRUE)
    d <- fmin - pred$mean
    sigma <- sqrt(pred$s2)

  } else if (model.from == "km"){
    colnames(x) <- colnames(model@X)

    if(is.null(fmin)){
      predq <- DiceKriging::predict.km(model, model@X+sqrt(.Machine$double.eps), type = "UK", light.return = TRUE)
      m <- which.min(predq$mean + qnorm(beta) * predq$sd)
      fmin <- predq$mean[m]
    }

    pred <- DiceKriging::predict.km(model,  x, type = "UK", light.return = TRUE)
    d <- fmin - pred$mean
    sigma <- pred$sd
  } else{
    stop("model.from has to be one of laGP or km")
  }

  dn <- d / sigma
  ei <- d * pnorm(dn) + sigma * dnorm(dn)

  maei <- ei * (1-sqrt(comp.model.noise.var/(sigma^2+comp.model.noise.var)))

  return(maei)
}






#' Expected Quantile Improvement
#'
#' @param x
#' @param model
#' @param model.from
#' @param comp.model.noise.var
#' @param beta
#' @param qmin
#' @param x.design
#' @param predf
#'
#' @return
#' @export
#'
#' @examples


EQI <- function(x,
                    model,
                    model.from = c("laGP", "km"),
                    comp.model.noise.var = 0,
                    beta = 0.9,
                    qmin = NULL,
                    x.design = NULL ,
                    predf = laGP::predGPsep){

    if(beta<0.5 | beta >= 1) stop("beta has to be between [0.5, 1)")
    if(is.null(nrow(x))) x <- matrix(x, nrow = 1)

    if(model.from == "laGP"){
      pred <- predf(model,  x, lite=TRUE, nonug = TRUE)

      if(is.null(qmin)){
        if(is.null(x.design)){
          return("x.design must be submitted when qmin is not given")
        } else{
          predq <- predf(model, x.design, lite = TRUE, nonug = TRUE)
          qmin <- min(predq$mean + qnorm(beta) * sqrt(predq$s2))
        }}

      mn <- pred$mean
      sn2 <- pred$s2
    } else if (model.from == "km"){
      colnames(x) <- colnames(model@X)
      pred <- DiceKriging::predict.km(model,  x, type = "UK", light.return = TRUE)

      if(is.null(qmin)){
        predq <- DiceKriging::predict.km(model, model@X+sqrt(.Machine$double.eps), type = "UK", light.return = TRUE)
        qmin <- min(predq$mean + qnorm(beta) * predq$sd)
      }

      mn <- pred$mean
      sn2 <- pred$sd^2
    } else{
      stop("model.from has to be one of laGP or km")
    }



    mq <- mn + qnorm(beta) * sqrt((comp.model.noise.var * sn2)/(comp.model.noise.var + sn2))

    sq <- sn2/sqrt(comp.model.noise.var + sn2)

    eqi <- (qmin - mq)*pnorm((qmin-mq)/sq) + sq*dnorm((qmin-mq)/sq)

    return(eqi)
  }


#' Modified Expected Improvement
#'
#' @param x matrix of data points at which MEI should be calculated at
#' @param model.noisy model of class km or laGP
#' @param model.det model of class km or laGP
#' @param model.type km or laGP
#' @param type kriging type, e.g UK or SK
#' @param x.design current design matrix
#' @param y.design observed values from the computer model
#' @param fmin kriging prediction at the point with the lowest sample mean
#' @param predf function for making predictions from the Gaussian process
#'
#' @return
#' @export
#'
#' @examples
MEI <- function(x,
                  model.noisy,
                  model.det = NULL,
                  model.type = "laGP",
                  type = "UK",
                  x.design = NULL ,
                  y.design = NULL,
                  fmin = NULL,
                  predf = laGP::predGPsep){

  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)
  if(is.null(nrow(y.design))) y.design <- matrix(x, ncol = 1)

  if(is.null(fmin)){
    if(is.null(x.design) | is.null(y.design)){
      return("x.design and y.design must be submitted when fmin is not given")
    } else{
      d <- ncol(x.design)
      colnames(x.design) <- paste0("x", 1:d)
      colnames(y.design) <- "y"
      design <- cbind(data.frame(x.design), data.frame(y.design))
      design <- design %>%
        dplyr::group_by(.dots = colnames(x.design)) %>%
        dplyr::summarise(sample.mean = mean(y)) %>%
        dplyr::ungroup() %>%
        dplyr::slice_min(sample.mean, n=1)
      if(model.type == "laGP"){
        fmin <- predf(model.noisy, matrix(design[1,1:d], nrow = 1), lite = TRUE, nonug = TRUE)$mean
      } else {
        fmin <- DiceKriging::predict.km(model.noisy, matrix(design[1,1:d], nrow = 1), type =  type)$mean
      }

    }}

  if(is.null(model.det)){
    if(is.null(x.design) | is.null(y.design)){
      return("x.design and y.design must be submitted when model.det is not given")
    } else{
      d <- ncol(x.design)
      colnames(x.design) <- paste0("x", 1:d)
      colnames(y.design) <- "y"
      design <- cbind(data.frame(x.design), data.frame(y.design))
      design <- design %>%
        dplyr::group_by(.dots = colnames(x.design)) %>%
        dplyr::summarise(sample.mean = mean(y)) %>%
        dplyr::ungroup()
      if(model.type == "laGP"){
        da <- laGP::darg(list(mle = TRUE, max = 0.5), lhs::randomLHS(1000, 2))
        model.det <- laGP::newGPsep(design[,1:d], design$sample.mean, d = da$start, g = 1e-6, dK = TRUE)
        laGP::mleGPsep(model.det, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$msg
      } else {
        model.det <- DiceKriging::km(design = design[,1:d], response = design$sample.mean, covtype="matern5_2", control=list(trace=FALSE) )
      }

    }
  }

  if(model.type == "laGP"){
    pred.noisy <- predf(model.noisy,  x, lite=TRUE, nonug = TRUE)
    pred.det <- predf(model.det, x, lite=TRUE, nonug=TRUE)
    sigma <- sqrt(pred.det$s2)
  } else {
    colnames(x) <- colnames(model.noisy@X)
    pred.noisy <- DiceKriging::predict.km(model.noisy, x, type =  type)
    colnames(x) <- colnames(model.det@X)
    pred.det <- DiceKriging::predict.km(model.det, x+sqrt(.Machine$double.eps), type =  type)
    sigma <- pred.det$sd
  }

  d <- fmin - pred.noisy$mean




  dn <- d / sigma
  mei <- d * pnorm(dn) + sigma * dnorm(dn)


  return(mei)
}



#' Correlated Knowledge Gradient
#'
#' @param x
#' @param model
#' @param new.noise.var
#' @param type
#'
#' @return
#' @export
#'
#' @examples
CKG <- function(x,
                  model,
                  new.noise.var = 0,
                  type = "SK"){

  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)


  colnames(x) <- colnames(model@X)
  pred <- DiceKriging::predict.km(model, x, type = "SK")
  mn.x <- pred$mean
  sn.x <- pred$sd
  k.x <- pred$c
  Tinvk.x <- pred$Tinv.c

  pred.cur <- DiceKriging::predict.km(model, model@X+sqrt(.Machine$double.eps), type = "SK")

  mn.cur <- pred.cur$mean
  Tinvk.cur <- pred.cur$Tinv.c

  A <- -c( mn.cur, mn.x)


  B <- k.x - t(Tinvk.cur)%*%Tinvk.x
  B <- c(B, sn.x^2)
  B <- B/sqrt(sn.x^2 + new.noise.var)

  Isort <- order(x = B, y = A)
  b <- B[Isort]
  a <- A[Isort]
  Iremove <- numeric()
  for (i in 1:(model@n)) {
    if (b[i + 1] == b[i]) {
      Iremove <- c(Iremove, i)
    }
  }
  if (length(Iremove) > 0) {
    b <- b[-Iremove]
    a <- a[-Iremove]
  }
  C <- c(-.Machine$double.xmax, .Machine$double.xmax)
  M <- length(a)
  A <- 1

  for(i in 2:M){
    C[i+1] <- .Machine$double.xmax
    loopdone <- FALSE
    while (loopdone == FALSE) {
      j <- A[length(A)]
      C[j+1] <- (a[j]-a[i])/(b[i]-b[j])
      if(length(A) != 1 && C[j+1] <= C[A[length(A)-1]+1]){
        A <- A[-length(A)]
      } else {
        A <- c(A, i)
        loopdone <- TRUE
      }

    }
  }

  agz <- a[A]
  bgz <- b[A]
  cgz <- c(C[A+1], .Machine$double.xmax)
  M <-  length(A)
  ckg <- 0
  for(i in 1:(M-1)){
    ckg <- ckg + (bgz[i+1]-bgz[i]) *(dnorm(-abs(cgz[i])) -abs(cgz[i])*pnorm(-abs(cgz[i])))
  }

  return(ckg)

}




#' Approximate Knowledge Gradient
#'
#' @param x
#' @param model
#' @param new.noise.var
#' @param type
#'
#' @return
#' @export
#'
#' @examples

AKG <- function (x,
          model,
          new.noise.var = 0,
          type = "SK"
          )
{
  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)

  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  tau2.new <- new.noise.var
  predx <- DiceKriging::predict.km(model, newdata = newdata, type = type,
                      checkNames = FALSE)
  mk.x <- predx$mean
  sk.x <- predx$sd
  c.x <- predx$c
  V.x <- predx$Tinv.c


    predX <- DiceKriging::predict.km(model, newdata = model@X+sqrt(.Machine$double.eps), type = type,
                        checkNames = FALSE)
    mk.X <- predX$mean
    V.X <- predX$Tinv.c
    m_min <- min(c(mk.X, mk.x))

    cn <- c.x - t(V.X) %*% V.x
    cn <- c(cn, sk.x^2)
    A <- c(mk.X, mk.x)
    B <- cn/sqrt(tau2.new + sk.x^2)
    sQ <- B[length(B)]
    A <- -A
    nobs <- model@n
    Isort <- order(x = B, y = A)
    b <- B[Isort]
    a <- A[Isort]
    Iremove <- numeric()
    for (i in 1:(nobs)) {
      if (b[i + 1] == b[i]) {
        Iremove <- c(Iremove, i)
      }
    }
    if (length(Iremove) > 0) {
      b <- b[-Iremove]
      a <- a[-Iremove]
    }
    nobs <- length(a) - 1
    C <- rep(0, nobs + 2)
    C[1] <- -1e+36
    C[length(C)] <- 1e+36
    A1 <- 0
    for (k in 2:(nobs + 1)) {
      nondom <- 1
      if (k == nobs + 1) {
        nondom <- 1
      }
      else if ((a[k + 1] >= a[k]) && (b[k] == b[k + 1])) {
        nondom <- 0
      }
      if (nondom == 1) {
        loopdone <- 0
        count <- 0
        while (loopdone == 0 && count < 1000) {
          count <- count + 1
          u <- A1[length(A1)] + 1
          C[u + 1] <- (a[u] - a[k])/(b[k] - b[u])
          if ((length(A1) > 1) && (C[u + 1] <= C[A1[length(A1) -
                                                    1] + 2])) {
            A1 <- A1[-length(A1)]
          }
          else {
            A1 <- c(A1, k - 1)
            loopdone <- 1
          }
        }
      }
    }
    at <- a[A1 + 1]
    bt <- b[A1 + 1]
    ct <- C[c(1, A1 + 2)]
    maxNew <- 0
    for (k in 1:length(at)) {
      maxNew <- maxNew + at[k] * (pnorm(ct[k + 1]) - pnorm(ct[k])) +
        bt[k] * (dnorm(ct[k]) - dnorm(ct[k + 1]))
    }
    AKG <- maxNew - (-m_min)

  return(AKG)
}
