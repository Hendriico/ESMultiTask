#' TSSO algorithm
#'
#' Two-stage sequential optimisation (TSSO) algorithm uses modified expected improvement (MEI) to find
#' the minimum of a noisy blackbox function
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
#' @param rmin
#' @param predf
#' @param optim.multi.start
#'
#' @return
#' @export
#'
#' @examples
TSSO <- function(
  Xinit,
  f,
  N = 100,
  noisy = TRUE,
  sigma = 0,
  model.type = c("laGP", "km"),
  type = "UK",
  lower = 0,
  upper = 1,
  B = 10,
  rmin=2,
  predf = laGP::predGPsep,
  optim.multi.start = 5) {

  X <- Xinit

  I <- ceiling((N-nrow(Xinit))/B)

  y <- matrix(f(Xinit, noisy = noisy, sigma = sigma), ncol=1)

  m <- ncol(X)

  colnames(X) <- paste0("x", 1:m)
  colnames(y) <- "y"
  design <- cbind(data.frame(X), data.frame(y))
  design <- design %>%
    dplyr::group_by(.dots = colnames(X)) %>%
    dplyr::summarise(sample.mean = mean(y),
                     sample.sd = sd(y)) %>%
    dplyr::ungroup()


  if(model.type == "laGP"){

    da <- laGP::darg(list(mle = TRUE, max = 0.5), lhs::randomLHS(1000, 2))

      ga <- laGP::garg(list(mle = TRUE, max = var(y)), y)
      gpi.noisy <- laGP::newGPsep(X, y, d = da$start, g = ga$start, dK = TRUE)
      laGP::mleGPsep(gpi.noisy, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))$msg

      gpi.det <- laGP::newGPsep(design[,1:m], design$sample.mean, d = da$start, g = 1e-6, dK = TRUE)
      laGP::mleGPsep(gpi.det, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$msg

  } else if(model.type == "km"){
    gpi.noisy <- DiceKriging::km(design = X, response = y, covtype="gauss", control=list(trace=FALSE), nugget.estim = TRUE )
    gpi.det <- DiceKriging::km(design = design[,1:m], response = design$sample.mean, covtype="gauss", control=list(trace=FALSE) )

  } else {
    stop("model.type has to be either km or laGP")
  }



  ra <- numeric(I+1)

  for(i in 1:I) {

    Tn <- N-nrow(X)
    if(Tn < B) B <- Tn

    xmin <- design%>%
      dplyr::slice_min(sample.mean, n=1)

    if(model.type == "laGP"){
      fmin <- predf(gpi.noisy, xmin[1,1:m], lite = TRUE)$mean
    } else {
      fmin <- DiceKriging::predict.km(object = gpi.noisy, newdata = xmin[1,1:m], light.return = TRUE, cov.compute = FALSE, type = type)$mean
    }


    start <- data.matrix(xmin[1,1:m])
    if(optim.multi.start > 1){
      start <- rbind(start, lhs::maximinLHS(optim.multi.start-1, m))
    }

    optim.res <- matrix(nrow = optim.multi.start, ncol = m+1)

    for(j in 1:optim.multi.start){
      out <- optim(start[j,],
                   ESMultiTask::MEI,
                   method = "L-BFGS-B",
                   lower = lower,
                   upper = upper,
                   model.noisy = gpi.noisy,
                   model.det = gpi.det,
                   model.type = model.type,
                   predf = predf,
                   fmin = fmin,
                   control = list(fnscale=-1)
      )
      optim.res[j,] <- c(out$par, out$value)
    }

    xnew.s <- optim.res[which.max(optim.res[,m+1]),1:m]

    ra[i+1] <- ra[i] + ceiling(min((B-rmin)/I, N-ncol(Xinit)-i*B))

    if(ra[i+1] > B-rmin){
      ra[i+1] <- max(B-rmin, 0)
    }

    if(N-ncol(Xinit)-i*B-ra[i+1]>0){
      rs <- B-ra[i+1]
    }

    xnew.s <- matrix(rep(xnew.s, rs), ncol = m, byrow = TRUE)

    ynew.s <- matrix(f(xnew.s, noisy = noisy, sigma = sigma), ncol=1)


    xnew.a <- ESMultiTask::OCBA(design, ra = ra[i+1])$X
    ynew.a <- matrix(f(xnew.a, noisy = noisy, sigma = sigma), ncol=1)

    xnew <- rbind(xnew.s, xnew.a)
    ynew <- rbind(ynew.s, ynew.a)

    X <- rbind(X, xnew)
    y <- rbind(y, ynew.s, ynew.a)


    design <- cbind(data.frame(X), data.frame(y))
    design <- design %>%
      dplyr::group_by(.dots = colnames(X)) %>%
      dplyr::summarise(sample.mean = mean(y),
                       sample.sd = sd(y)) %>%
      dplyr::ungroup()

    if(model.type == "laGP"){
      laGP::updateGPsep(gpi.noisy, xnew, ynew)
        mle <- laGP::mleGPsep(gpi.noisy, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, 1), ab = c(da$ab, ga$ab))

        gpi.det <- laGP::newGPsep(design[,1:m], design$sample.mean, d = da$start, g = 1e-6, dK = TRUE)
        laGP::mleGPsep(gpi.det, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$msg

    } else{
      gpi.noisy <- DiceKriging::update(gpi.noisy, newX = xnew.s, newy = ynew.s, nugget.reestim = noisy, newX.alreadyExist = FALSE)
      gpi.noisy <- DiceKriging::update(gpi.noisy, newX = xnew.a, newy = ynew.a, nugget.reestim = noisy, newX.alreadyExist = TRUE)

    }



  }


  xmin <- design%>%
    dplyr::slice_min(sample.mean, n=1)




  return(list(xmin = xmin[1,1:m], X = X, y = y, ytrue = f(xmin[1,1:m], noisy = FALSE)))

}
