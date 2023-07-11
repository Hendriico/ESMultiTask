#' eTSSO algorithm
#'
#' Extended Two-stage sequential optimisation (eTSSO) algorithm uses modified expected improvement (MEI),
#' as well as varying number of replicated samples to find
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
#' @param rmin
#' @param predf
#' @param optim.multi.start
#'
#' @return
#' @export
#'
#' @examples
eTSSO <- function(
  Xinit,
  f,
  N = 100,
  noisy = TRUE,
  sigma = 0,
  model.type = c("laGP", "km"),
  type = "UK",
  lower = 0,
  upper = 1,
  rmin=2,
  predf = laGP::predGPsep,
  optim.multi.start = 5) {



  X <- Xinit

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

  N0 <- nrow(Xinit)/nrow(design)

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


  B <- numeric()
  B[1] <- rmin

  i <- 2
  while(nrow(X) < N) {


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

    xnew.s <- matrix(optim.res[which.max(optim.res[,m+1]),1:m], nrow=1)

    var.ocba <- ESMultiTask::OCBA(design, ra = 1)$sd^2
    if(model.type == "laGP"){
      sd2.mei <- predf(gpi.det, xnew.s, lite = TRUE, nonug = TRUE)$s2
    } else{
      colnames(xnew.s) <- colnames(gpi.det@X)
      sd2.mei <- predict.km(gpi.det, xnew.s, light.return = TRUE,type = type)$sd^2
    }

    B[i] <- max(floor(B[i-1]*(1+var.ocba/(var.ocba+sd2.mei))), N0+i-1)
    Tn <- N - nrow(X)
    if(Tn < B[i]) B[i] <- Tn
    if(Tn < rmin) rmin <- Tn

    xnew.s <- matrix(rep(xnew.s, rmin), ncol = m, byrow = TRUE)

    ynew.s <- matrix(f(xnew.s, noisy = noisy, sigma = sigma), ncol=1)


    xnew.a <- ESMultiTask::OCBA(design, ra = B[i]-rmin)$X
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
    i <- i+1



  }


  xmin <- design%>%
    dplyr::slice_min(sample.mean, n=1)




  return(list(xmin = xmin[1,1:m], X = X, y = y, ytrue = f(xmin[1,1:m], noisy = FALSE)))

}
