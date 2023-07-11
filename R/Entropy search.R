



#' GP distance
#'
#' Calculates the Gaussian correlation matrix between two matrices
#'
#' @param X1
#' @param X2
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
GPdistance <- function(X1,
                       X2 = NULL,
                       theta){

  if(is.null(ncol(X1))) X1 <- matrix(X1, ncol=1)
  if(is.null(X2)) X2 <- X1
  d <- ncol(X1)

  D <- 0

  for(i in 1:d){
    D <- D + laGP::distance(X1[,i], X2[,i])/theta[i]
  }

  return(D)
}

#' Covariance matrix
#'
#' Function to calculate a variance-covariance matrix for a Gaussian process
#' Supported correlation functions are Gaussian, matern3_2 and matern5_2
#'
#' @param X1
#' @param X2
#' @param theta
#' @param correlation
#'
#' @return
#' @export
#'
#' @examples

covMatrix <- function(X1,
                      X2 = NULL,
                      theta,
                      correlation = "gauss"){

  if(is.null(ncol(X1))) X1 <- matrix(X1, ncol=1)
  if(is.null(X2)) X2 <- X1
  d <- ncol(X1)



  if(correlation == "gauss"){
    D <- 0

    for(i in 1:d){
      D <- D + laGP::distance(X1[,i], X2[,i])/theta[i]
    }
    out <- exp(-D)
  } else if(correlation == "matern_3_2" | correlation == "matern3_2"){

    out <- 1
    for(i in 1:d){
      D <- sqrt(laGP::distance(X1[,i], X2[,i]))
      out <- out*(1+sqrt(3)*D/theta[i])*exp(-sqrt(3)*D/theta[i])
    }


  }else if(correlation == "matern_5_2" | correlation == "matern5_2"){

    out <- 1
    for(i in 1:d){
      D <- laGP::distance(X1[,i], X2[,i])
      Dsqrt <- sqrt(D)
      out <- out*(1+sqrt(5)*Dsqrt/theta[i]+(5/3)*(Dsqrt/theta[i])^2)*exp(-sqrt(5)*Dsqrt/theta[i])
    }


  }else{
    stop("correlation has to be either 'gauss', 'matern3_2'/'matern_3_2' or 'matern5_2'/'matern_5_2'")
  }


  return(out)
}


#' GP update
#'
#' Calculates mean and variance after adding in a point to a Gaussian Process
#'
#'
#' @param X
#' @param X.cur
#' @param y.cur
#' @param X.add
#' @param y.add
#' @param theta
#' @param mean.cur
#' @param Sigma
#' @param sigma
#' @param mean.add
#'
#' @return
#' @export
#'
#' @examples
GPupdate <- function(X,
                     X.cur,
                     y.cur = NULL,
                     X.add,
                     y.add = NULL,
                     theta,
                     mean.cur,
                     Sigma,
                     sigma = NULL,
                     mean.add = NULL,
                     correlation = "gauss"){

  if(is.null(nrow(X.add))) X.add <- matrix(X.add, nrow = 1)
  g <- theta[length(theta)]

  if(is.null(sigma)){
    if(is.null(y.cur)){
      stop("y.cur must be given if sigma is NULL")
    }else{
      sigma <- drop(t(y.cur) %*%  solve(covMatrix(X.cur, X.cur, theta, correlation) + diag(g, nrow = nrow(X.cur))) %*% y.cur / (length(y.cur)))
    }

  }

  SxX <-  (covMatrix(X, X.add, theta, correlation) -
             covMatrix(X, X.cur, theta, correlation)%*%
             solve( covMatrix(X.cur, X.cur, theta, correlation)+ diag(g, nrow = nrow(X.cur)) ) %*%
             covMatrix(X.cur, X.add, theta, correlation))

  SXX <-  solve( covMatrix(X.add, X.add, theta, correlation)  -
                   covMatrix(X.add, X.cur, theta, correlation)%*%
                   solve(covMatrix(X.cur, X.cur, theta, correlation) + diag(g, nrow = nrow(X.cur)) ) %*%
                   covMatrix(X.cur, X.add, theta, correlation)  + g)

  addmean <-  SxX %*% SXX
  newSigma <- Sigma - sigma*(SxX %*%SXX%*% t(SxX))
  res <- list(addmean = addmean, Sigma = newSigma)

  if(!is.null(y.add) & !is.null(mean.add)){
    newmean <- mean.cur + addmean *(y.add - mean.add)
    list$mean <- newmean
  }
  return(res)
}



#' Entropy Search
#'
#' @param X
#' @param f
#' @param N
#' @param lower
#' @param upper
#' @param M
#' @param K
#' @param particles
#' @param g
#' @param update.particles
#' @param temps.start
#' @param discrete.type
#' @param prop.EI
#' @param direction
#' @param updateGP
#' @param track.progress
#' @param rePT
#' @param walk.sd 
#'
#' @return
#' @export
#'
#' @examples
entropy_search <- function(X,
                           f,
                           N,
                           lower = 0,
                           upper = 1,
                           M = 100,
                           K=1000,
                           particles = 5000,
                           g = NULL,
                           update.particles = 5,
                           temps.start = NULL,
                           walk.sd = NULL,
                           discrete.type = "EI",
                           prop.EI = 0.8,
                           direction = "both",
                           updateGP = 1,
                           track.progress = TRUE,
                           rePT= TRUE){

  if(is.null(ncol(X))) X <- matrix(X, ncol = 1)
  m <- ncol(X)
  d <- ifelse(direction == "both", 2*m, m)
  t <- round(4*sqrt(m)) + 1
  if(is.null(temps.start) | length(temps.start) != t) temps.start <- (100^(1/(t-1)))^(0:(t-1))
  if(is.null(walk.sd) | length(walk.sd) != t) walk.sd <- (1)/(4*sqrt(max(temps.start)))*sqrt(temps.start)#*(upper-lower)
  v <- 100/t
  t0 <- 1000/t

  if(track.progress){
    pb <- progress::progress_bar$new(
      format = " Evaluating [:bar] :percent Time left: :eta Elapsed :elapsedfull",
      total = N-nrow(X), clear = FALSE)
  }

 

  y <- f(X)
  y.mod <- (y-mean(y))/sd(y)

  da <- laGP::darg(d = NULL, X = lhs::randomLHS(3000, m)*(upper-lower)+lower) # Setting a prior for theta
  if(is.null(g)){
    ga <- laGP::garg(list(mle = TRUE, max = var(y.mod)), y.mod)
    gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = ga$start, dK = TRUE)
    theta <- laGP::mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))$theta
  } else{
    gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = g, dK = TRUE)
    theta <- c(laGP::mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$d, g)
  }

  predcur <- laGP::predGPsep(gpi, X, lite = TRUE, nonug = TRUE)

  if(direction %in% c("max", "both") ){
    curmax <- max(predcur$mean)
  }

  if(direction %in% c("min", "both")){
    curmin <- min(predcur$mean)
  }
  EI.M <- ceiling(M * prop.EI)
  xparticles <- matrix(nrow = particles, ncol = m)
  xdiscrete <- matrix(nrow = EI.M, ncol = m)
  responses <- numeric(particles)

  if(direction == "both"){
    EI.Mmax <- round(0.5*EI.M)
    EI.Mmin <- EI.M - EI.Mmax
    particles.Mmax <- round(0.5*particles)
    particles.Mmin <- particles - particles.Mmax
  } else if(direction == "max"){
    EI.Mmax <- EI.M
    EI.Mmin <- 0
    particles.Mmax <- particles
    particles.Mmin <- 0
  } else if(direction == "min"){
    EI.Mmin <- EI.M
    EI.Mmax <- 0
    particles.Mmax <- 0
    particles.Mmin <- particles
  }


  if(direction %in% c("max", "both") ){

    maxsample <- ESMultiTask::parallel_tempering_sample(fn = ESMultiTask::EI,
                                                      temps.start = temps.start,
                                                      N = particles.Mmax,
                                                      p = m,
                                                      bounded = TRUE,
                                                      lower = lower,
                                                      upper = upper,
                                                      model = gpi,
                                                      fmax = curmax)

    xparticles[1:particles.Mmax, ] <- maxsample$PT_sample
    responses[1:particles.Mmax] <- maxsample$responses

    weights <- responses[1:particles.Mmax]/sum(responses[1:particles.Mmax])
    xdiscrete[1:EI.Mmax,] <- xparticles[sample(1:particles.Mmax, EI.Mmax, prob = weights),]

  }
  if(direction %in% c("min", "both") ){

    minsample <- ESMultiTask::parallel_tempering_sample(fn = ESMultiTask::EI,
                                                      temps.start = temps.start,
                                                      N = particles.Mmin,
                                                      p = m,
                                                      bounded = TRUE,
                                                      lower = lower,
                                                      upper = upper,
                                                      model = gpi,
                                                      fmin = curmin,
                                                      maximisation = FALSE)

    xparticles[(1:particles.Mmin)+particles.Mmax, ] <- minsample$PT_sample
    responses[(1:particles.Mmin)+particles.Mmax] <- minsample$responses

    weights <- responses[(1:particles.Mmin)+particles.Mmax]/sum(responses[(1:particles.Mmin)+particles.Mmax])
    xdiscrete[(1:EI.Mmin)+EI.Mmax,] <- xparticles[sample((1:particles.Mmin)+particles.Mmax, EI.Mmin, prob = weights),]
  }
  if(prop.EI < 1){
    xdiscrete <- rbind(xdiscrete, lhs::maximinLHS(M-EI.M, m)*(upper-lower) + lower)
  }


  for (i in nrow(X):(N-1)) {
    pred <- laGP::predGPsep(gpi, xdiscrete, nonug = TRUE)
    Sigma <- pred$Sigma
    sigma <- NULL

    entropies <- numeric(M)

    entropies <- foreach::foreach(j= 1:M, .packages = "ESMultiTask", .combine = "c") %dopar% {

      prednew <- ESMultiTask::GPupdate(X = xdiscrete, X.cur = X, y.cur = y.mod, X.add = xdiscrete[j,], theta = theta, mean.cur = pred$mean, Sigma = Sigma, sigma = sigma)

      simres <- as.data.frame(ESMultiTask::xdist(xdiscrete, pred$mean, prednew$Sigma, prednew$addmean, Sigma[j,j], direction, K))
      colnames(simres) <- paste0("x", 1:d)
      simres <- simres %>%
        dplyr::group_by(.dots = colnames(simres)) %>%
        dplyr::summarise(counts = dplyr::n()/K, .groups = 'drop')

      -sum(simres$counts*log(simres$counts))
    }

    xnew <- matrix(xdiscrete[which.min(entropies),], nrow = 1)
    ynew <- f(xnew)

    X <- rbind(X, xnew)
    y <- c(y, ynew)
    y.mod <- (y-mean(y))/sd(y)

    if(is.null(g)){
      ga <- laGP::garg(list(mle = TRUE, max = var(y.mod)), y.mod)
      gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = ga$start, dK = TRUE)
      theta <- laGP::mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))$theta
    } else{
      gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = g, dK = TRUE)
      theta <- c(laGP::mleGPsep(gpi, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)$d, g)
    }

    if(nrow(X) != N){

      predcur <- laGP::predGPsep(gpi, X, lite = TRUE, nonug = TRUE)


      if(direction %in% c("max", "both") ){
        curmax <- max(predcur$mean)
      }

      if(direction %in% c("min", "both")){
        curmin <- min(predcur$mean)
      }

      if((nrow(X) + 2) %% update.particles == 0 & rePT){

        if(direction %in% c("max", "both") ){

          maxsample <- parallel_tempering_sample(fn = ESMultiTask::EI,
                                                 temps.start = temps.start,
                                                 N = particles.Mmax,
                                                 p = m,
                                                 bounded = TRUE,
                                                 lower = lower,
                                                 upper = upper,
                                                 model = gpi,
                                                 fmax = curmax)

          xparticles[1:particles.Mmax, ] <- maxsample$PT_sample
          responses[1:particles.Mmax] <- maxsample$responses

          weights <- responses[1:particles.Mmax]/sum(responses[1:particles.Mmax])
          xdiscrete[1:EI.Mmax,] <- xparticles[sample(1:particles.Mmax, EI.Mmax, prob = weights),]

        }
        if(direction %in% c("min", "both") ){

          minsample <- parallel_tempering_sample(fn = ESMultiTask::EI,
                                                 temps.start = temps.start,
                                                 N = particles.Mmin,
                                                 p = m,
                                                 bounded = TRUE,
                                                 lower = lower,
                                                 upper = upper,
                                                 model = gpi,
                                                 fmin = curmin,
                                                 maximisation = FALSE)

          xparticles[(1:particles.Mmin)+particles.Mmax, ] <- minsample$PT_sample
          responses[(1:particles.Mmin)+particles.Mmax] <- minsample$responses

          weights <- responses[(1:particles.Mmin)+particles.Mmax]/sum(responses[(1:particles.Mmin)+particles.Mmax])
          xdiscrete[(1:EI.Mmin)+EI.Mmax,] <- xparticles[sample((1:particles.Mmin)+particles.Mmax, EI.Mmin, prob = weights),]
        }



      } else{
        xparticles <- xparticles + rnorm(m*particles, sd = walk.sd[1])
        xparticles[xparticles < lower] <- lower
        xparticles[xparticles > upper] <- upper
        if(direction %in% c("max", "both") ){
          EI.new <- ESMultiTask::EI(matrix(xparticles[1:particles.Mmax,], nrow = particles.Mmax), gpi, curmax)
          weights <- EI.new/sum(EI.new)
          xparticles[1:particles.Mmax,] <- xparticles[sample(1:particles.Mmax, replace = TRUE, prob = weights),]
          xpartun <- matrix(unique(xparticles[(1:particles.Mmax),]), ncol = m)
          unt <- nrow(xpartun)
          if(unt < EI.Mmax){
            xdiscrete[1:EI.Mmax,] <- xpartun[sample(1:unt, EI.Mmax, replace = TRUE),]
          } else{
            xdiscrete[1:EI.Mmax,] <- xpartun[sample(1:unt, EI.Mmax),]
          }


        }

        if(direction %in% c("min", "both")){
          EI.new <- ESMultiTask::EI(matrix(xparticles[(1:particles.Mmin) + particles.Mmax,], nrow = particles.Mmin), gpi, fmin = curmin, maximisation = FALSE)
          weights <- EI.new/sum(EI.new)
          xparticles[(1:particles.Mmin) + particles.Mmax,] <- xparticles[sample((1:particles.Mmin) + particles.Mmax, replace = TRUE, prob = weights),]
          xpartun <- matrix(unique(xparticles[(1:particles.Mmin)+particles.Mmax,]), ncol = m)
          unt <- nrow(xpartun)
          if(unt < EI.Mmin){
            xdiscrete[(1:EI.Mmin)+EI.Mmax,] <- xpartun[sample(1:unt, EI.Mmin, replace = TRUE),]

          } else{
            xdiscrete[(1:EI.Mmin)+EI.Mmax,] <- xpartun[sample(1:unt, EI.Mmin),]

          }

        }

      }
      if(prop.EI < 1){
        xdiscrete[(EI.M+1):M,] <- lhs::maximinLHS(M-EI.M, m)*(upper-lower) + lower
      }

    }

    if(track.progress){
      pb$tick()
    }

  }
  res <- list(X=X, y = y)

  return(res)
}




