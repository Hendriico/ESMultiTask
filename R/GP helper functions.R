#' Make predictions from a Gaussian process, under different packages
#'
#' @param gpi
#' @param X
#' @param model.type
#'
#' @return
#' @export
#'
#' @examples

predGP_out <- function(gpi, X, model.type = "laGP"){

  if(is.null(nrow(X))) X <- matrix(X, nrow = 1)
  if(model.type == "laGP"){
    yout <- laGP::predGPsep(gpi, X, lite = TRUE)$mean
  } else if(model.type == "robust"){
    yout <- RobustGaSP::predict.rgasp(gpi, X)$mean
  } else if(model.type == "km"){
    yout <- DiceKriging::predict.km(gpi, X, "SK", se.compute = FALSE, light.return = TRUE)$mean
  }
  return(yout)
}



#' Find the optimum point of a Gaussian Process, given the design
#'
#' @param X
#' @param y
#' @param model.type
#' @param correlation
#' @param lower
#' @param upper
#' @param optimisation.type
#' @param optim.start
#' @param lhs.size
#'
#' @return
#' @export
#'
#' @examples

find_GP_optim <- function(X,
                          y,
                          model.type = "laGP",
                          correlation = "gauss",
                          lower = 0,
                          upper = 1,
                          optimisation.type = "optim",
                          optim.start = 10,
                          lhs.size = 1000){

  if(is.null(nrow(X))) X <- matrix(X, ncol=1)
  y.mod <- (y-mean(y))/sd(y)

  if(!optimisation.type %in% c("optim", "lhs")){
    stop("optimisation.type has to be either 'optim' or 'lhs'")
  }


  if(model.type == "laGP" & correlation != "gauss"){
    stop("Only 'gauss' correlation function is supported with laGP")
  }

  if(!correlation %in% c("gauss", "matern_3_2", "matern_5_2", "matern3_2", "matern5_2")){
    stop("correlation has to be either 'gauss', 'matern_3_2' or 'matern_5_2'")
  }


  d <- ncol(X)

  if(model.type == "laGP"){
    da <- laGP::darg(d = NULL, lhs::randomLHS(3000, d)*(upper-lower) + lower)
    ga <- laGP::garg(list(mle = TRUE, max = var(y.mod)), y.mod)
    gpi <- laGP::newGPsep(X, y.mod, d = da$start, g = ga$start, dK = TRUE)
    theta <- laGP::mleGPsep(gpi, param = "both", tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))$theta

  } else if(model.type == "robust"){
    sink("NUL")
    if(correlation == "gauss"){
      gpi <- RobustGaSP::rgasp(X, y, nugget.est = TRUE, kernel_type = "pow_exp", alpha = rep(2, d))
    } else{
      gpi <- RobustGaSP::rgasp(X, y, nugget.est = TRUE, kernel_type = correlation, num_initial_values = 2*d)
      if(gpi@nugget < 0.001) gpi <- RobustGaSP::rgasp(X, y, nugget = 0.001, kernel_type = correlation, num_initial_values = 2*d)
      theta <- c(1/gpi@beta_hat, gpi@nugget)
    }
    sink()
  } else if(model.type == "km"){
    gpi <- DiceKriging::km(design = X, response = y, covtype = correlation, nugget.estim = TRUE, control=list(trace=FALSE), multistart = 5*d)
    colnames(X) <- colnames(gpi@X)
    theta <- c(gpi@covariance@range.val, gpi@covariance@nugget)
  }else {
    stop("model.type has to be either robust, km or laGP")
  }




  if(optimisation.type == "lhs"){
    if(lhs.size > 2000){
      test <- rbind(lhs::randomLHS(lhs.size,d)*(upper-lower) + lower, X)
    } else{
      test <- rbind(lhs::maximinLHS(lhs.size,d)*(upper-lower) + lower, X)

    }

    if(model.type == "laGP"){
      yout <- laGP::predGPsep(gpi, test, lite = TRUE)$mean
    } else if(model.type == "robust"){
      yout <- RobustGaSP::predict.rgasp(gpi, test)$mean
    } else if(model.type == "km"){
      yout <- DiceKriging::predict.km(gpi, test, "SK", se.compute = FALSE, light.return = TRUE)$mean
    }
    xmax <- test[which.max(yout),]
    xmin <- test[which.min(yout),]
    res = list(xmax =xmax, xmin = xmin)
  } else{

    start <- rbind(X, lhs::maximinLHS(optim.start, d))

    optim.res <- matrix(nrow = nrow(start), ncol = d+1)



    for(i in 1:nrow(start)){
      out <- optim(start[i,],
                   predGP_out,
                   method = "L-BFGS-B",
                   lower = lower,
                   upper = upper,
                   gpi = gpi,
                   model.type = model.type
      )
      optim.res[i,] <- c(out$par, out$value)
    }

    xmin = optim.res[which.min(optim.res[,d+1]),1:d]
    res = list(xmin = xmin)

  }

  if(model.type == "laGP"){
    laGP::deleteGPsep(gpi)
    }

  return(res)
}

