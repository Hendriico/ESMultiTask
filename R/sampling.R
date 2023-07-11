
#' Parallel tempering
#'
#' Parallel tempering algorithm to sample from any given (unnormalised) PDF
#'
#'
#' @param fn
#' @param temps.start
#' @param start.loc
#' @param n.chains
#' @param N
#' @param p
#' @param proposal.dist
#' @param walk.sd
#' @param bounded
#' @param lower
#' @param upper
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'


parallel_tempering_sample <- function(fn,
                                      temps.start = NULL,
                                      start.loc = NULL,
                                      n.chains = NULL,
                                      N = 1000,
                                      p = NULL,
                                      proposal.dist = NULL,
                                      walk.sd = NULL,
                                      bounded = FALSE,
                                      lower = NULL,
                                      upper = NULL,
                                      ...){



  if(is.null(start.loc) & is.null(p)){
    stop("Must specify at least one of start.loc or p")
  }

  if(is.null(p)) p <- ncol(start.loc)

  if(is.null(temps.start) & is.null(start.loc) & is.null(n.chains)) n.chains <- 4

  if(is.null(temps.start) & is.null(n.chains)) n.chains <- nrow(start.loc)

  if(is.null(start.loc) & is.null(n.chains)) n.chains <- length(temps.start)

  if(is.null(temps.start)) temps.start <- (100^(1/(n.chains-1)))^(0:(n.chains-1))




  if(is.null(walk.sd) & is.null(proposal.dist)){

    proposal.dist <- rnorm

    walk.sd <- (1)/(4*sqrt(max(temps.start)))*sqrt(temps.start)

    if(bounded) walk.sd <- walk.sd*(upper-lower)
  }



  if(is.null(start.loc)){
    start.loc <- lhs::maximinLHS(n.chains, p)
    if(bounded) start.loc <- start.loc*(upper-lower) + lower
  }

  v <- 100/n.chains
  t0 <- 1000/n.chains

  temps <- temps.start
  proposalold <- start.loc
  proposalold_value <- fn(proposalold, ...)^(1/temps)
  swapcount <- numeric(n.chains-1)
  swapaccept <- numeric(n.chains-1)

  PT_sample <- matrix(nrow = N, ncol = p)
  responses <- numeric(N)

  for(k in 1:N){

    proposal <- proposalold + matrix(proposal.dist(p*n.chains, sd = walk.sd), ncol = p)
    if(bounded){
      proposal[proposal < lower] <- lower + abs(rnorm(1, sd = walk.sd[1]))
      proposal[proposal > upper] <- upper - abs(rnorm(1, sd = walk.sd[1]))
    }


    proposaltrue <- fn(proposal, ...)
    proposal_value <- proposaltrue^(1/temps)
    accept <-  runif(n.chains) < proposal_value/proposalold_value
    accept <- replace(accept, is.na(accept), FALSE)
    proposalold[accept,] <- proposal[accept,]
    proposalold_value[accept] <- proposal_value[accept]

    swap <- sample(1:(n.chains-1), 1)
    swapcount[swap] <- swapcount[swap] + 1
    swap <- c(swap, swap+1)

    Ti <- max(2, swap[1])
    Ai <- swapaccept[Ti]/(swapcount[Ti]+0.01)
    Ai1 <- swapaccept[Ti-1]/(swapcount[Ti-1]+0.01)
    kt <- 1/v*t0/(k+t0)

    test <- #0.25/(Ai+0.01)*
      min(prod(fn(matrix(proposalold[swap,], nrow = 2), ...)^(1/temps[rev(swap)]))/prod(proposalold_value[swap]), 1)
    test <- replace(test, is.na(test),0)
    if(runif(1) < test){
      proposalold[swap, ] <- proposalold[rev(swap),]
      proposalold_value[swap] <-  fn(matrix(proposalold[swap,], nrow = 2),...)^(1/temps[swap])
      swapaccept[swap[1]] <- swapaccept[swap[1]] + 1
    }
    PT_sample[k,] <- proposalold[1,]
    responses[k] <- proposalold_value[1]

    temps[Ti] <- exp(kt*(Ai1-Ai)  + log(temps[Ti] - temps[Ti-1]) )+ temps[Ti-1]



  }

  return(list(PT_sample = PT_sample, responses = responses))

}
