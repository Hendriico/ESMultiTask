
#' Optimal computing budget allocation algorithm
#'
#' @param design
#' @param sample.means
#' @param sample.sd
#' @param ra
#'
#' @return
#' @export
#'
#' @examples
OCBA <- function(
  design,
  sample.means = design$sample.mean,
  sample.sd = design$sample.sd,
  ra
){
  d <- ncol(design)-2
  m <- which.min(sample.means)

  sample.mean.min <- sample.means[m]
  sample.means <- sample.means - sample.mean.min

  N <- rep(1, length(sample.means))

  default <- ifelse(m==1, 2, 1)

  for (i in (default+1):length(sample.means)) {
    if(i !=m) N[i] <- N[default]/((sample.sd[default]/sample.means[default])/(sample.sd[i]/sample.means[i]))^2
  }

  N[m] <- 0

  for (i in 1:length(sample.means)) {
    if(i !=m){
      N[m] <- N[m] + N[i]^2/sample.sd[i]^2

    }
  }
  N[m] <- sqrt(N[m])*sample.sd[m]

  N <- N/(sum(N)/ra)
  roundN <- round(N)
  npos <- which(roundN>0)
  if(sum(roundN) < ra){
    extra <- ra - sum(roundN)
    nextra <- order(N, decreasing = T)[1:extra]
    npos <- c(rep(npos, round(N)[npos]), nextra)
  } else{
    npos <- rep(npos, round(N)[npos])
  }

  res <- data.matrix(design[npos,1:d])
  res.sd <- design[npos, d+2]
  colnames(res) <- NULL
  return(list(X = res, sd = res.sd))
}
