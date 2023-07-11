
library(doFuture)



langer <- function(X,noisy = TRUE, sigma = 0.1, m=5, cvec, A)
{
  X <- X*10
  if(is.null(nrow(X))) X <- matrix(X, nrow=1)
  res <- numeric(nrow(X))
  
  d <- ncol(X)
  
  if (missing(cvec)) {
    if (m == 5){
      cvec <- c(1,2,5,2,3)
    }
    else {
      stop('Value of the m-dimensional vector cvec is required.')
    }
  }
  
  if (missing(A)) {
    if (m==5 && d==2) {
      A <- matrix(c(3,5,5,2,2,1,1,4,7,9),5,2,byrow=TRUE)
    } else if(m==5 & d%%2==0){
      A <- matrix(rep(c(3,5,2,1,7,5,2,1,4,9), d/2),5,d,byrow=FALSE)
    }else {
      stop('Value of the (mxd)-dimensional matrix A is required.')
    }
  }
  
  for (i in 1:nrow(X)) {
    xx <- X[i,]
    xxmat <- matrix(rep(xx,times=m), m, d, byrow=TRUE)    
    inner <- rowSums((xxmat-A[,1:d])^2)	
    res[i] <- sum(cvec * exp(-inner/pi) * cos(pi*inner))
  }
  if(noisy) res <- res + rnorm(nrow(X), sd = sigma)
  return(res)
}



f <- langer

results <- tibble::tibble()
upper <- 1
lower <- 0

set.seed(10)
registerDoFuture()
plan(multisession)

for(i in 1:20){
  
  Xinit <- lhs::maximinLHS(20,4)*(upper-lower)+lower
  tryCatch(
    {
      tic <- proc.time()
      res.EI <- ESMultiTask::multitask_EI(Xinit, f, 120, g=NULL, lower, upper)
      res <- ESMultiTask::entropy_search(Xinit, f, 120, lower, upper, g=NULL,  K=8000, M=120, particles = 8000, prop.EI = 1, update.particles = 2, direction = "both")#, walk.sd = c(0.05000000, 0.06667607, 0.08891397, 0.11856869, 0.15811388, 0.21084825, 0.28117066, 0.37494710, 0.50000000)/2)
      
      
      for (n in seq(20, 120, by=5)) {
        new <- tibble::tibble(i =i, N = n, algorithm = "ES",
                              # xmin = f(ESMultiTask::find_GP_optim(res$X[1:n,], res$y[1:n], lower = lower, upper = upper, optimisation.type = "optim", lhs.size = 10000)$xmin, 0),
                              # xmax = f(ESMultiTask::find_GP_optim(res$X[1:n,], -res$y[1:n], lower = lower, upper = upper, optimisation.type = "optim", lhs.size = 10000)$xmin, 0))
                              xmin = min(f(res$X[1:n,], noisy = FALSE)),
                              xmax = max(f(res$X[1:n,], noisy = FALSE))
                              # xmin = min(res$y[1:n]),
                              # xmax = max(res$y[1:n])
        )
        
        results <- rbind(results, new)
      }
      
      for (n in seq(20, 120, by=5)) {
        
        new <- tibble::tibble(i = i, N = n, algorithm = "EI",
                              # xmin = f(ESMultiTask::find_GP_optim(res.EI$X[1:n,], res.EI$y[1:n], lower = lower, upper = upper, optimisation.type = "optim", lhs.size = 10000)$xmin, 0),
                              # xmax = f(ESMultiTask::find_GP_optim(res.EI$X[1:n,], -res.EI$y[1:n], lower = lower, upper = upper, optimisation.type = "optim", lhs.size = 10000)$xmin, 0)
                              xmin = min(f(res.EI$X[1:n,], noisy = FALSE)),
                              xmax = max(f(res.EI$X[1:n,], noisy = FALSE))
                              # xmin = min(res.EI$y[1:n]),
                              # xmax = max(res.EI$y[1:n])
        )
        results <- rbind(results, new)
      }
      
    },
    error = function(e) print(e)
  )
  
  
  print(paste("i = ", i, " time was ", proc.time()[3]-tic[3]))
  
}

library(ggplot2)
library(dplyr)

results %>%
  group_by(algorithm, N) %>%
  summarise( xmin = mean(xmin),
             xmax = mean(xmax)) %>%
  ggplot(aes(N, xmin, col = algorithm)) +
  geom_line() +
  geom_line(aes(y=xmax))

results %>%
  ggplot(aes(algorithm, xmin, col = algorithm)) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~N)

write.csv(results, "multitask_res_4d.csv", row.names = FALSE)

results <- vroom::vroom("multitask_res1.csv")
