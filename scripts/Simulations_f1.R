library(doFuture)
f <- function(x){
  if(is.null(nrow(x))) x <- matrix(x, ncol=2)

  (5*x[,1]*sin(x[,2]*2)*exp(-(x[,1])^6-(x[,2]+1.5)^6 ) + 2*cos(x[,2]/5)*exp(x[,1]/5-(x[,2]-1)^6)-0.2*x[,1] + 0.2*x[,2]^2- 0.4*x[,2])-1.5
}



results <- tibble::tibble()
upper <- 3
lower <- -3

set.seed(10)
registerDoFuture()
plan(multisession)

for(i in 1:20){
  
  Xinit <- lhs::maximinLHS(20,2)*(upper-lower)+lower
  tryCatch(
    {
      tic <- proc.time()
      res.EI <- ESMultiTask::multitask_EI(Xinit, f, 100, g=1e-06, lower, upper, add = "both")
      res <- ESMultiTask::entropy_search(Xinit, f, 100, lower, upper, g=1e-06,  K=5000, M=80, particles = 5000, prop.EI = 1)
      
      
      for (n in seq(20, 100, by=5)) {
        new <- tibble::tibble(i =i, N = n, algorithm = "ES",
                              # xmin = f(ESMultiTask::find_GP_optim(res$X[1:n,], res$y[1:n], lower = lower, upper = upper, optimisation.type = "optim", lhs.size = 10000)$xmin, 0),
                              # xmax = f(ESMultiTask::find_GP_optim(res$X[1:n,], -res$y[1:n], lower = lower, upper = upper, optimisation.type = "optim", lhs.size = 10000)$xmin, 0))
                              # xmin = min(f(res$X[1:n,])),
                              # xmax = max(f(res$X[1:n,]))
                              xmin = min(res$y[1:n]),
                              xmax = max(res$y[1:n])
        )
        
        results <- rbind(results, new)
      }
      
      for (n in seq(20, 100, by=5)) {
        
        new <- tibble::tibble(i = i, N = n, algorithm = "EI",
                              # xmin = f(ESMultiTask::find_GP_optim(res.EI$X[1:n,], res.EI$y[1:n], lower = lower, upper = upper, optimisation.type = "optim", lhs.size = 10000)$xmin, 0),
                              # xmax = f(ESMultiTask::find_GP_optim(res.EI$X[1:n,], -res.EI$y[1:n], lower = lower, upper = upper, optimisation.type = "optim", lhs.size = 10000)$xmin, 0)
                              # xmin = min(f(res.EI$X[1:n,])),
                              # xmax = max(f(res.EI$X[1:n,]))
                              xmin = min(res.EI$y[1:n]),
                              xmax = max(res.EI$y[1:n])
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

