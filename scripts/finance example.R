
library(quantmod)
library(ggplot2)
library(forecast)
library(doFuture)
library(doRNG)

toDate <- as.Date("2023-06-01")

getSymbols("^GSPC",src="yahoo",from="2017-01-01",to = toDate)
getSymbols("AMZN",src="yahoo",from="2017-01-01",to = toDate)
getSymbols("MSFT",src="yahoo",from="2017-01-01",to = toDate)
getSymbols("TSLA",src="yahoo",from="2017-01-01",to = toDate)
getSymbols("AAPL",src="yahoo",from="2017-01-01",to = toDate)
getSymbols("BTC-USD",src="yahoo",from="2017-01-01",to = toDate)
getSymbols("MU",src="yahoo",from="2017-01-01",to = toDate)
getSymbols("GOOGL",src="yahoo",from="2017-01-01",to = toDate)
getSymbols("SHOP",src="yahoo",from="2017-01-01",to = toDate)


modelfit <- auto.arima(GSPC$GSPC.Close, lambda = "auto")
modelfit1 <-auto.arima(AMZN$AMZN.Close, lambda = "auto")
modelfit2 <- auto.arima(MSFT$MSFT.Close, lambda = "auto")
modelfit3 <- auto.arima(TSLA$TSLA.Close, lambda = "auto")
modelfit4 <- auto.arima(AAPL$AAPL.Close, lambda = "auto")
modelfit5 <- auto.arima(`BTC-USD`$`BTC-USD.Close`, lambda = "auto")
modelfit6 <- auto.arima(MU$MU.Close, lambda = "auto")
modelfit7 <- auto.arima(GOOGL$GOOGL.Close, lambda = "auto")
modelfit8 <- auto.arima(SHOP$SHOP.Close, lambda = "auto", d = 1)



getsimulations <-  function(model, nsim = 10000, days = 60){
  buy_value <- model$x[length(model$x)]
  sim <- foreach::foreach(i = 1:nsim, .combine = "c", .packages = "forecast") %dorng% {
    close_value <- simulate(model, nsim = days)[days]
    close_value/buy_value
  }
  return(sim)
}


risk_score <- function(val){
  q <- quantile(val, prob=c(0.05, 0.95))
  
  q1 <- max(1-q[1],0)
  q2 <- q[2]-1
  if(q1 > 0.2){
    
    out <- ((q2-0.2)^2*40-(q1-0.2)^2*40)*exp(ifelse(q1>0.4, 0.4-q1, 0))^8*min(q2-q1+1,1)^40
  } else{
    out <- -(q2)*(0.2-q1)*100
  }
  return(out)
  
}

simres <- cbind(getsimulations(modelfit),getsimulations(modelfit1),getsimulations(modelfit2),
                getsimulations(modelfit3), getsimulations(modelfit4), getsimulations(modelfit5),
                getsimulations(modelfit6), getsimulations(modelfit7), getsimulations(modelfit8),
                rep(1.05, 100000))

f <- function(x){
  if(is.null(nrow(x))) x <- matrix(x, nrow=1)
  
  out <- numeric(nrow(x))
  
  for(i in 1:nrow(x)){
    x.temp <- x[i,]
    x.temp <- x.temp/sum(x.temp)

    out[i] <- risk_score(rowSums(t(t(simres)*x.temp)))
    
  }
  return(out)
}


results <- tibble::tibble()
upper <- 1
lower <- 0

set.seed(100)
registerDoFuture()
plan(multisession)

for(i in 1:20){
  
  Xinit <- lhs::maximinLHS(20,10)*(upper-lower)+lower
  tryCatch(
    {
      tic <- proc.time()
      res.EI <- ESMultiTask::multitask_EI(Xinit, f, 300, g=NULL, lower, upper, add = "both", optim.multi.start = 50)
      proc.time() - tic
      res <- ESMultiTask::entropy_search(Xinit, f, 300, lower, upper, g=NULL,  K=10000, M=120, particles = 6000, prop.EI = 1, updateGP = 1,update.particles = 3)
      
      
      for (n in seq(20, 300, by=5)) {
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
      
      for (n in seq(20, 300, by=5)) {
        
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

