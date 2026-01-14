data_gen_gcomp <- function(n_indiv, beta_ttt){
  
  x1 <- rnorm(n = n_indiv, mean = 1, sd = 0.5)
  x2 <- rnorm(n = n_indiv ,mean = 2, sd = 1.5)
  
  x3 <- rbinom(n = n_indiv, size = 1, prob = 0.4)
  x4 <- rbinom(n = n_indiv, size = 2, prob = 0.8)
  
  x_ttt <- rbinom(n = n_indiv, size = 1, prob = 0.5)
  
  beta_1 <- 4
  beta_2 <- 0.5
  
  beta_3 <- 1.5
  beta_4 <- -2.5
  
  beta_ttt <- beta_ttt
  
  lp <- beta_ttt*x_ttt + beta_1*x1 + beta_2*x2 + beta_3*x3 + beta_4*x4
  
  y_continuous <- rnorm(n = n_indiv, mean = lp)
  y_count <- rpois(n = n_indiv, lambda = exp(lp))
  
  data <- data.frame(x1 = x1,
                     x2 = x2,
                     x3 = x3,
                     x4 = x4,
                     ttt = x_ttt,
                     y_continuous = y_continuous,
                     y_count = y_count)
  
  return(data)
  
}
