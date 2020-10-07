################################################

#  The function below determines how well predict 
# function forecast n steps ahead for 
# a time series simulated from arma(p q)

###############################################


guess_arima <- function(min_pq,max_pq,p,q,training_size) { 
 
###################### Make sure there is a casual solution ########################### 
  p_roots <- c(1:p)
  while (min(p_roots)<=1) {
    phi <- runif(p, min=-1, max=1) # randomly generated phi for ARMA(p q)
    p_roots <- polyroot(c(1,-phi))
    p_roots <- Re(p_roots)^2+Im(p_roots)^2
  }
  
  t_roots <- c(1:q)
  while (min(t_roots)<=1) {
    theta <- runif(q, min=-1, max=1)# randomly generated theta for ARMA(p q)
    t_roots <- polyroot(c(1,-theta))
    t_roots <- Re(t_roots)^2+Im(t_roots)^2
  }
  
##################### Simulate Data for ARMA(p q) ###################################
  
  sim <- arima.sim(n = (training_size+10), list(ar=phi, ma = theta),sd = sqrt(1))
  training <- sim[1:training_size] #Trainig Data
  predict_d <- sim[(training_size+1):(training_size+10)] #Test Data h = 10
 
############## Find fit model with smallest AIC & Calculate MSE ####################
  guess_p = 1
  guess_q = 1
  min_aic = 10^6
  for (i in min_pq:max_pq) {
    for (j in min_pq:max_pq) {
      m_aic = try(AIC(arima(x = training, order = c(i, 0, j))),silent = T) #get AIC for each group of p q
      if ((typeof(m_aic)=="double")&&(m_aic < min_aic)) {
        min_aic = m_aic
        guess_p = i
        guess_q = j
      }
    }
  }

  if (min_aic < 10^6) {
    estimate <- predict(arima(training, order=c(guess_p, 0,guess_q)), n.ahead=10)# Predict 10 steps ahead using predict function 
    mse1 <- mean((predict_d[1] - estimate$pred[1])^2) # Calculate mean square error for predictin 1 step ahead
    mse2 <- mean((predict_d[1:2] - estimate$pred[1:2])^2)# Calculate mean square error for predictin 2 step ahead
    mse5 <- mean((predict_d[1:5] - estimate$pred[1:5])^2)# Calculate mean square error for predictin 5 step ahead
    mse10 <- mean((predict_d[1:10] - estimate$pred[1:10])^2)# Calculate mean square error for predictin 10 step ahead
    
    correct_pq = ((p==guess_p)&&(q==guess_q)) #Determine whether p q we guess are equal to actual p q
    
    return (c(correct_pq,mse1,mse2,mse5,mse10))
  } else {
    return (guess_arima(min_pq,max_pq,p,q,training_size))
  }

}


########## Run guess_arima n times to see how many times p q we guess = actual p q ################

experiment_arima <- function(n) {
  
  count500 = 0
  mse500_1 = vector()
  mse500_2 = vector()
  mse500_5 = vector()
  mse500_10 = vector()
  
  count20 = 0
  mse20_1 = vector()
  mse20_2 = vector()
  mse20_5 = vector()
  mse20_10 = vector()
  for (i in 1:n) {
    min_pq = 1
    max_pq = 6
    p <- round(runif(1, min=min_pq, max=max_pq))
    q <- round(runif(1, min=min_pq, max=max_pq))
    result_500 = guess_arima(min_pq,max_pq,p,q,500)
    result_20 = guess_arima(min_pq,max_pq,p,q,20)
    
    count500 = count500 + result_500[1]
    mse500_1 <- c(mse500_1,result_500[2])
    mse500_2 <- c(mse500_2,result_500[3])
    mse500_5 <- c(mse500_5,result_500[4])
    mse500_10 <- c(mse500_10,result_500[5])
    
    count20 = count20 + result_20[1]
    mse20_1 <- c(mse20_1,result_20[2])
    mse20_2 <- c(mse20_2,result_20[3])
    mse20_5 <- c(mse20_5,result_20[4])
    mse20_10 <- c(mse20_10,result_20[5])
  }
  
  writeLines("count500: ")
  print(count500)
  writeLines("mse500_1: ")
  print(mse500_1)
  writeLines("mse500_2: ")
  print(mse500_2)
  writeLines("mse500_5: ")
  print(mse500_5)
  writeLines("mse500_10: ")
  print(mse500_10)
  
  writeLines("count20: ")
  print(count20)
  writeLines("mse20_1: ")
  print(mse20_1)
  writeLines("mse20_2: ")
  print(mse20_2)
  writeLines("mse20_5: ")
  print(mse20_5)
  writeLines("mse20_10: ")
  print(mse20_10)
}

