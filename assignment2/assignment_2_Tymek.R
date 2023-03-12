rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# options(scipen=999)
options(scipen=0)
# dev.off()

library(MASS)
library(dplyr)
# library(tsibble)
library(forecast)
# library(matlib)
library(nlme)


######## 1. STABILITY ##################
n = 200
parameters_1 = c("phi_1" = 0.8, "theta_1" = 0.8, "theta_2" = -0.5, "std_noise" = 0.4)
simulations_number = 10
results = matrix(NA, nrow = n, ncol = simulations_number)
for (i in 1:simulations_number){
  results[,i] = arima.sim(n = n, list(ar = parameters_1[1], ma = parameters_1[2:3]),
                          sd = parameters_1[4])
}
average_model <- apply(results, MARGIN = 1, mean)
acf_plot <- acf(average_model)
plot(acf_plot)

MAX_LAG = 20
par(mfrow=c(1,1))
matplot(results, type = "l", col = "grey")
lines(average_model, col = "red", lwd = 2)

acf_average <- acf(average_model, na.actioin = na.pass, plot=FALSE,
                   lag.max = MAX_LAG)
pacf_average <- pacf(average_model, na.actioin = na.pass, plot=FALSE,
                     lag.max = MAX_LAG)

acf_full <- acf(results, na.action=na.pass, plot=FALSE, lag.max = MAX_LAG)
acf_matrix <- t(apply(acf_full$acf, MARGIN = 1, diag))

pacf_full <- pacf(results, na.action=na.pass, plot=FALSE, lag.max = MAX_LAG)
pacf_matrix <- t(apply(pacf_full$acf, MARGIN = 1, diag))

par(mfrow=c(1,2))
value = 0.15
matplot(acf_matrix, type = "l", col = "grey")
abline(h = 0, col = "black")
abline(h = value, col = "blue", lty = 2)
abline(h = -value, col = "blue", lty = 2)
lines(acf_average$acf, col = "red", lwd = 2)

value = 0.15
matplot(pacf_matrix, type = "l", col = "grey")
abline(h = 0, col = "black")
abline(h = value, col = "blue", lty = 2)
abline(h = -value, col = "blue", lty = 2)
lines(pacf_average$acf, col = "red", lwd = 2)

autoplot(fit) # fits the roots in the complex unit
# https://otexts.com/fpp2/arima-r.html

######## 2. APARTMENTS ###############
df <- read.table("A2_sales.txt", sep="",header=TRUE)
df <- df %>% 
  mutate(year = as.numeric(substr(Quarter, 1, 4)),
         quarter = as.numeric(substr(Quarter, 6, 6)),
         index = seq(1, dim(df)[1])) %>% 
  mutate(time = as.numeric(year + (quarter - 1)/4)) 
mu = 2070
period = 4
noise_variance = 36963

par(mfrow=c(1,1))
plot(df$index, df$Sales - mu)


parameters_2 = c("phi_1" = 1.04, "phi_2" = 0.2, "Phi_1" = 0.86,
               "Theta_1" = -0.42, "mu" = mu)
model <- arima(df$Sales, order = c(2, 0, 0),
               seasonal = list(order = c(1, 0, 1), period = period),
               include.mean = TRUE,
               fixed = parameters_2)
predictions <- predict(model, n.ahead = 2, se.fit = TRUE)

model$sigma2

model_our_sigma <- model
model_our_sigma$sigma2 <- noise_variance
str(model_our_sigma)

predictions_our_sigma <- predict(model_our_sigma, n.ahead = 2, se.fit = TRUE)

# ? how to set sigma in the arima function?

# using forecast package
df$Sales %>%
  Arima(order = c(2, 0, 0),
        seasonal = list(order = c(1, 0, 1), period = period),
        include.mean = TRUE,
        fixed = parameters_2) %>%
  forecast(h=2) %>%
  autoplot

######## 3. ARMA SIMULATIONS ###############
n_3 = 300
simulations_number_3 = 100
phi_2_list = c(-0.52, -0.98)
sigma_squared_list = c(0.1^2, 5^2)
parameters_3 = c("phi_1" = 1.5, "phi_2" = NA)
results_3 <- array(NA, dim=c(length(phi_2_list), length(sigma_squared_list),
                             n_3, simulations_number_3))
# dim(results_3): (2, 2, 300, 100) 
# phi, sigma, n_3, simulations_number_3

# simulations
for(i in 1:length(phi_2_list)){
  for(j in 1:length(sigma_squared_list)){
    results_tmp = matrix(NA, nrow = n_3, ncol = simulations_number_3)
    phi_2 = phi_2_list[i]
    sigma_2 = sigma_squared_list[j] 
    for (index in 1:simulations_number_3){
      results_tmp[,index] = arima.sim(n = n_3,
                                  list(ar = c(parameters_3[1], phi_2),
                                  sd = sigma_2))
    }
    results_3[i, j, , ] = as.array(results_tmp)
  }
}

# estimation_phi_2 <- array(NA, dim=c(4, simulations_number_3))
estimation_full <- array(NA, dim=c(2, 2, 2, simulations_number_3))
# phi, sigma, n_3, estimation_parameters, simulations_number_3

for(i in 1:length(phi_2_list)){
  for(j in 1:length(phi_2_list)){
    for(simulation_number in 1:simulations_number_3){
      series = results_3[i,j, ,simulations_number]
      model_tmp_full <- arima(series, order = c(2, 0, 0),
                              include.mean = FALSE)
      estimation_full[i,j, ,simulation_number] =  model_tmp_full$coef[c("ar1", "ar2")]
    }
  }
}
  
    
    
    
for(i in 1:dim(results_3)[1]){
  for(j in 1:dim(results_3)[3]){
    series = results_3[i, ,j]
    # model_tmp = arima(series, order = c(2, 0, 0),
    #                   include.mean = FALSE,
    #                   fixed = parameters_3)
    # estimation_phi_2[i, j] = model_tmp$coef["ar2"]
    
    model_tmp_full <- arima(series, order = c(2, 0, 0),
                            include.mean = FALSE)
    estimation_full[i, ,j] =  model_tmp_full$coef[c("ar1", "ar2")]
  }
}

n_bins = 20
par(mfrow=c(2,2))
hist(t(estimation_phi_2)[, 1], breaks = n_bins)
hist(t(estimation_phi_2)[, 2], breaks = n_bins)
hist(t(estimation_phi_2)[, 3], breaks = n_bins)
hist(t(estimation_phi_2)[, 4], breaks = n_bins)

par(mfrow=c(2,2))
plot(estimation_full[1, 1, ], estimation_full[1, 2, ])
plot(estimation_full[2, 1, ], estimation_full[2, 2, ])
plot(estimation_full[3, 1, ], estimation_full[3, 2, ])
plot(estimation_full[4, 1, ], estimation_full[4, 2, ])

# Conclusionis
# variance doesn't influence the stability
# 






