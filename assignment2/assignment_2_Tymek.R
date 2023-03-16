rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# options(scipen=999)
options(scipen=0)
# dev.off()
Sys.setenv(LANG = "en")

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

ARMAacf(ar = parameters_1[1], ma =parameters_1[2:3] , lag.max = 10)

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


parameters_2 = c("phi_1" = 1.04, "phi_2" = -0.2, "Phi_1" = 0.86,
               "Theta_1" = -0.42, "mu" = mu)
model <- arima(df$Sales, order = c(2, 0, 0),
               seasonal = list(order = c(1, 0, 1), period = period),
               include.mean = TRUE,
               fixed = parameters_2)
predictions <- predict(model, n.ahead = 2, se.fit = TRUE)

parameters_sarima <- c()
ARMAtoMA(model)

model$sigma2

model_our_sigma <- model
model_our_sigma$sigma2 <- noise_variance
str(model_our_sigma)

predictions_our_sigma <- predict(model_our_sigma, n.ahead = 2, se.fit = TRUE)

stuff <- as.numeric(predictions_our_sigma$pred)
std <- as.numeric(predictions_our_sigma$se)
c(stuff[1] - qnorm(0.975)*std[1], stuff[1] + qnorm(0.975)*std[1])
c(stuff[2] - qnorm(0.975)*std[2], stuff[2] + qnorm(0.975)*std[2])


# ? how to set sigma in the arima function?

# using forecast package
df$Sales %>%
  Arima(order = c(2, 0, 0),
        seasonal = list(order = c(1, 0, 1), period = period),
        include.mean = TRUE,
        fixed = parameters_2) %>%
  forecast(h=2) %>%
  autoplot

######## 3. ARMA SIMULATIONS #####################################3
set.seed(400)

n_3 = 300
simulations_number_3 = 100
# phi_2_list = -c(seq(0.52, 0.98, 0.1), 0.98)
# sigma_list =c(0.1, 1, 3, 5)

phi_2_list = c(-0.52, -0.98)
# phi_2_list = -seq(0.52, 0.98, 0.01)
sigma_list = c(0.1, 5)

parameters_3 = c("phi_1" = 1.5, "phi_2" = NA)
results_3 <- array(NA, dim=c(length(phi_2_list), length(sigma_list),
                             n_3, simulations_number_3))
# dim(results_3): (2, 2, 300, 100) 
# phi, sigma, n_3, simulations_number_3

# simulations
for(i in 1:length(phi_2_list)){
  for(j in 1:length(sigma_list)){
    results_tmp = matrix(NA, nrow = n_3, ncol = simulations_number_3)
    phi_2 = phi_2_list[i]
    sigma = sigma_list[j] 
    for (index in 1:simulations_number_3){
      results_tmp[,index] = arima.sim(n = n_3,
                                  list(ar = c(parameters_3[1], phi_2)),
                                  sd = sigma, n.start=2)
    }
    results_3[i, j, , ] = as.array(results_tmp)
  }
}

estimation_full <- array(NA, dim=c(length(phi_2_list), length(sigma_list), 2, simulations_number_3), 
                         dimnames = list(phi_2_list, sigma_list,
                                         c("phi_1_estimated", "phi_2_estimated")))
# phi, sigma, estimation_parameters, simulations_number_3

for(i in 1:length(phi_2_list)){
  for(j in 1:length(sigma_list)){
    print(c(phi_2_list[i], sigma_list[j]))
    for(simulation_number in 1:simulations_number_3){
      series = results_3[i,j, ,simulation_number]
      model_tmp_full <- arima(series, order = c(2, 0, 0),
                              include.mean = FALSE)
      estimation_full[i,j, ,simulation_number] =  model_tmp_full$coef[c("ar1", "ar2")]
    }
  }
}

estimation_full = -estimation_full

par(mfrow = c(1, 1))
plot(estimation_full[1, 1,1,], estimation_full[1, 1,2,], col = "red",
     xlim = c(-1.7, -1.3), ylim = c(0.4, 1.2),
     xlab="phi_1", ylab = "phi_2",
     main = "Estimates over simulations along with the stable region")
lines(estimation_full[1, 2,1, ], estimation_full[1, 2,2, ], type = "p", col = "green")
lines(estimation_full[2, 1,1, ], estimation_full[2, 1,2, ], type = "p", col = "black")
lines(estimation_full[2, 2,1, ], estimation_full[2, 2,2, ], type = "p", col = "blue")
polygon(x = c(-2, 2, 0),                           # X-Coordinates of polygon
        y = c(1, 1, -1),                             # Y-Coordinates of polygon
        col = rgb(0, 0, 0.5, alpha = 0.1))
abline(v = -1.5, lty = 2)
abline(h = -phi_2_list,lty = 2)
legend("topright", 95,
       legend=c("phi_2 = 0.52, sigma=0.1",
                "phi_2 = 0.52, sigma=5",
                "phi_2 = 0.98, sigma=0.1",
                "phi_2 = 0.98, sigma=5",
                "Stable region"),
       col=c("red", "green", "black", "blue", rgb(0, 0, 0.5, alpha = 0.1)),
       lty = c(NA,NA,NA,NA, 1), cex=0.8,
       inset = 0.02, pch=c(1,1,1,1, NA))

results <- matrix(NA, nrow = 4, ncol = 5)
par(mfrow=c(length(sigma_list), length(phi_2_list)))

phi_2_small_xlim = c(0.35, 0.65)
phi_2_large_xlim =  c(0.9, 1)

i = 1
for(sigma_2_index in 1:length(sigma_list)){
  for(phi_2_index in 1:length(phi_2_list)){
    
    sd_res <- sd(results_3[phi_2_index, sigma_2_index, , ])
    ifelse(phi_2_list[phi_2_index] == -0.52,
           hist(estimation_full[phi_2_index, sigma_2_index, 2, ], breaks = 25, xlab="phi_2_estimated",
                main = paste("phi_2=", -phi_2_list[phi_2_index],
                             " | sigma=",  sigma_list[sigma_2_index],
                             " | sd_phi_2=", round(sd(estimation_full[phi_2_index, sigma_2_index, 2, ]), 3)),
                
                ylim = c(0, 20), xlim = phi_2_small_xlim),
           hist(estimation_full[phi_2_index, sigma_2_index, 2, ], breaks = 25, xlab="phi_2_estimated",
                main = paste("phi_2=", -phi_2_list[phi_2_index],
                             " | sigma=",  sigma_list[sigma_2_index],
                             " | sd_phi_2=", round(sd(estimation_full[phi_2_index, sigma_2_index, 2, ]), 3)),
                
                ylim = c(0, 20), xlim = phi_2_large_xlim)
           )

    phi_2_estimated_mean = mean(estimation_full[phi_2_index, sigma_2_index, 2, ])
    phi_2_estimated_sd = sd(estimation_full[phi_2_index, sigma_2_index, 2, ])
    CI_lower = phi_2_estimated_mean - qnorm(0.975)*phi_2_estimated_sd/sqrt(simulations_number_3)
    CI_upper = phi_2_estimated_mean + qnorm(0.975)*phi_2_estimated_sd/sqrt(simulations_number_3)
    quantile_empricial = quantile(estimation_full[phi_2_index, sigma_2_index, 2, ], 0.95)
    
    abline(v = phi_2_estimated_mean, col = "red", lty = 1, lwd = 2)
    abline(v = CI_lower, lwd = 2, lty = 2, col = "blue")
    abline(v = CI_upper, lwd = 2, lty = 2, col = "blue")
    abline(v = quantile_empricial, lwd = 2, lty = 2, col = "black")
    
    results[i, ] <- c(phi_2_estimated_mean, phi_2_estimated_sd, CI_lower, CI_upper, quantile_empricial)
    i = i + 1
  }
}
round(results, 4)

plot(phi_2_list_plotting, df_sd_1,
     xlab="phi_2 true value", ylab = "standard deviation of phi_2 estimates",
     main = "Standard deviation of the phi_2 estimates as a function of true value of phi_2")
lines(phi_2_list_plotting, df_sd_2, col = "red", type = "p")
lines(phi_2_list_plotting, fitted.values(model_1), col = "black")
lines(phi_2_list_plotting, fitted.values(model_2), col = "red")
legend("topright", 95,
       legend=c("sigma=0.1", "polynomial regression; degree=2",
                "sigma=5", "polynomial regression; degree=2"),
       col=c("black", "black", "red", "red"), lty = c(NA,1,NA,1), cex=0.8,
       inset = 0.02, pch=c(1,NA,1,NA))


par(mfrow=c(length(sigma_list), length(phi_2_list)))
# par(mfrow=c(1, 1))
for(sigma_2_index in 1:length(sigma_list)){
  for(phi_2_index in 1:length(phi_2_list)){
    sd_res <- sd(results_3[phi_2_index, sigma_2_index, , ])
    
    matplot(results_3[phi_2_index, sigma_2_index, , 1:4], type = "l",
          lty = 1:1,
         main = paste("phi_2=", -phi_2_list[phi_2_index],
                      " | sigma=",  sigma_list[sigma_2_index],
                      " | sd_phi_2=",
                      round(sd(estimation_full[phi_2_index, sigma_2_index, 2, ]), 3)),
         xlab="", ylab = "X_t",
         
    )
  }
}

par(mfrow=c(length(sigma_list), length(phi_2_list)))
for(sigma_2_index in 1:length(sigma_list)){
  for(phi_2_index in 1:length(phi_2_list)){
    sd_res <- sd(results_3[phi_2_index, sigma_2_index, , ])
    
    plot(estimation_full[phi_2_index, sigma_2_index, 1, ],
         estimation_full[phi_2_index, sigma_2_index, 2, ],
         xlab="phi_1_estimated", ylab = "phi_2_estimated",
         xlim = c(-1.6, -1.4), 
         main = paste("phi_2=", -phi_2_list[phi_2_index],
                      " | sigma=",  sigma_list[sigma_2_index],
                      " | sd_phi_2=",
                      round(sd(estimation_full[phi_2_index, sigma_2_index, 2, ]), 3)))
  }
}
