# File Description -------------------------------------------------------------
#
#   Tymoteusz Barcinski - s221937
#   
#   Time Series Analysis
#   Assignment 4 - Kalman Filter
#
#_______________________________________________________________________________
rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
options(scipen=0)
Sys.setenv(LANG = "en")

library(MASS)
library(dplyr)
# library(tsibble)
library(forecast)
# library(matlib)
library(nlme)
library(urca)
library(stats)
library(lmtest)

# reading data #####
df <- data.frame(read.table("A4_Kulhuse.csv", sep=",",header=TRUE))
head(df)
dim(df)
df <- df %>% 
  mutate(DateTime = as.POSIXlt(DateTime))

# time conversion
time <- as.POSIXct(dplyr::pull(df, DateTime))

# NA values
apply(is.na(df), 2, sum)
seq(length(df$Sal))[is.na(df$Sal)]

# Question 4.1: Presenting data #----------------------------------------------
plot(df$Sal, type = "p")
# there are a 4 outliers
outliers_sal = seq(length(df$Sal))[df$Sal < mean(df$Sal, na.rm = T) -
                      3*sd(df$Sal, na.rm = T)] %>% na.omit() %>% as.numeric()
outliers_sal

y_no_outliers = df$Sal
y_no_outliers[outliers_sal] = NA
plot(y_no_outliers)

plot(df$ODO)

plot(df$Sal, ylim = c(16.5, 21.5), pch = 20,
     ylab = "Salinity [PSU]",
     xlab = "Time [30 min.]")
# SVR to detect outliers in a signal

plot(df$ODO, ylim = c(7, 11.5), pch = 20,
     ylab = "Disolved oxygen [mg/L]",
     xlab = "Time [30 min.]")


# Question 4.2: Random walk state-space model of salinity #---------------------
# Xt = X_{t-1} + et
# Yt = Xt + epsilon_t

## ARIMA stuff #----------------------------------------------------------------
m1 <- Arima(y_no_outliers, order = c(0, 1, 1))
summary(m1)
tsdisplay(m1$residuals)

# m_auto = auto.arima(y_no_outliers)
# summary(m_auto)

# Question 4.3: Pure Kalman filter # -------------------------------------------
source("kalman_Tymek.R")

A = matrix(1)
C = matrix(1)

Y = df$Sal
X0 = Y[1]
Sigma1 = matrix(0.01)
Sigma2 = matrix(0.005)

k_pure <- kalman(Y=Y, A=A, C=C, Sigma.1=Sigma1, Sigma.2=Sigma2,
             Xhat0=X0, verbose=TRUE)
names(k_pure)

# # plotting reconstructions
# par(mfrow = c(1, 1))
# sd <- sqrt(k1$Sigma.xx.rec[1,1,])
# matlines(k1$rec[,1] + qnorm(0.975)*cbind(0,c(0,sd), c(0,-sd)),
#          col= "red",lty=c(1,2,2))

## plotting predictions all ######
residuals = Y - k_pure$pred[1:length(Y), 1]
residuals_sd = (residuals - mean(residuals, na.rm = T))/
  sqrt(k_pure$Sigma.yy.pred[1, 1,1:length(Y)])

par(mfrow = c(1, 2))
plot(Y, pch = 20,
     ylab = "Salinity [PSU]")
lines(outliers_sal, Y[outliers_sal], type = "p", col = "blue",
      lwd = 2, pch = 20)
lines(outliers_sal, Y[outliers_sal], type = "p", col = "blue",
      pch = 20)
sd <- sqrt(k_pure$Sigma.yy.pred[1,1,])
matlines(k_pure$pred[,1] + qnorm(0.975)*cbind(0,sd,-sd),
         col= c("green", "red", "red"),lty=c(1,2,2))
# par(mfrow = c(1, 1))
plot(residuals_sd, ylab = "Standarized residuals", pch = 20)
lines(outliers_sal, residuals_sd[outliers_sal], type = "p", col = "blue",
      lwd = 2, pch = 20)

par(mfrow = c(1, 1))
plot(Y, ylab = "Salinity [PSU]", ylim = c(16.5, 21.5), pch = 20)
lines(outliers_sal, Y[outliers_sal], type = "p", col = "blue",
      lwd = 2)
lines(outliers_sal, Y[outliers_sal], type = "p", col = "blue")
sd <- sqrt(k_pure$Sigma.xx.pred[1,1,])
matlines(k_pure$pred[,1] + qnorm(0.975)*cbind(0,sd,-sd),
         col= c("green", "red", "red"),lty=c(1,2,2))

# plot(y_no_outliers)
# matlines(k_pure$pred[,1] + qnorm(0.975)*cbind(0,sd,-sd),
#          col= "red",lty=c(1,2,2))
# 
# 
# par(mfrow = c(1, 2))
# Acf(residuals_sd, na.action = na.pass)
# Pacf(residuals_sd, na.action = na.pass)


## zooming in ####
begining_zoom = 800
end_zoom = 950
range_zoom = seq(from = begining_zoom, to = end_zoom)
outliers_sal_zoom = outliers_sal[outliers_sal %in% range_zoom]

par(mfrow = c(1, 2))
plot(range_zoom, Y[range_zoom],
     ylab = "Salinity [PSU]",
     xlab = "Index", pch = 20)
lines(outliers_sal_zoom, Y[outliers_sal_zoom],
      type = "p", col = "blue", lwd = 2)

sd <- sqrt(k_pure$Sigma.yy.pred[1,1,])
matlines(range_zoom,
         k_pure$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd[range_zoom],-sd[range_zoom]),
         col= c("green", "red", "red"), lty=c(1,2,2),
         pch = 20)

plot(range_zoom, residuals_sd[range_zoom],
     ylab = "Standarized residuals",
     xlab = "Index", pch = 20)
lines(outliers_sal_zoom, residuals_sd[outliers_sal_zoom],
      type = "p", col = "blue", lwd = 2, pch = 20)

# Report the values that defines the 
# final state of the filter (at observation 5000).
# RECONSTRUCTION
k_pure$rec[length(Y),1]
# 20.75815

kalman_gain = k_pure$Sigma.xx.pred / k_pure$Sigma.yy.pred

par(mfrow = c(2, 1))
plot(kalman_gain)
plot(df$Sal)

# Question 4.4: Skipping outliers when filtering # -----------------------------

k_skipping_outliers <- kalman(Y=Y, A=A, C=C, Sigma.1=Sigma1, Sigma.2=Sigma2,
             Xhat0=X0, verbose=TRUE, rule_6sd = TRUE)

residuals = Y - k_skipping_outliers$pred[1:length(Y), 1]
residuals_sd_skipped = (residuals - mean(residuals, na.rm = T))/
  sqrt(k_pure$Sigma.yy.pred[1, 1,1:length(Y)])

sum(!k_skipping_outliers$within_range_vector, na.rm = T)
# 10 observetaions were classified as outliers
# indexes of classified observations
skipped_observations = seq(length(Y))[!k_skipping_outliers$within_range_vector] %>%
  na.omit() %>% as.numeric()
skipped_observations

## Plotting #######

### Zooming in ######
no_outliers = k_skipping_outliers$within_range_vector[range_zoom]
begining_zoom = 800
end_zoom = 950
range_zoom = seq(from = begining_zoom, to = end_zoom)
outliers_sal_zoom = skipped_observations[skipped_observations %in% range_zoom]

par(mfrow = c(1, 2))
plot(range_zoom, Y[range_zoom],
     ylab = "Salinity [PSU]",
     xlab = "Index", pch = 20)

sd <- sqrt(k_skipping_outliers$Sigma.yy.pred[1,1,])
matlines(range_zoom,
         k_skipping_outliers$pred[range_zoom, 1]
         + qnorm(0.975)*
           cbind(0,sd[range_zoom],
                 -sd[range_zoom]),
         col= c("green", "red","red"), lty=c(1,2,2), pch = 20)
abline(v = outliers_sal_zoom, col = "blue")

plot(range_zoom[no_outliers], residuals_sd_skipped[range_zoom][no_outliers],
     ylab = "Standarized residuals",
     xlab = "Index", pch = 20)
abline(h = 0, col = "red")
abline(v = outliers_sal_zoom, col = "blue")

### Full plot ######
par(mfrow = c(1, 2))
plot(Y, ylab = "Salinity [PSU]",
     xlab = "Index")

sd <- sqrt(k_skipping_outliers$Sigma.xx.pred[1,1,])
matlines(seq(length(Y)),
         k_skipping_outliers$pred[seq(length(Y)), 1]
         + qnorm(0.975)*
           cbind(0,sd[seq(length(Y))],
                 -sd[seq(length(Y))]),
         col= "red",lty=c(1,2,2))
abline(v = outliers_sal_zoom, col = "green")

plot(residuals_sd_skipped[range_zoom][no_outliers],
     ylab = "Standarized residuals",
     xlab = "Index")
abline(h = 0, col = "blue")
abline(v = outliers_sal_zoom, col = "green")

residuals_true = residuals_sd_skipped
residuals_true[seq(length(residuals_sd_skipped)) %in%
                 skipped_observations] = NA

plot(residuals_true)

par(mfrow = c(1, 2))
Acf(residuals_true, na.action = na.pass, main = "")
Pacf(residuals_true, na.action = na.pass, main = "")

par(mfrow = c(1, 1))
plot(residuals_true, type = "p")
arima_residuals = Arima(residuals_true, order = c(4, 0, 1))
tsdisplay(arima_residuals$residuals)

## Plotting predictions ####
par(mfrow = c(1, 1))
plot(Y)
sd <- sqrt(k_skipping_outliers$Sigma.xx.pred[1,1,])
matlines(k_skipping_outliers$pred[,1] + qnorm(0.975)*cbind(0,sd,-sd),
         col= "red",lty=c(1,2,2))

investigation = 
  
par(mfrow = c(1, 2))
plot(k_skipping_outliers$pred)
lines(k_skipping_outliers$rec, type = "p", col = "red")
plot(k_skipping_outliers$pred - k_skipping_outliers$rec)

Acf(k_skipping_outliers$pred - k_skipping_outliers$rec,
    na.action = na.pass)

## Zooming in #####
plot(range_zoom, Y[range_zoom])
matlines(range_zoom,
         k_skipping_outliers$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd[range_zoom],-sd[range_zoom]),
         col= "red",lty=c(1,2,2))

# RECONSTRUCTION
k_skipping_outliers$rec[length(Y),1]
# 20.75815

k_pure$

# Question 4.5: Optimizing the variances # -------------------------------------

source("kalman_Tymek.R")
library(numDeriv)

A = matrix(1)
C = matrix(1)
Y = df$Sal
X0 = Y[1]

plot(seq(end_zoom), Y[seq(end_zoom)])
abline(v = begining_zoom)
# isn't the variance increasing?
sd_800 = sd(Y[seq(begining_zoom)])



begining_optimizatioin = 1
end_optimization = 800
# end_optimization = length(Y)

range_optimization = seq(begining_optimizatioin, end_optimization)
Ypart = Y[range_optimization]
sigma2_observation_empirical = var(Ypart, na.rm = T)

theta_initial = c("sigma2_system" = log(0.01),
                  "sigma2_observation" = log(sigma2_observation_empirical))

kalman_test = kalman(Y=Ypart, A=A, C=C,
                     Sigma.1 = matrix(exp(theta_initial["sigma2_system"])),
                     Sigma.2 = matrix(exp(theta_initial["sigma2_observation"])),
                     Xhat0=X0, verbose=TRUE, rule_6sd = TRUE,
                     rule_6sd_optimization_value = sd_800)

kalman_test$sd_vector
plot(sqrt(kalman_test$Sigma.yy.pred[1, 1,]))


objective <- function(theta){
  
  k_optimization <- kalman(Y=Ypart, A=A, C=C,
                           Sigma.1 = matrix(exp(theta["sigma2_system"])),
                           Sigma.2 = matrix(exp(theta["sigma2_observation"])),
               Xhat0=X0, verbose=TRUE, rule_6sd = TRUE,
               rule_6sd_optimization_value = sd_800)
  print(
    c(k_optimization$Sigma.yy.rec[1, 1,dim(k_optimization$Sigma.yy.rec)[3]],
      k_optimization$Sigma.yy.pred[1, 1,dim(k_optimization$Sigma.yy.pred)[3]],
      k_optimization$Sigma.xx.rec[1, 1,dim(k_optimization$Sigma.xx.rec)[3]],
      k_optimization$Sigma.xx.pred[1, 1,dim(k_optimization$Sigma.xx.pred)[3]],
      sum(!k_optimization$within_range_vector))
    )
  
  nepso <- (Ypart[-1] - k_optimization$pred[-c(1, length(Ypart) + 1),1])^2 /
    k_optimization$Sigma.yy.pred[1,1,-c(1, length(Ypart) + 1)]
  
  # return the negative log likelihood
  return(0.5 * sum(nepso + log(k_optimization$Sigma.yy.pred[1,1,-c(1, length(Ypart) + 1)]),
                   na.rm = TRUE))
}



theta_initial = c("sigma2_system" = log(0.01),
                  "sigma2_observation" = log(0.01))
# sanity check
# objective(theta_initial)

opt = nlminb(theta_initial, objective)
opt_2 = optim(theta_initial, objective)

opt

opt$par
exp(opt$par)
round(exp(opt$par), 6)

objective_forced <- function(theta_forced){
  
  k_optimization <- kalman(Y=Ypart, A=A, C=C,
                           Sigma.1 = matrix(exp(theta_forced["sigma2_system"])),
                           Sigma.2 = matrix(0),
                           Xhat0=X0, verbose=TRUE, rule_6sd = TRUE,
                           rule_6sd_optimization_value = sd_800)
  # print(
  #   c(k_optimization$Sigma.yy.rec[1, 1,dim(k_optimization$Sigma.yy.rec)[3]],
  #     k_optimization$Sigma.yy.pred[1, 1,dim(k_optimization$Sigma.yy.pred)[3]],
  #     k_optimization$Sigma.xx.rec[1, 1,dim(k_optimization$Sigma.xx.rec)[3]],
  #     k_optimization$Sigma.xx.pred[1, 1,dim(k_optimization$Sigma.xx.pred)[3]])
  # )
  # 
  nepso <- (Ypart[-1] - k_optimization$pred[-c(1, length(Ypart) + 1),1])^2 /
    k_optimization$Sigma.yy.pred[1,1,-c(1, length(Ypart) + 1)]
  
  # return the negative log likelihood
  return(0.5 * sum(nepso + log(k_optimization$Sigma.yy.pred[1,1,-c(1, length(Ypart) + 1)]),
                   na.rm = TRUE))
}
theta_initial_forced = c("sigma2_system" = log(0.01))

opt_forced = nlminb(theta_initial_forced, objective_forced)
opt_forced$par
exp(opt_forced$par)
round(exp(opt_forced$par), 6)

# the likelihoods are the same - hence we can set variance of sigma 
# observations to 0. Value: -2106.553
opt_forced$objective
opt$objective

Sigma1_ml = matrix(exp(opt_forced$par))
Sigma2_ml = matrix(0)

Sigma1_ml_report = matrix(1.885*10^(-3))
Sigma2_ml_report = matrix(6.765*10^(-6))

Sigma1_ml = Sigma1_ml_report
Sigma2_ml = Sigma2_ml_report
  
k_final = kalman(Y=Y, A=A, C=C,
                 Sigma.1 = Sigma1_ml,
                 Sigma.2 = Sigma2_ml,
                 Xhat0=X0, verbose=TRUE, rule_6sd = TRUE)

residuals_final = Y - k_final$pred[1:length(Y), 1]
residuals_final_skipped = (residuals_final - mean(residuals_final, na.rm = T))/
  sqrt(k_final$Sigma.xx.pred[1, 1,1:length(Y)])

## zooming in ####
# range_zoom = seq(length(Y))

layout(matrix(c(1, 1, 2), nrow=1, byrow=TRUE))

range_zoom = seq(800, 950)
plot(range_zoom, Y[range_zoom], ylim = c(17, 18.5),
     ylab = "Y", xlab = "Index", main="")
grid()
sd_ml <- sqrt(k_final$Sigma.xx.pred[1,1,])
matlines(range_zoom,
         k_final$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd_ml[range_zoom],-sd_ml[range_zoom]),
         col= "red",lty=c(1,2,2))
sd_3 <- sqrt(k_skipping_outliers$Sigma.xx.pred[1,1,])
matlines(range_zoom,
         k_skipping_outliers$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd_3[range_zoom],-sd_3[range_zoom]),
         col= "blue",lty=c(1,2,2))
abline(v = c(887, 894), col = "green")
legend("top",
       legend=c(
         "Kalman ML optimized", "95 PI",
         "Kalman skipping outliers", "95 PI"),
       col=c("red","red", "blue", "blue"), lty=c(1,2, 1, 2))

range_zoom = seq(880, 900)
plot(range_zoom, Y[range_zoom], ylim = c(17, 18.5),
     ylab = "Y", xlab = "Index", main="")
grid()
sd_ml <- sqrt(k_final$Sigma.xx.pred[1,1,])
matlines(range_zoom,
         k_final$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd_ml[range_zoom],-sd_ml[range_zoom]),
         col= "red",lty=c(1,2,2))
sd_3 <- sqrt(k_skipping_outliers$Sigma.xx.pred[1,1,])
matlines(range_zoom,
         k_skipping_outliers$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd_3[range_zoom],-sd_3[range_zoom]),
         col= "blue",lty=c(1,2,2))
abline(v = c(887, 894), col = "green")
legend("top",
       legend=c(
         "Kalman ML optimized", "95 PI",
         "Kalman skipping outliers", "95 PI"),
       col=c("red","red", "blue", "blue"), lty=c(1,2, 1, 2))


residuals_sd_skipped[range_zoom][no_outliers]
residuals_final_skipped

par(mfrow = c(2, 2))
plot(residuals_sd_skipped[range_zoom][no_outliers])
plot(residuals_final_skipped[range_zoom][no_outliers],
      type = "p", col = "red")
plot(residuals[range_zoom][no_outliers])
plot(residuals_final[range_zoom][no_outliers],
     type = "p", col = "red")

## full plots ####
plot(Y)
sd <- sqrt(k_final$Sigma.xx.pred[1,1,])
matlines(k_final$pred[,1] + qnorm(0.975)*
           cbind(0,sd,-sd),
         col= "red",lty=c(1,2,2))

# ## residual analysis ####
# residuals_final_standarized = Y - k_final$pred[seq(end_zoom), 1]
# 
# 
# residuals_skipping_standarized = (residuals - mean(residuals, na.rm = T))/
#   k_skipping_outliers$Sigma.xx.pred[1, 1,seq(end_zoom)]
# plot(residuals_sd)
# lines(residuals_ml)
# 
# sum(residuals_ml^2, na.rm = T)
# # 1257.6
# sum((Y - k3$pred[1:length(Y), 1])^2, na.rm = T)
# # 995.9976
# # worse results bc we estimated the variance of the noise only on the
# # first 800 which is not representative
# 
# ## Results with the optimization based on all observations ####
# sum(residuals_ml^2, na.rm = T)
# # 983.8123
# sd(residuals_ml, na.rm = T)
# # 0.4485359
# 
# residuals_regular = Y - k3$pred[1:length(Y), 1]
# sum((residuals_regular)^2, na.rm = T)
# # 995.9976
# sd(residuals_regular, na.rm = T)
# # 0.4512983

# RECONSTRUCTION
k_final$rec[length(Y),1]
# 20.77

# Question 4.6: Model for dissolved oxygen # ----------------------------------
plot(df$ODOsat)
abline(h = 100, col = "red")
acf(df$ODOsat)



