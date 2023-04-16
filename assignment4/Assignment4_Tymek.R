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

#reading data #####
df <- data.frame(read.table("A4_Kulhuse.csv", sep=",",header=TRUE))
head(df)
dim(df)
df <- df %>% 
  mutate(DateTime = as.POSIXlt(DateTime))

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
# SVR to detect outliers in a signal

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
par(mfrow = c(1, 1))
plot(Y)
sd <- sqrt(k_pure$Sigma.xx.pred[1,1,])
matlines(k_pure$pred[,1] + qnorm(0.975)*cbind(0,sd,-sd),
         col= "red",lty=c(1,2,2))

plot(y_no_outliers)
matlines(k_pure$pred[,1] + qnorm(0.975)*cbind(0,sd,-sd),
         col= "red",lty=c(1,2,2))

residuals = Y - k_pure$pred[1:length(Y), 1]
residuals_sd = (residuals - mean(residuals, na.rm = T))/
  k_pure$Sigma.xx.pred[1, 1,1:length(Y)]
plot(residuals_sd)

## zooming in ####
begining_zoom = 800
end_zoom = 950
range_zoom = seq(from = begining_zoom, to = end_zoom)

plot(Y[range_zoom])
matlines(k_pure$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd[range_zoom],-sd[range_zoom]),
         col= "red",lty=c(1,2,2))

plot(y_no_outliers[range_zoom])
matlines(k_pure$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd[range_zoom],-sd[range_zoom]),
         col= "red",lty=c(1,2,2))

plot(range_zoom, residuals_sd[range_zoom])

# Report the values that defines the 
# final state of the filter (at observation 5000).
# RECONSTRUCTION
k_pure$rec[length(Y),1]
# 20.75815

# ## no outliers primitive ##########
# k2 <- kalman(Y=y_no_outliers, A=A, C=C, Sigma.1=Sigma1, Sigma.2=Sigma2,
#              Xhat0=X0, verbose=TRUE)
# 
# par(mfrow = c(1, 1))
# plot(y_no_outliers)
# sd <- sqrt(k2$Sigma.xx.rec[1,1,])
# matlines(k2$rec[,1] + qnorm(0.975)*cbind(0,c(0,sd), c(0,-sd)),
#          col= "red",lty=c(1,2,2))


# Question 4.4: Skipping outliers when filtering # -----------------------------
k_skipping_outliers <- kalman(Y=Y, A=A, C=C, Sigma.1=Sigma1, Sigma.2=Sigma2,
             Xhat0=X0, verbose=TRUE, rule_6sd = TRUE)

## Plotting predictions ####
par(mfrow = c(1, 1))
plot(Y)
sd <- sqrt(k_skipping_outliers$Sigma.xx.pred[1,1,])
matlines(k_skipping_outliers$pred[,1] + qnorm(0.975)*cbind(0,sd,-sd),
         col= "red",lty=c(1,2,2))

## Zooming in #####
plot(range_zoom, Y[range_zoom])
matlines(range_zoom,
         k_skipping_outliers$pred[range_zoom,1] + qnorm(0.975)*
           cbind(0,sd[range_zoom],-sd[range_zoom]),
         col= "red",lty=c(1,2,2))

# plot(y_no_outliers[range_zoom])
# matlines(k_skipping_outliers$pred[range_zoom,1] + qnorm(0.975)*
#            cbind(0,sd[range_zoom],-sd[range_zoom]),
#          col= "red",lty=c(1,2,2))

sum(!k_skipping_outliers$within_range_vector, na.rm = T)
# 10 observetaions were classified as outliers
# indexes of classified observations
seq(length(Y))[!k_skipping_outliers$within_range_vector] %>%
  na.omit() %>% as.numeric()

# RECONSTRUCTION
k_skipping_outliers$rec[length(Y),1]
# 20.75815

# Question 4.5: Optimizing the variances # -------------------------------------
plot(seq(end_zoom), Y[seq(end_zoom)])
abline(v = begining_zoom)
# isn't the variance increasing?

source("kalman_Tymek.R")
library(numDeriv)

A = matrix(1)
C = matrix(1)
Y = df$Sal
X0 = Y[1]

begining_optimizatioin = 1
end_optimization = 800
# end_optimization = length(Y)

range_optimization = seq(begining_optimizatioin, end_optimization)
Ypart = Y[range_optimization]

objective <- function(theta){
  
  k_optimization <- kalman(Y=Ypart, A=A, C=C,
                           Sigma.1 = matrix(exp(theta["sigma2_system"])),
                           Sigma.2 = matrix(exp(theta["sigma2_observation"])),
               Xhat0=X0, verbose=TRUE, rule_6sd = TRUE)
  
  nepso <- (Ypart[-1] - k_optimization$pred[-c(1, length(Ypart) + 1),1])^2 /
    k_optimization$Sigma.yy.pred[1,1,-c(1, length(Ypart) + 1)]
  
  # return the negative log likelihood
  return(0.5 * sum(nepso + log(k_optimization$Sigma.yy.pred[1,1,-c(1, length(Ypart) + 1)]),
                   na.rm = TRUE))
}


sigma2_observation_empirical = var(Ypart, na.rm = T)
theta_initial = c("sigma2_system" = log(0.01),
                  "sigma2_observation" = log(sigma2_observation_empirical))
# sanity check
objective(theta_initial)

opt = nlminb(theta_initial, objective)
opt
round(exp(opt$par), 6)

Sigma1_ml = matrix(exp(opt$par)["sigma2_system"])
Sigma2_ml = matrix(exp(opt$par)["sigma2_observation"])

k_final = kalman(Y=Y, A=A, C=C,
                 Sigma.1 = Sigma1_ml,
                 Sigma.2 = Sigma2_ml,
                 Xhat0=X0, verbose=TRUE, rule_6sd = TRUE)

## zooming in ####
# range_zoom = seq(length(Y))
plot(range_zoom, Y[range_zoom], ylim = c(17, 19))
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
legend("topright",
       legend=c(
         "Kalman ML optimized", "95 CI",
         "Kalman skipping outliers", "95 CI"),
       col=c("red","red", "blue", "blue"), lty=c(1,2, 1, 2))


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



