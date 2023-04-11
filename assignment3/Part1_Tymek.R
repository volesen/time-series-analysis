# File Description -------------------------------------------------------------
#
#   Tymoteusz Barcinski - s221937
#   
#   Time Series Analysis
#   Assignment 3 - Part 1
#
#_______________________________________________________________________________
rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
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

#reading data #####
df <- data.frame(read.table("A3Data.csv", sep=",",header=TRUE))
head(df)
dim(df)
df <- df %>% 
  mutate(year = as.numeric(substr(X, 1, 4)),
         quarter = as.numeric(substr(X, 6, 6)),
         index = seq(1, dim(df)[1])) %>% 
  mutate(time = as.numeric(year + (quarter - 1)/4)) 
sum(is.na(df))
df_train <- df[1:122, ]

#plotting stuff ##########

par(mfrow = c(3, 4))
plot(df_train$time, df_train$Denmark)
acf_1 <- acf(df_train$Denmark, plot=F)
pacf_1 <- pacf(df_train$Denmark, plot=F)
plot(acf_1)
plot(pacf_1)
hist(df_train$Denmark, breaks = 20)

plot(df_train$time[-1], diff(df_train$Denmark, differences = 1))
acf_2 <- acf(diff(df_train$Denmark, differences = 1), plot=F)
pacf_2 <- pacf(diff(df_train$Denmark, differences = 1), plot=F)
plot(acf_2)
plot(pacf_2)
hist(diff(df_train$Denmark, differences = 1), breaks = 20)

plot(df_train$time[3:length(df_train$time)], diff(df_train$Denmark, differences = 2))
acf_2 <- acf(diff(df_train$Denmark, differences = 2), plot=F)
pacf_2 <- pacf(diff(df_train$Denmark, differences = 2), plot=F)
plot(acf_2)
plot(pacf_2)
hist(diff(df_train$Denmark, differences = 2), breaks = 20)

## transformations #######
par(mfrow = c(1, 1))
boxcox(df_train$Denmark ~ 1)

par(mfrow = c(3, 4))
plot(df_train$time[-1], diff(df_train$Denmark, differences = 1))
acf_2 <- acf(diff(df_train$Denmark, differences = 1), plot=F)
pacf_2 <- pacf(diff(df_train$Denmark, differences = 1), plot=F)
plot(acf_2)
plot(pacf_2)
hist(diff(df_train$Denmark, differences = 1), breaks = 20)

differenced_log = diff(log(df_train$Denmark), differences = 1)
plot((df_train$time[-1]), differenced_log)
acf_2 <- acf(differenced_log, plot=F)
pacf_2 <- pacf(differenced_log, plot=F)
plot(acf_2)
plot(pacf_2)
hist(differenced_log, breaks = 20)

differenced_sqrt = diff(sqrt(df_train$Denmark), differences = 1)
plot((df_train$time[-1]), differenced_sqrt)
acf_2 <- acf(differenced_sqrt, plot=F)
pacf_2 <- pacf(differenced_sqrt, plot=F)
plot(acf_2)
plot(pacf_2)
hist(differenced_sqrt, breaks = 20)

# Modelling ##########

results <- data.frame(matrix(ncol = 3, nrow = 0))
colnames_results <- c("p", "d", "q", "P", "D", "Q", "Frequency",
                      "AIC", "BIC")
results <- data.frame(matrix(ncol = length(colnames_results), nrow = 0))
colnames(results) <- colnames_results

add_model_resulst <- function(model){
  tmp <- c(as.numeric(arimaorder(model)), AIC(model), BIC(model))
  tmp_df <- data.frame(t(tmp))
  colnames(tmp_df) <- colnames(results)
  results <<- rbind(results, tmp_df)
  print(results)
}
  
attach(df_train)

model_auto <- auto.arima(diff(log(Denmark), difference = 1))
summary(model_auto)
# add_model_resulst(model_auto)

m1 <- arima(log(Denmark), order = c(0, 1, 0),
            seasonal = list(order = c(0, 1, 0), period = 4))
summary(m1)
# checkresiduals(m1)
# ggtsdisplay(log(Denmark))
tsdisplay(m1$residuals)
add_model_resulst(m1)

m2 <- arima(log(Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(0, 1, 0), period = 4))
summary(m2)
tsdisplay(m2$residuals)
add_model_resulst(m2)

m3 <- Arima(log(Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m3)
tsdisplay(m3$residuals)
add_model_resulst(m3)

m3_ar_seasonal <- Arima(log(Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(1, 1, 1), period = 4))
summary(m3_ar_seasonal)
tsdisplay(m3_ar_seasonal$residuals)
add_model_resulst(m3_ar_seasonal)
lrt(m3, m3_ar_seasonal)

anna_model <- Arima(log(Denmark), order = c(2, 1, 0),
                        seasonal = list(order = c(1, 1, 1), period = 4))
summary(anna_model)
tsdisplay(anna_model$residuals)
add_model_resulst(anna_model)
lrt(m3, anna_model)


m4 <- arima(log(Denmark), order = c(1, 1, 1),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m4)
tsdisplay(m4$residuals)
add_model_resulst(m4)

auto.arima(m4$residuals)

m5 <- arima(log(Denmark), order = c(1, 1, 2),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m5)
tsdisplay(m5$residuals)
add_model_resulst(m5)

m6 <- arima(log(Denmark), order = c(2, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m6)
c("AIC" = AIC(m6), "BIC" = BIC(m6))
tsdisplay(m6$residuals)
add_model_resulst(m6)

m7 <- arima(log(Denmark), order = c(1, 0, 1),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m7)
# c("AIC" = AIC(m6), "BIC" = BIC(m6))
tsdisplay(m7$residuals)
add_model_resulst(m7)

lrt(m3, m7)
lrt(m7, m3)

#### investing box cox

m_boxcox <- Arima(log(Denmark), order = c(2, 1, 0),
                  seasonal = list(order = c(0, 1, 1), period = 4),
                  lambda = "auto")
summary(m_boxcox)

# m7 <- Arima(log(Denmark), order = c(2, 1, 1),
#             seasonal = list(order = c(0, 1, 1), period = 4),
#             xreg = )
# summary(m7)
# c("AIC" = AIC(m7), "BIC" = BIC(m7))
# tsdisplay(m7$residuals)
# add_model_resulst(m7)

#### paramters: estimated AR, MA + intercept (ignore when differencing)
### + noise_sigma


########### model reduction ############
lrt <- function(model_small, model_big){
  test_stat = abs(as.numeric(-2*(logLik(model_small) - logLik(model_big))))
  df_difference = abs(length(model_big$coef) - length(model_small$coef))
  return(1 - pchisq(test_stat, df = df_difference))
}

lrt(m2, m3)
lrt(m3, m4)
lrt(m4, m5)
# condlucion: we should stay with the model m3
summary(m3)

par(mfrow = c(1, 2))
plot(Denmark)
lines(exp(m3$fitted), col = "red")

m3_test <- Arima(log(Denmark[1:100]), order = c(1, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m3_test)
tsdisplay(m3_test$residuals)
predictions <- predict(m3_test, n.ahead = 10)

par(mfrow = c(1, 2))
plot(Denmark)
lines(exp(m3$fitted), col = "red")

plot(Denmark)
lines(c(exp(m3_test$fitted), exp(predictions$pred)), col = "red")


############## Additional regressors ####################
par(mfrow = c(2, 2))
# plot(InflationRate)
plot(diff(InflationRate, differences = 1))
plot(log(InflationRate))
plot(diff(sqrt(InflationRate), differences = 1))
plot(diff(log(InflationRate), differences = 1))

# boxcox(InflationRate~1)
# model <- auto.arima(sqrt(InflationRate))
# par(mfrow = c(1, 1))
# plot(InflationRate)
# lines(model$fitted^2, col = "red")

ccf(InflationRate, Denmark) # WRONG
# I don't really check that
# the direct cross-correlation function between the input and response
# series gives a misleading indication of the relation between the input
# and response series.

#### WRONG !!!
# m3_xreg <- Arima(log(Denmark), order = c(1, 1, 0),
#             seasonal = list(order = c(0, 1, 1), period = 4),
#             xreg = cbind(df_train$InflationRate, df_train$InterestRate))
# summary(m3_xreg)
# # 
# # lrt(m3, m3_xreg)
# # tsdisplay(m3$residuals)

model_inflation <- auto.arima(diff(InflationRate, differences = 1))
summary(model_inflation)
checkresiduals(model_inflation)
tsdisplay(model_inflation)

arima_reg <- Arima(InflationRate, order = c(3, 1, 1),
                   seasonal = list(order = c(0, 0, 0), period = 0))
summary(arima_reg)
# checkresiduals(m1)
# ggtsdisplay(log(Denmark))
tsdisplay(arima_reg$residuals)

############### prewhitning
pwy_Denmark <- Arima(Denmark, model=arima_reg)
pwx <- arima_reg$residuals

ccf(pwx, pwy_Denmark$residuals, na.action=na.omit)
acf(cbind(pwx, pwy_Denmark$residuals), na.action=na.omit)

xreg_df = cbind(lag9x = stats::lag(InflationRate,-9))

m3_xreg <- Arima(log(Denmark), order = c(1, 1, 0),
                 seasonal = list(order = c(0, 1, 1), period = 4),
                 xreg = xreg_df)
summary(m3_xreg)
 
lrt(m3, m3_xreg)
# tsdisplay(m3$residuals)

#### Interest Rates 
model_interest_rates <- auto.arima(diff(InterestRate, differences = 1))
summary(model_interest_rates)
checkresiduals(model_interest_rates)
# tsdisplay(model_interest_rates$residuals)

# auto.arima(InterestRate)

arima_reg_interest_rates <- Arima(InterestRate, order = c(2, 1, 2),
                   seasonal = list(order = c(0, 0, 0), period = 0))
summary(arima_reg_interest_rates)
checkresiduals(arima_reg_interest_rates)

pwy_Denmark_interest_rates <- Arima(Denmark, model=arima_reg_interest_rates)
pwx <- arima_reg_interest_rates$residuals

ccf(pwx, pwy_Denmark$residuals, na.action=na.omit)
acf(cbind(pwx, pwy_Denmark_interest_rates$residuals), na.action=na.omit)

par(mfrow = c(2, 2))
# plot(InflationRate)
plot(diff(InterestRate, differences = 1))
plot(InterestRate)
plot(diff(sqrt(InterestRate), differences = 1))
plot(diff(log(InterestRate), differences = 1))


#################
# Y <- cumsum(y)
# Y12 <- diffinv(y,12)[-(1:12)]
# Y_1_12 <-cumsum(diffinv(y,12)[-(1:12)])
# (test1 <- arima(y,order=c(1,0,0),xreg=cbind(x1,x2)))
# (test2 <- arima(Y,order=c(1,1,0),xreg=cbind(cumsum(x1),cumsum(x2))))
# (test3 <- arima(Y12,order=c(1,0,0),seasonal=list(order=c(0,1,0),period=12),
#                 xreg=cbind(diffinv(x1,12)[-(1:12)],diffinv(x2,12)[-(1:12)])))
# 
# (testWRONG <- arima(Y,order=c(1,1,0),xreg=cbind(x1,x2)))
# (test4 <- arima(Y_1_12,order=c(1,1,0),seasonal=list(order=c(0,1,0),period=12),
#                 xreg=cbind(cumsum(diffinv(x1,12)[-(1:12)]),cumsum(diffinv(x2,12)[-(1:12)]))))
# 
# 
# s <- 1:10
# d <- diff(s)
# diffinv(d, 5)
# diffinv(s)
# cumsum(s)[-(1:5)]

