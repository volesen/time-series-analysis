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
library(urca)
library(stats)
library(lmtest)

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

df[, "InterestRate"][is.na(df[, "InterestRate"])] = df$InterestRate[124]
df[, "InflationRate"][is.na(df[, "InflationRate"])] = df$InflationRate[124]

# plotting stuff ##########
## plots of the differencing #####
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
grid()
acf_2 <- acf(diff(df_train$Denmark, differences = 1), plot=F)
pacf_2 <- pacf(diff(df_train$Denmark, differences = 1), plot=F)
plot(acf_2)
plot(pacf_2)
hist(diff(df_train$Denmark, differences = 1), breaks = 20)

differenced_log = diff(log(df_train$Denmark), differences = 1)
plot((df_train$time[-1]), differenced_log)
grid()
acf_2 <- acf(differenced_log, plot=F)
pacf_2 <- pacf(differenced_log, plot=F)
plot(acf_2)
plot(pacf_2)
hist(differenced_log, breaks = 20)

differenced_log_2 = diff(log(df_train$Denmark), differences = 2)
plot((df_train$time[-(1:2)]), differenced_log_2)
grid()
acf_2 <- acf(differenced_log_2, plot=F)
pacf_2 <- pacf(differenced_log_2, plot=F)
plot(acf_2)
plot(pacf_2)
hist(differenced_log_2, breaks = 20)

# differenced_sqrt = diff(sqrt(df_train$Denmark), differences = 1)
# plot((df_train$time[-1]), differenced_sqrt)
# acf_2 <- acf(differenced_sqrt, plot=F)
# pacf_2 <- pacf(differenced_sqrt, plot=F)
# plot(acf_2)
# plot(pacf_2)
# hist(differenced_sqrt, breaks = 20)

## only seasonal differencing ###########
par(mfrow = c(2, 4))
seasonal_4 = diff(df_train$Denmark, differences = 1, lag = 4)
plot(df_train$time[-(1:4)], seasonal_4)
grid()
acf_2 <- acf(seasonal_4, plot=F)
pacf_2 <- pacf(seasonal_4, plot=F)
plot(acf_2)
plot(pacf_2)
hist(seasonal_4, breaks = 20)

seasonal_4_log = diff(log(df_train$Denmark), differences = 1, lag = 4)
plot(df_train$time[-(1:4)], seasonal_4_log)
grid()
acf_2 <- acf(seasonal_4_log, plot=F)
pacf_2 <- pacf(seasonal_4_log, plot=F)
plot(acf_2)
plot(pacf_2)
hist(seasonal_4_log, breaks = 20)

seasonal_4_sqrt = diff(sqrt(df_train$Denmark), differences = 1, lag = 4)
plot(df_train$time[-(1:4)], seasonal_4_sqrt)
grid()
acf_2 <- acf(seasonal_4_sqrt, plot=F)
pacf_2 <- pacf(seasonal_4_sqrt, plot=F)
plot(acf_2)
plot(pacf_2)
hist(seasonal_4_sqrt, breaks = 20)

test_unit_roots = ur.kpss(seasonal_4_log)
summary(test_unit_roots)
summary(ur.kpss(seasonal_4_sqrt))

par(mfrow = c(1, 2))
plot(Denmark)
acf(Denmark)

## final comparison ###
par(mfrow = c(3, 4))

differenced_log = diff(log(df_train$Denmark), differences = 1)
plot((df_train$time[-1]), differenced_log)
grid()
acf_2 <- acf(differenced_log, plot=F)
pacf_2 <- pacf(differenced_log, plot=F)
plot(acf_2)
plot(pacf_2)
hist(differenced_log, breaks = 20)

seasonal_4_log = diff(log(df_train$Denmark), differences = 1, lag = 4)
plot(df_train$time[-(1:4)], seasonal_4_log)
grid()
acf_2 <- acf(seasonal_4_log, plot=F)
pacf_2 <- pacf(seasonal_4_log, plot=F)
plot(acf_2)
plot(pacf_2)
hist(seasonal_4_log, breaks = 20)

summary(ur.kpss(seasonal_4_log))

m1 <- arima(log(Denmark), order = c(0, 1, 0),
            seasonal = list(order = c(0, 1, 0), period = 4))
summary(m1)
plot(df_train$time, m1$residuals)
grid()
acf_2 <- acf(m1$residuals, plot=F)
pacf_2 <- pacf(m1$residuals, plot=F)
plot(acf_2)
plot(pacf_2)
hist(m1$residuals, breaks = 20)

summary(ur.kpss(m1$residuals))
# decision to do the differencing two times

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

# model_boxcox_auto <- Arima((Denmark), order = c(1, 1, 0),
#                         seasonal = list(order = c(0, 1, 1), period = 4),
#                         lambda = "auto")
# model_boxcox_auto
# lambda= -0.06814993 
# support for the log transformation
# 
# model_auto <- auto.arima(diff(log(Denmark), difference = 1))
# summary(model_auto)
# add_model_resulst(model_auto)

### differencing regular and seasonal #######
m0 <- Arima(log(Denmark), order = c(0, 1, 0),
            seasonal = list(order = c(0, 0, 0), period = 4))
summary(m0)
tsdisplay(m0$residuals)

m0_v2 <- Arima(log(Denmark), order = c(0, 1, 0),
            seasonal = list(order = c(1, 0, 0), period = 4))
summary(m0_v2)
tsdisplay(m0_v2$residuals)

m0_v3 <- Arima(log(Denmark), order = c(1, 1, 0),
               seasonal = list(order = c(1, 0, 1), period = 4))
summary(m0_v3)
tsdisplay(m0_v3$residuals)
summary(ur.kpss(m0_v3$residuals))

par(mfrow = c(1, 1))
qqnorm(m0_v3$residuals)
qqline(m0_v3$residuals)
shapiro.test(m0_v3$residuals)



summary(ur.kpss(m0$residuals))
# data can be considered stationary
m0 <- Arima(log(Denmark), order = c(0, 1, 0),
            seasonal = list(order = c(0, 0, 0), period = 4))
summary(m0)
tsdisplay(m0$residuals)

###################################################################
m1 <- Arima(log(Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(1, 0, 0), period = 4))
summary(m1)
# checkresiduals(m1)
# ggtsdisplay(log(Denmark))
tsdisplay(m1$residuals)
plot(m1$residuals, type = "p")
kpss.test(m1$residuals)



add_model_resulst(m1)

m2 <- Arima(log(Denmark), order = c(0, 1, 1),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m2)
tsdisplay(m2$residuals)
add_model_resulst(m2)

m3 <- Arima(log(Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m3)
tsdisplay(m3$residuals)
add_model_resulst(m3)

library(tseries)
kpss.test(m3$residuals)
plot(m3$residuals, type = "p")

m4 <- Arima(log(Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(1, 1, 1), period = 4))
summary(m4)
tsdisplay(m4$residuals)
add_model_resulst(m4)

m5 <- Arima(log(Denmark), order = c(1, 1, 1),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m5)
tsdisplay(m5$residuals)
add_model_resulst(m5)

m6 <- Arima(log(Denmark), order = c(2, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m6)
c("AIC" = AIC(m6), "BIC" = BIC(m6))
tsdisplay(m6$residuals)
add_model_resulst(m6)

m7 <- Arima(log(Denmark), order = c(1, 1, 1),
            seasonal = list(order = c(1, 1, 1), period = 4))
summary(m7)
# c("AIC" = AIC(m6), "BIC" = BIC(m6))
tsdisplay(m7$residuals)
add_model_resulst(m7)

m10 <- Arima(log(Denmark), order = c(0, 1, 1),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m10)
# c("AIC" = AIC(m6), "BIC" = BIC(m6))
tsdisplay(m10$residuals)
add_model_resulst(m10)

options(scipen=999)
t(round(rbind(
  as.numeric(lrtest(m3, m1)$'Pr(>Chisq)')[2],
  as.numeric(lrtest(m3, m2)$'Pr(>Chisq)')[2],
  as.numeric(lrtest(m3, m4)$'Pr(>Chisq)')[2],
  as.numeric(lrtest(m3, m5)$'Pr(>Chisq)')[2],
  as.numeric(lrtest(m3, m6)$'Pr(>Chisq)')[2],
  as.numeric(lrtest(m3, m7)$'Pr(>Chisq)')[2]
), 4))
options(scipen=0)


## AR vs differencing ###################################
m1log = Arima(log(Denmark), order = c(0, 1, 0),
              seasonal = list(order = c(1, 0, 0), period = 4))
summary(m1log)
tsdisplay(m1log$residuals)
add_model_resulst(m1log)

m2log = Arima(log(Denmark), order = c(1, 1, 0),
              seasonal = list(order = c(1, 0, 0), period = 4))
summary(m2log)
tsdisplay(m2log$residuals)
add_model_resulst(m2log)

# par(mfrow = c(1, 1))
# qqnorm(m2log$residuals)
# qqline(m2log$residuals)
# shapiro.test(m2log$residuals)

m3log = Arima(log(Denmark), order = c(1, 1, 0),
              seasonal = list(order = c(1, 0, 1), period = 4))
summary(m3log)
tsdisplay(m3log$residuals)
add_model_resulst(m3log)


autoplot(m3log)
autoplot(m3)


par(mfrow = c(1, 1))
qqnorm(m2log$residuals)
qqline(m2log$residuals)
shapiro.test(m2log$residuals)

m4log = Arima(log(Denmark), order = c(1, 1, 0),
              seasonal = list(order = c(1, 0, 1), period = 4))
summary(m4log)
tsdisplay(m4log$residuals)
add_model_resulst(m4log)

### log differencing twice ####

m1log2 = Arima(log(Denmark), order = c(0, 2, 0),
              seasonal = list(order = c(1, 0, 0), period = 4))
summary(m1log2)
tsdisplay(m1log2$residuals)
add_model_resulst(m1log2)

m2log2 = Arima(log(Denmark), order = c(1, 2, 0),
              seasonal = list(order = c(1, 0, 0), period = 4))
summary(m2log2)
tsdisplay(m2log2$residuals)
add_model_resulst(m2log2)

# par(mfrow = c(1, 1))
# qqnorm(m2log$residuals)
# qqline(m2log$residuals)
# shapiro.test(m2log$residuals)

m3log2 = Arima(log(Denmark), order = c(0, 2, 2),
              seasonal = list(order = c(1, 0, 0), period = 4))
summary(m3log2)
tsdisplay(m3log2$residuals)
add_model_resulst(m3log2)

par(mfrow = c(1, 1))
qqnorm(m3log$residuals)
qqline(m3log$residuals)
shapiro.test(m3log$residuals)

m4log = Arima(log(Denmark), order = c(1, 1, 0),
              seasonal = list(order = c(0, 1, 1), period = 4))
summary(m4log)
tsdisplay(m4log$residuals)
add_model_resulst(m4log)

## investing box cox ###########
# 
# m_boxcox <- Arima((Denmark), order = c(1, 1, 0),
#                   seasonal = list(order = c(0, 1, 1), period = 4),
#                   lambda = "auto")
# summary(m_boxcox)
# 
# m3_ar_seasonal <- Arima(log(Denmark), order = c(1, 1, 0),
#             seasonal = list(order = c(1, 1, 1), period = 4))
# summary(m3_ar_seasonal)
# tsdisplay(m3_ar_seasonal$residuals)
# add_model_resulst(m3_ar_seasonal)
# lrt(m3, m3_ar_seasonal)
# 
# anna_model <- Arima(log(Denmark), order = c(2, 1, 0),
#                         seasonal = list(order = c(1, 1, 1), period = 4))
# summary(anna_model)
# tsdisplay(anna_model$residuals)
# add_model_resulst(anna_model)
# lrt(m3, anna_model)
# 
# m7 <- Arima(log(Denmark), order = c(2, 1, 1),
#             seasonal = list(order = c(0, 1, 1), period = 4),
#             xreg = )
# summary(m7)
# c("AIC" = AIC(m7), "BIC" = BIC(m7))
# tsdisplay(m7$residuals)
# add_model_resulst(m7)

#### paramters: estimated AR, MA + intercept (ignore when differencing)
### + noise_sigma

m8 <- Arima(log(Denmark), order = c(2, 1, 0),
            seasonal = list(order = c(0, 0, 1), period = 4))
summary(m8)
# c("AIC" = AIC(m7), "BIC" = BIC(m7))
tsdisplay(m8$residuals)
add_model_resulst(m8)

m9 <- Arima(log(Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(1, 0, 1), period = 4))
summary(m9)
# c("AIC" = AIC(m7), "BIC" = BIC(m7))
tsdisplay(m9$residuals)
add_model_resulst(m9)

par(mfrow = c(1, 2))
qqnorm(m3$residuals)
qqline(m3$residuals)

qqnorm(m9$residuals)
qqline(m9$residuals)

shapiro.test(m9$residuals)
shapiro.test(m3$residuals)

binom.test(
  sum((m3$residuals[-1] * m3$residuals[-length(m3$residuals)]) < 0),
  length(m3$residuals) - 1
)

lrt(m9, m3)

predictions_m3 = predict(m3, n.ahead = 6)
exp(predictions_m3$pred)
exp(predictions_m3$se)

predictions_forcest = forecast(m3, h = 6)
exp(predictions_forcest$mean)
Box.test(m3$residuals, type = "Ljung-Box")

# model reduction ############
lrt <- function(model_small, model_big){
  test_stat = abs(as.numeric(-2*(logLik(model_small) - logLik(model_big))))
  df_difference = abs(length(model_big$coef) - length(model_small$coef))
  return(1 - pchisq(test_stat, df = df_difference))
}

lrt(m2, m3)
lrt(m3, m4)
lrt(m4, m5)
results
# condlucion: we should stay with the model m3
summary(m3)

# par(mfrow = c(1, 2))
# plot(Denmark)
# lines(exp(m3$fitted), col = "red")
# 
# m3_test <- Arima(log(Denmark[1:100]), order = c(1, 1, 0),
#             seasonal = list(order = c(0, 1, 1), period = 4))
# summary(m3_test)
# tsdisplay(m3_test$residuals)
# predictions <- predict(m3_test, n.ahead = 10)
# 

## Predictions with the final model ######
predictions <- predict(m3, n.ahead = 6)
predictions

m3

par(mfrow = c(1, 1))
plot(df$time, c(Denmark, exp(as.numeric(predictions$pred))))
abline(v = 122)
lines(, col = "red")

df$Sales %>%
  Arima(order = c(2, 0, 0),
        seasonal = list(order = c(1, 0, 1), period = period),
        include.mean = TRUE,
        fixed = parameters_2) %>%
  forecast(h=2) %>%
  autoplot

m3 %>% forecast(h=6) %>% autoplot
# model is trained for one step ahead forecast


# Additional regressors ####################
par(mfrow = c(2, 2))
# plot(InflationRate)
plot(diff(InflationRate, differences = 1))
plot(InflationRate)
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

plot(diff(df_train$InflationRate, differences = 1))
summary(ur.kpss(diff(df_train$InflationRate, differences = 1)))

acf(df$InflationRate)
pacf(df$InflationRate)

arima_reg <- Arima(InflationRate, order = c(3, 1, 1),
                   seasonal = list(order = c(0, 0, 1), period = 4))
summary(arima_reg)
# checkresiduals(m1)
# ggtsdisplay(log(Denmark))
tsdisplay(arima_reg$residuals)

## prewhitning ############
pwy_Denmark <- Arima(df_train$Denmark, model=arima_reg)
pwx <- arima_reg$residuals

ccf(pwx, pwy_Denmark$residuals, na.action=na.omit)
acf(cbind(pwx, pwy_Denmark$residuals), na.action=na.omit)

xreg_df = cbind(lag9x = stats::lag(ts(InflationRate),-9))

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

## 8.5.3 - doing it with the linear regression ##############################
# Notice, that it is assumed
# that both input and output series are mean-corrected.

# AICc comparisons must have the same orders
# of differencing. But RMSE test set comparisons
# can involve any models
# https://robjhyndman.com/uwafiles/8-Seasonal-ARIMA.pdf


attach(df)
library(stats)
lags_max = 7

InterestRate_ts = ts(df$InterestRate) - mean(df$InterestRate)
df_InterestRate = cbind(InterestRate_ts,
          stats::lag(InterestRate_ts,-1),
          stats::lag(InterestRate_ts,-2),
          stats::lag(InterestRate_ts,-3),
          stats::lag(InterestRate_ts,-4),
          stats::lag(InterestRate_ts,-5),
          stats::lag(InterestRate_ts,-6),
          stats::lag(InterestRate_ts,-7))

InflationRate_ts = ts(df$InflationRate) - mean(df$InflationRate)
df_InflationRate = cbind(InflationRate_ts,
          stats::lag(InflationRate_ts,-1),
          stats::lag(InflationRate_ts,-2),
          stats::lag(InflationRate_ts,-3),
          stats::lag(InflationRate_ts,-4),
          stats::lag(InflationRate_ts,-5),
          stats::lag(InflationRate_ts,-6),
          stats::lag(InflationRate_ts,-7))

df_xreg_identification = data.frame(cbind(df_InterestRate, df_InflationRate)) %>% na.omit()

# y = df$Denmark[(lags_max + 1):(nrow(df_xreg_identification) +lags_max)]
# lm_model = lm(y ~ ., data = df_xreg_identification)
# summary(lm_model)
# par(mfrow = c(2, 2))
# plot(lm_model)
# shapiro.test(lm_model$residuals)
# summary(lm_model, correlation = T)
# NO correlation of coefficients

# y_differenced = m1$fitted[(lags_max + 1):(length(m1$fitted) + 2)]
# # bc 2 additional observations in the independent regressors
# # wrong! it returns after the differencing
# 
# y_differenced = diff(log(df$Denmark), differences = 1, lag = 1)[(lags_max + 1):
#                                                                   (nrow(df_xreg_identification) +lags_max)]
# lm_model_2 = lm(y_differenced ~ ., data = df_xreg_identification)
# summary(lm_model_2)
# par(mfrow = c(2, 2))
# plot(lm_model_2)
# par(mfrow = c(1, 2))
# acf(lm_model_2$residuals)
# pacf(lm_model_2$residuals)
# model can be trused - inlude only InflationRate_ts without any lags               

# arima_full <- Arima(log(Denmark)[(lags_max + 1):(nrow(df_xreg_identification) +lags_max)],
#                     order = c(1, 1, 0),
#                     seasonal = list(order = c(0, 1, 1), period = 4),
#                     xreg = as.matrix(df_xreg_identification))
# summary(arima_full)
# HMM ...
y_testing <- diff(diff(log(df$Denmark), differences = 1, lag = 1), lag = 4)[(lags_max + 1 - 4):
                                                                              (nrow(df_xreg_identification) + 4 - 1)]
mean(y_testing[1:114])

lm_model_3 = lm(y_testing ~ ., data = df_xreg_identification)
summary(lm_model_3)

acf(lm_model_3$residuals)

drop1(lm_model_3, test = "F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InterestRate.stats..lag.InterestRate_ts...2.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InflationRate.stats..lag.InflationRate_ts...6.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . -df_InterestRate.stats..lag.InterestRate_ts...6.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InterestRate.stats..lag.InterestRate_ts...4.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InflationRate.stats..lag.InflationRate_ts...4.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InterestRate.stats..lag.InterestRate_ts...3.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InflationRate.stats..lag.InflationRate_ts...5.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InflationRate.stats..lag.InflationRate_ts...7.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InterestRate.stats..lag.InterestRate_ts...7.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InterestRate.stats..lag.InterestRate_ts...5.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InflationRate.stats..lag.InflationRate_ts...2.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InflationRate.stats..lag.InflationRate_ts...1.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InterestRate.stats..lag.InterestRate_ts...1.)
drop1(lm_model_3, test="F")
lm_model_3 <- update(lm_model_3, . ~ . - df_InterestRate.InterestRate_ts)
drop1(lm_model_3, test="F")

summary(lm_model_3)


options(scipen=0)
summary(lm_model_3)
# summary(lm_model_3, correlation = T)
# no colineartiy - no need to use the Ridge regression

par(mfrow = c(2, 2))
plot(lm_model_3)
par(mfrow = c(1, 2))
acf(lm_model_3$residuals)
pacf(lm_model_3$residuals)

# xreg_properly_differenced = cumsum(diffinv(df_train$InflationRate[1:length(df_train$Denmark)],4)[-(1:4)])
InterestRate_tmp = ts(cumsum(diffinv(df$InterestRate,4)[-(1:4)]))
InflationRate_tmp = ts(cumsum(diffinv(df$InflationRate,4)[-(1:4)]))

xreg_properly_differenced = cbind(
  "InterestRate_lag0" = InterestRate_tmp,
  "InterestRate_lag1" = stats::lag(InterestRate_tmp,-1),
  "InterestRate_lag2" = stats::lag(InterestRate_tmp,-2),
  "InterestRate_lag3" = stats::lag(InterestRate_tmp,-3),
  "InterestRate_lag4" = stats::lag(InterestRate_tmp,-4),
  "InterestRate_lag5" = stats::lag(InterestRate_tmp,-5),
  "InterestRate_lag6" = stats::lag(InterestRate_tmp,-6),
  "InterestRate_lag7" = stats::lag(InterestRate_tmp,-7),
  "InterestRate_lag8" = stats::lag(InterestRate_tmp,-8),
  "InterestRate_lag9" = stats::lag(InterestRate_tmp,-9),
  "InterestRate_lag10" = stats::lag(InterestRate_tmp,-10),
  "InterestRate_lag11" = stats::lag(InterestRate_tmp,-11),
  "InterestRate_lag12" = stats::lag(InterestRate_tmp,-12),
  "InflationRate_lag0" = InflationRate_tmp,
  "InflationRate_lag1" = stats::lag(InflationRate_tmp,-1),
  "InflationRate_lag2" = stats::lag(InflationRate_tmp,-2),
  "InflationRate_lag3" = stats::lag(InflationRate_tmp,-3),
  "InflationRate_lag4" = stats::lag(InflationRate_tmp,-4),
  "InflationRate_lag5" = stats::lag(InflationRate_tmp,-5),
  "InflationRate_lag6" = stats::lag(InflationRate_tmp,-6),
  "InflationRate_lag7" = stats::lag(InflationRate_tmp,-7),
  "InflationRate_lag8" = stats::lag(InflationRate_tmp,-8),
  "InflationRate_lag9" = stats::lag(InflationRate_tmp,-9),
  "InflationRate_lag10" = stats::lag(InflationRate_tmp,-10),
  "InflationRate_lag11" = stats::lag(InflationRate_tmp,-11),
  "InflationRate_lag12" = stats::lag(InflationRate_tmp,-12)
)

cat("InterestRate_lag", 1, sep = "")
max_lag = 13
c(toString(cat("InterestRate_lag", 1, sep = "")) = 5)

cat("InterestRate_lag", 1, sep = "")

u = 5
names(u) = "stuff"
u

df_proper_lags = array(NA, dim = )
for(i in seq(max_lag)){
  x_creating_stuff = stats::lag(InterestRate_tmp,-(i - 1))
  # names(x_creating_stuff) = cat("InterestRate_lag", (i - 1), sep = "")
  df_proper_lags = cbind(x_creating_stuff, df_proper_lags)
  # colnames(df_proper_lags) = c(cat("InterestRate_lag", (i - 1), sep = ""),
  #                              colnames(df_proper_lags))
}

data.frame(df_proper_lags)

test_stuff = ts(cumsum(diffinv(df$InflationRate,4)[-(1:4)]))
cbind(test_stuff, stats::lag(test_stuff,-1))

data.frame(df_InflationRate)

arima_xreg = Arima(log(df_train$Denmark), order = c(1, 1, 0),
                   seasonal = list(order = c(0, 1, 1), period = 4),
                   xreg = xreg_properly_differenced[1:nrow(df_train), ])
summary(arima_xreg)

arima_xreg_inflation = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)],
                             order = c(1, 1, 0),
                   seasonal = list(order = c(0, 1, 1), period = 4),
                   xreg = xreg_properly_differenced[lags_number:nrow(df_train),
                                                    "InflationRate_lag0"])
summary(arima_xreg_inflation)

lrt(arima_xreg, arima_xreg_inflation)
lrt(m3, arima_xreg_inflation)
lrtest(m3, arima_xreg_inflation)

lags_number = 6
################################################################
arima_xreg_inflation_lag02 = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
                             seasonal = list(order = c(0, 1, 1), period = 4),
                             xreg = xreg_properly_differenced[lags_number:nrow(df_train),
                                                              c("InflationRate_lag0")])
summary(arima_xreg_inflation_lag02)

# lrtest(m3, arima_xreg_inflation_lag02)

# lrtest(arima_xreg_inflation, arima_xreg_inflation_lag01)
# c(AIC(arima_xreg_inflation), AIC(arima_xreg_inflation_lag01))

##############################################################
lags_number = 1
arima_xreg_inflation_lag023 = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
                                   seasonal = list(order = c(0, 1, 1), period = 4),
                                   xreg = xreg_properly_differenced[lags_number:nrow(df_train),
                                                                    c("InflationRate_lag0",
                                                                      "InflationRate_lag2")])
summary(arima_xreg_inflation_lag023)
# lrtest(arima_xreg_inflation_lag02, arima_xreg_inflation_lag023)


arima_xreg_inflation_lag023_2 = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
                                    seasonal = list(order = c(0, 1, 1), period = 4),
                                    xreg = xreg_properly_differenced[lags_number:nrow(df_train),
                                                                     c("InflationRate_lag0",
                                                                       "InflationRate_lag2",
                                                                       "InflationRate_lag3")])
tsdisplay(arima_xreg_inflation_lag023_2$residuals)
tsdisplay(m3$residuals)

par(mfrow = c(1, 1))
qqnorm(arima_xreg_inflation_lag023_2$residuals)
qqline(arima_xreg_inflation_lag023_2$residuals)

shapiro.test(arima_xreg_inflation_lag023_2$residuals)

##############################################################################3

list_regressors_input = c(
  "InterestRate_lag0", "InterestRate_lag1", "InterestRate_lag2", "InterestRate_lag3", 
  "InterestRate_lag4",  "InterestRate_lag5", "InflationRate_lag0",
  "InflationRate_lag1", "InflationRate_lag2", "InflationRate_lag3",
  "InflationRate_lag4", "InflationRate_lag5"
)

# list_regressors_input_2 = c(
#   "InterestRate_lag0", "InterestRate_lag1", "InterestRate_lag2",   "InterestRate_lag3",  
#   "InterestRate_lag4","InterestRate_lag5",  "InterestRate_lag6",  "InterestRate_lag7",  
#   "InterestRate_lag8", "InterestRate_lag9", "InterestRate_lag10", "InterestRate_lag11", 
#   "InterestRate_lag12", "InflationRate_lag0", "InflationRate_lag1", "InflationRate_lag2", 
#   "InflationRate_lag3",  "InflationRate_lag4",  "InflationRate_lag5",  "InflationRate_lag6", 
#   "InflationRate_lag7",  "InflationRate_lag8",  "InflationRate_lag9",  "InflationRate_lag10",
#   "InflationRate_lag11", "InflationRate_lag12"
# )

arima_reduction <- function(list_regressors){
  results_regressors <- array(NA, dim = c(length(list_regressors), 1))
  rownames(results_regressors) = list_regressors
  lags_number = 6
  for(regressor in list_regressors){
    tmp_list_regressors = list_regressors[list_regressors != regressor]
    arima_xreg_full_tmp = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
                                seasonal = list(order = c(0, 1, 1), period = 4),
                                xreg = xreg_properly_differenced_standarized[lags_number:nrow(df_train), 
                                                                 list_regressors])
    arima_xreg_full_tmp_reduction = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
                                          seasonal = list(order = c(0, 1, 1), period = 4),
                                          xreg = xreg_properly_differenced_standarized[lags_number:nrow(df_train), 
                                                                           tmp_list_regressors])
    test_tmp = lrtest(arima_xreg_full_tmp, arima_xreg_full_tmp_reduction)
    # print(as.numeric(test_tmp$`Pr(>Chisq)`)[2])
    results_regressors[rownames(results_regressors) == regressor] = as.numeric(test_tmp$`Pr(>Chisq)`)[2]
  }
  testing_stuff = c(results_regressors)
  names(testing_stuff) = rownames(results_regressors)
  return(testing_stuff)
}

wrapper_function <- function(list_regressors, significant = 0.05){
  
  tmp_list_regressors = list_regressors
  results_LRT = arima_reduction(list_regressors)
  sorted_pvalues = sort(results_LRT, decreasing = T)
  i = 1
  
  while(sorted_pvalues[1] > significant){
    print(i)
    new_values = sorted_pvalues[names(sorted_pvalues) != names(sorted_pvalues[1])]
    tmp_list_regressors <- names(new_values)
    
    results_LRT = arima_reduction(tmp_list_regressors)
    sorted_pvalues = sort(results_LRT, decreasing = T)
    i = i + 1
  }
  return(sorted_pvalues)
}

results_pvalues_sorted_significant = wrapper_function(list_regressors_input)
results_pvalues_sorted_significant

# wrapper_function(list_regressors_input_2)
# 
# lags_final_model = 1
# list_regressors_input_24 = c("InterestRate_lag10", "InterestRate_lag11",
#                              "InflationRate_lag5", "InflationRate_lag0")
# 
# final_arima_model_24 = Arima(log(df_train$Denmark)[lags_final_model:nrow(df_train)], order = c(1, 1, 0),
#                           seasonal = list(order = c(0, 1, 1), period = 4),
#                           xreg = xreg_properly_differenced[lags_final_model:nrow(df_train), 
#                                                            list_regressors_input_24])
# summary(final_arima_model_24)
# tsdisplay(final_arima_model_24$residuals)

Box.test(final_arima_model_24$residuals, type = "Ljung-Box")


list_regressors_input = c(
  "InterestRate_lag1", 
  "InterestRate_lag4",
  "InflationRate_lag1",
  "InflationRate_lag4"
)



m3_testing <- Arima(log(Denmark)[lags_final_model:nrow(df_train)], order = c(1, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m3_testing)
tsdisplay(m3$residuals)
add_model_resulst(m3)

lrtest(m3_testing, final_arima_model)

# final_arima_model_testing = Arima(log(df_train$Denmark), order = c(1, 1, 0),
#                           seasonal = list(order = c(0, 1, 1), period = 4),
#                           xreg = xreg_properly_differenced[1:nrow(df_train), 
#                                                            list_regressors_input])
# summary(final_arima_model)

lrtest(final_arima_model, final_arima_model_testing)


# regressor = "InterestRate_lag0"
# arima_reduction(list_regressors_input)
# 
# list_regressors_input_2 = c(
#   "InterestRate_lag0", "InterestRate_lag1", "InterestRate_lag2",   "InterestRate_lag3",  
#   "InterestRate_lag4","InterestRate_lag5",  "InterestRate_lag6",  "InterestRate_lag7",  
#   "InterestRate_lag8", "InterestRate_lag9", "InterestRate_lag10", "InterestRate_lag11", 
#   "InterestRate_lag12", "InflationRate_lag0", "InflationRate_lag1", "InflationRate_lag2", 
#   "InflationRate_lag3",  "InflationRate_lag4",  "InflationRate_lag5",  "InflationRate_lag6", 
#   "InflationRate_lag7",  "InflationRate_lag8",  "InflationRate_lag9",  "InflationRate_lag10",
#   "InflationRate_lag11", "InflationRate_lag12"
# )




# lags_number = 12
tmp_list_regressors = list_regressors_input[list_regressors_input != regressor]

arima_xreg_full_tmp = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
                            seasonal = list(order = c(0, 1, 1), period = 4),
                            xreg = xreg_properly_differenced[lags_number:nrow(df_train), 
                                                             list_regressors_input])

arima_xreg_full_tmp_reduction = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
                                      seasonal = list(order = c(0, 1, 1), period = 4),
                                      xreg = xreg_properly_differenced[lags_number:nrow(df_train), 
                                                                       tmp_list_regressors])

test_tmp = lrtest(arima_xreg_full_tmp, arima_xreg_full_tmp_reduction)
lrt(arima_xreg_full_tmp, arima_xreg_full_tmp_reduction)
# print(as.numeric(test_tmp$`Pr(>Chisq)`)[2])
results_regressors[rownames(results_regressors) == regressor] = as.numeric(test_tmp$`Pr(>Chisq)`)[2]


# lags_number = 6
# arima_xreg_full = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
#                                       seasonal = list(order = c(0, 1, 1), period = 4),
#                                       xreg = xreg_properly_differenced[lags_number:nrow(df_train), ])
# 
# 
# arima_xreg_full_reduction = Arima(log(df_train$Denmark)[lags_number:length(df_train$Denmark)], order = c(1, 1, 0),
#                         seasonal = list(order = c(0, 1, 1), period = 4),
#                         xreg = xreg_properly_differenced[lags_number:nrow(df_train), 
#                                                          c("InflationRate_lag5")])




#################################################################################

lrtest(arima_xreg_inflation_lag023, arima_xreg_inflation_lag023_2)



m3 <- Arima(log(df_train$Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m3)
tsdisplay(m3$residuals)


lrt(m3, arima_xreg)
lrtest(m3, arima_xreg)


c(AIC(arima_xreg), BIC(arima_xreg))
# AIC is better, LRT doesn't reject the model
# interest rates is lagged vs inflationa rate ???
acf(cbind("Inflation_rate" = df$InflationRate,
          "interest_rate" = df$InterestRate), na.action=na.omit)
ccf(df$InflationRate, df$InterestRate, na.action=na.omit)
fccf(df$InterestRate, df$InflationRate, na.action=na.omit)

predictions_xreg <- predict(arima_xreg, n.ahead = 6,
                            xreg = tail(df[, ], 6))

# forecast.arima,h = 5, xreg=tf[114:length(tf)]

predictions <- a
exp(predictions)

arima_xreg %>%
  forecast(h=6, xreg = tail(xreg_properly_differenced, 6)) %>% 
  autoplot

%>%
  autoplot


m3 %>% forecast(h=6) %>% autoplot

plot(exp(c(predictions$fitted, predictions$mean)))
abline(v = 122)

plot(log(df_train$Denmark), type = "p")
lines(m3$fitted, col = "red")

acf(cbind(df$InflationRate, df$InterestRate))



xreg_predictions

##############################

xreg_properly_differenced_standarized = (xreg_properly_differenced - mean(xreg_properly_differenced,
                                                                          na.rm = T))/
  sd(xreg_properly_differenced, na.rm = T)

lags_final_model = 6
final_arima_model = Arima(df_train$Denmark[lags_final_model:nrow(df_train)], order = c(1, 1, 0),
                          seasonal = list(order = c(0, 1, 1), period = 4),
                          xreg = xreg_properly_differenced_standarized[lags_final_model:nrow(df_train), 
                                                           names(results_pvalues_sorted_significant)],
                          lambda = 0)
summary(final_arima_model)
tsdisplay(final_arima_model$residuals)
checkresiduals(final_arima_model)

final_arima_model_reducted = Arima(df_train$Denmark[lags_final_model:nrow(df_train)], order = c(1, 1, 0),
                          seasonal = list(order = c(0, 1, 1), period = 4),
                          xreg = xreg_properly_differenced_standarized[lags_final_model:nrow(df_train), 
                                                                       c(
                                                                         "InterestRate_lag4",
                                                                         "InflationRate_lag1",
                                                                         "InflationRate_lag4"
                                                                       )],
                          lambda = 0)
summary(final_arima_model_reducted)

lrtest(final_arima_model_reducted, final_arima_model)


list_regressors_input = c(
  "InterestRate_lag1", 
  "InterestRate_lag4",
  "InflationRate_lag1",
  "InflationRate_lag4"
)

m3 <- Arima(log(df_train$Denmark)[lags_final_model:nrow(df_train)], order = c(1, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m3)
lrtest(m3, final_arima_model)


####################################3

arima_dummy = Arima(df_train$Denmark[1:nrow(df_train)], order = c(1, 1, 0),
                          seasonal = list(order = c(0, 1, 1), period = 4),
                          xreg = xreg_properly_differenced_standarized[1:nrow(df_train), 
                                                           c("InterestRate_lag0", "InflationRate_lag0")],
                          lambda = 0)

summary(arima_dummy)
tsdisplay(arima_dummy$residuals)
checkresiduals(arima_dummy)

xreg_predictions <- final_arima_model %>%
  forecast(h=6, xreg = tail(xreg_properly_differenced_standarized[lags_final_model:(nrow(df_train)+6),
                                                      names(results_pvalues_sorted_significant)], 6))
predictions_dummy <- arima_dummy %>%
  forecast(h=6, xreg = tail(xreg_properly_differenced_standarized[1:(nrow(df_train)+6),
                                                      c("InterestRate_lag0", "InflationRate_lag0")], 6))

options(scipen=999)
results_pvalues_sorted_significant

xreg_predictions$mean
xreg_predictions$upper[,2]
xreg_predictions$lower[,2]

predictions_dummy$mean
predictions_dummy$upper[,2]
predictions_dummy$lower[,2]

##################333
n = nrow(df)
n_uni_res = nrow(df_train)
df$Time <- as.factor(df$X)

plot.new()
plot.window(xlim = c(1, n), ylim = c(0, max(c(df$InterestRate, df$InflationRate), na.rm=T)))
points(df$Time, df$InterestRate, col = "chartreuse3", type = "s")
points(df$Time, df$InflationRate, col="cadetblue", type="s")
axis(4)

plot.window(xlim=c(1, n), ylim=range(df$Denmark, na.rm=T))
lines(lags_final_model:(n_uni_res+1), 
      c(xreg_predictions$fitted, xreg_predictions$mean[1]), col = "black", lwd= 3)
lines((n_uni_res+1):(n_uni_res+6), xreg_predictions$mean[1:6], col = "black", lwd = 3)
lines((n_uni_res+1):(n_uni_res+6), xreg_predictions$lower[1:6], col = "black", lty = 2, lwd = 1)
lines((n_uni_res + 1):(n_uni_res + 6), xreg_predictions$upper[1:6], col = "black", lty = 2, lwd = 1)


lines(1:(n_uni_res+1), 
      c(predictions_dummy$fitted, predictions_dummy$mean[1]), col = "red", lwd= 3)
lines((n_uni_res+1):(n_uni_res+6), predictions_dummy$mean[1:6], col = "red", lwd = 3)
lines((n_uni_res+1):(n_uni_res+6), predictions_dummy$lower[1:6], col = "red", lty = 2, lwd = 1)
lines((n_uni_res + 1):(n_uni_res + 6), predictions_dummy$upper[1:6], col = "red", lty = 2, lwd = 1)

abline(v = n_uni_res + 1, lty = 1, col = "black")

points(df$Time, df$Denmark, col = "blue", pch = 20, type = "b")
axis(1, at=1:n, labels=df$Time, las=2)
axis(2)
box()
grid()

title("Average sales price, interest rate, and inflation rate", adj = 0)
mtext("Rate (%)", side = 4, las = 3, line = 3)
mtext("Average sales price (1000 DKK)", side = 2, las=3, line=3)
legend("top", legend=c("Forecast final model", "95% CI final model",
                       "Forecast naive model", "95% CI naive model",
                       "Denmark", "Interest rate", "Inflation rate"),
       col=c("black", "black", "red", "red", "blue",
             "chartreuse3", "cadetblue"),
       lty=c(1, 2, 1, 2, 1, 1, 1), cex=0.8)
grid()

