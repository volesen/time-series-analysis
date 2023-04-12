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

par(mfrow = c(4, 4))
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

differenced_log_2 = diff(log(df_train$Denmark), differences = 2)
plot((df_train$time[-(1:2)]), differenced_log_2)
acf_2 <- acf(differenced_log_2, plot=F)
pacf_2 <- pacf(differenced_log_2, plot=F)
plot(acf_2)
plot(pacf_2)
hist(differenced_log_2, breaks = 20)

differenced_sqrt = diff(sqrt(df_train$Denmark), differences = 1)
plot((df_train$time[-1]), differenced_sqrt)
acf_2 <- acf(differenced_sqrt, plot=F)
pacf_2 <- pacf(differenced_sqrt, plot=F)
plot(acf_2)
plot(pacf_2)
hist(differenced_sqrt, breaks = 20)

## only seasonal differencing ###########
par(mfrow = c(3, 4))
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

model_boxcox_auto <- Arima((Denmark), order = c(1, 1, 0),
                        seasonal = list(order = c(0, 1, 1), period = 4),
                        lambda = "auto")
model_boxcox_auto
# lambda= -0.06814993 
# support for the log transformation

model_auto <- auto.arima(diff(log(Denmark), difference = 1))
summary(model_auto)
# add_model_resulst(model_auto)

### differencing regular and seasonal #######
m0 <- Arima(log(Denmark), order = c(0, 1, 0),
            seasonal = list(order = c(0, 0, 0), period = 4))
summary(m0)
summary(ur.kpss(m0$residuals))

m1 <- Arima(log(Denmark), order = c(0, 1, 0),
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

m4 <- arima(log(Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(1, 1, 1), period = 4))
summary(m4)
tsdisplay(m4$residuals)
add_model_resulst(m4)

lrt(m4, m3)

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

## AR vs differencing #####
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
qqnorm(m3log$residuals)
qqline(m3log$residuals)
shapiro.test(m3log$residuals)

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

arima_reg <- Arima(InflationRate, order = c(3, 1, 1),
                   seasonal = list(order = c(0, 0, 0), period = 0))
summary(arima_reg)
# checkresiduals(m1)
# ggtsdisplay(log(Denmark))
tsdisplay(arima_reg$residuals)

## prewhitning ############
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

## 8.5.3 - doing it with the linear regression #####
attach(df)
library(stats)
lags_max = 7

InterestRate_ts = ts(df$InterestRate)
df_InterestRate = cbind(InterestRate_ts,
          stats::lag(InterestRate_ts,-1),
          stats::lag(InterestRate_ts,-2),
          stats::lag(InterestRate_ts,-3),
          stats::lag(InterestRate_ts,-4),
          stats::lag(InterestRate_ts,-5),
          stats::lag(InterestRate_ts,-6),
          stats::lag(InterestRate_ts,-7))

InflationRate_ts = ts(df$InflationRate)
df_InflationRate = cbind(InflationRate_ts,
          stats::lag(InflationRate_ts,-1),
          stats::lag(InflationRate_ts,-2),
          stats::lag(InflationRate_ts,-3),
          stats::lag(InflationRate_ts,-4),
          stats::lag(InflationRate_ts,-5),
          stats::lag(InflationRate_ts,-6),
          stats::lag(InflationRate_ts,-7))


df_xreg_identification = data.frame(cbind(df_InterestRate, df_InflationRate)) %>% na.omit()


y = df$Denmark[(lags_max + 1):(nrow(df_xreg_identification) +lags_max)]
lm_model = lm(y ~ ., data = df_xreg_identification)
summary(lm_model)
par(mfrow = c(2, 2))
plot(lm_model)
shapiro.test(lm_model$residuals)
# summary(lm_model, correlation = T)
# NO correlation of coefficients

y_differenced = m1$fitted[(lags_max + 1):(length(m1$fitted) + 2)]
# bc 2 additional observations in the independent regressors
# wrong! it returns after the differencing

y_differenced = diff(log(df$Denmark), differences = 1, lag = 1)[(lags_max + 1):
                                                                  (nrow(df_xreg_identification) +lags_max)]
lm_model_2 = lm(y_differenced ~ ., data = df_xreg_identification)
summary(lm_model_2)
par(mfrow = c(2, 2))
plot(lm_model_2)
par(mfrow = c(1, 2))
acf(lm_model_2$residuals)
pacf(lm_model_2$residuals)
# model can be trused - inlude only InflationRate_ts without any lags               

# arima_full <- Arima(log(Denmark)[(lags_max + 1):(nrow(df_xreg_identification) +lags_max)],
#                     order = c(1, 1, 0),
#                     seasonal = list(order = c(0, 1, 1), period = 4),
#                     xreg = as.matrix(df_xreg_identification))
# summary(arima_full)
# HMM ...
y_testing <- diff(diff(log(df$Denmark), differences = 1, lag = 1), lag = 4)[(lags_max + 1 - 4):
                                                                              (nrow(df_xreg_identification) + 4 - 1)]
lm_model_3 = lm(y_testing ~ ., data = df_xreg_identification)
summary(lm_model_3)
# summary(lm_model_3, correlation = T)
# no colineartiy - no need to use the Ridge regression

par(mfrow = c(2, 2))
plot(lm_model_3)
par(mfrow = c(1, 2))
acf(lm_model_3$residuals)
pacf(lm_model_3$residuals)

# xreg_properly_differenced = cumsum(diffinv(df_train$InflationRate[1:length(df_train$Denmark)],4)[-(1:4)])
xreg_properly_differenced = cbind(
  cumsum(diffinv(df$InterestRate,4)[-(1:4)]),
  cumsum(diffinv(df$InflationRate,4)[-(1:4)])
  )
arima_xreg = Arima(log(df_train$Denmark), order = c(1, 1, 0),
                   seasonal = list(order = c(0, 1, 1), period = 4),
                   xreg = xreg_properly_differenced[1:nrow(df_train), ])
summary(arima_xreg)

m3 <- Arima(log(df_train$Denmark), order = c(1, 1, 0),
            seasonal = list(order = c(0, 1, 1), period = 4))
summary(m3)
tsdisplay(m3$residuals)


lrt(m3, arima_xreg)
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

predictions <- arima_xreg %>%
  forecast(h=6, xreg = tail(xreg_properly_differenced, 6))
exp(predictions)

%>%
  autoplot


m3 %>% forecast(h=6) %>% autoplot

plot(exp(c(predictions$fitted, predictions$mean)))
abline(v = 122)


