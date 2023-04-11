library(zoo)
library(MASS)
library(forecast)
library(tseries)
library(car)

data <- read.csv("./A3Data.csv")
data$Date <- as.Date(as.yearqtr(data$X, format="%YQ%q"))
train <- na.omit(data)

# 
Y <- data$Denmark


# Q1.1
par(mfrow=c(1, 3))
plot(train$Date, train$Denmark, type="b")
plot(train$Date, train$InterestRate, type="b")
plot(train$Date, train$InflationRate, type="b")

# Q1.2
dev.off()
# Split Y into groups of 4
Y_groups <- split(Y, ceiling(seq_along(Y)/8))

# Calculate the mean of each group
Y_groups_mean <- sapply(Y_groups, mean)

# Calculate the range of each group, i.e. max - min
Y_groups_max <- sapply(Y_groups, max)
Y_groups_min <- sapply(Y_groups, min)
Y_groups_range <- Y_groups_max - Y_groups_min


# Plot the mean range
plot(Y_groups_mean, Y_groups_range, ylim=c(0, 500))
text(Y_groups_mean, Y_groups_range, labels=1:length(Y_groups_mean), pos=4)
abline(lm(Y_groups_range ~ Y_groups_mean))

BoxCox.lambda(data$Denmark, method = "guerrero")
BoxCox.lambda(data$Denmark, method = "loglik")

library(MASS)
boxcox(lm(Y ~ 1))
boxcox(Y)


Y_transformed <- log(Y)

par(mfrow=c(1, 2))
acf(Y_transformed)
pacf(Y_transformed)

# sales are not stationary, lets try diff'ing once
kpss.test(Y_transformed)


lag1 <- diff(Y_transformed, lag=1)
par(mfrow=c(1, 3))
plot(lag1)
acf(lag1)
pacf(lag1)

checkresiduals(ts(lag1))

kpss.test(lag1)

# That worked, somewhat, the mean is constant now
# There is clearly a season with period 4

Y_season <- diff(lag1, lag=4)
par(mfrow=c(1, 3))
plot(Y_season)
pacf(Y_season)
acf(Y_season)

# Q1.3

Y_transformed %>% diff(lag=1) %>% ggtsdisplay()
Y_transformed %>% diff(lag=1) %>% diff(lag=4)  %>% ggtsdisplay()

# https://otexts.com/fpp2/seasonal-arima.html
# The significant spike at lag 1 in the ACF suggests a non-seasonal MA(1) component,
# and the significant spike at lag 4 in the ACF suggests a seasonal MA(1) component. 

# We start with ARIMA(0,1,1)(0,1,1)_4
# Q1.4
auto.arima(Y_transformed, D=4)

fit <- arima(Y_transformed, order=c(0,1,1), seasonal=list(order=c(0,1,1), period=4))
fit
checkresiduals(fit)

fit %>% residuals %>% ggtsdisplay()

# QQ plot of residuals
dev.off()
qqPlot(residuals(fit))


qqnorm(residuals(fit))
qqline(residuals(fit))

# Test the number of sign changes
sign_changes <- function(x) {
  sum(diff(sign(x)) != 0)
}

sign_changes(residuals(fit))
length(residuals(fit))

library(sarima)
whiteNoiseTest(acf(residuals(fit)), h0="iid", n=length(residuals(fit)))


#
fit <- auto.arima(log(train$Capital), allowdrift=FALSE)
fit
autoplot(forecast(fit))

checkresiduals(fit)

est <- exp(log(train$Capital) - residuals(fit))

plot(train$Date, train$Capital)
lines(train$Date, est, col="orange")


# Marima
model.marima <- define.model(
  kvar = 4,
  ar = 1:4,
  ma=0,
  #reg.var = c(5, 6),
  indep=c(1, 2, 3, 4),
  #rem.var = c(5, 6)
)

data.marima <- 
  t(log(subset(train, select = c(Capital, Sealand, MidJutland, Rural))))
  #,train$InterestRate,
  #,train$InflationRate


data.marima

difference = matrix(c(1, 1,
                      2, 1, 
                      3, 1, 
                      4, 1), nrow=2)
difference

data.diffed <- define.dif(data.marima, difference=difference)

N <- ncol(data.marima)
N

nsteps <- 6

nandf = data.frame(
  Capital = c(NaN, NaN, NaN, NaN, NaN, NaN),
  Sealand =  c(NaN, NaN, NaN, NaN, NaN, NaN),
  MidJutland =  c(NaN, NaN, NaN, NaN, NaN, NaN),
  Rural =  c(NaN, NaN, NaN, NaN, NaN, NaN)
)

data.pred <- cbind(data.marima, t(nandf))

data.pred

fit.marima <- marima(
  data.diffed$y.dif,
  ar.pattern = model.marima$ar.pattern,
  ma.pattern = model.marima$ma.pattern,
  Plot="log.det",
  Check=TRUE
)

forecasts <- arma.forecast(
  series = data.pred,
  marima = fit.marima, 
  dif.poly = data.diffed$dif.poly,
  nstart=N,
  nstep=nsteps,
  check=TRUE
)

plot(forecasts$forecasts[1,], type="l")
lines(log(train$Capital))


