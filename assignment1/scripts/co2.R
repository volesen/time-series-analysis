rm(list=ls()) 
print(utils::getSrcDirectory(function(){}))
print(utils::getSrcFilename(function(){}, full.names = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MASS)
library(dplyr)
# library(tsibble)
library(forecast)

df <- read.table("A1_co2.txt", header = TRUE,
                 sep = " ", dec = ".")
df <- data.frame(index = seq(1, dim(df)[1]), df)

df_train <- df[1:(dim(df)[1]-20),]
df_test <- tail(df, n = 20)
c(dim(df_train), dim(df_test))
p_month <- 12

### 1.2 OLS ###  
lm_basic_train <- lm(co2 ~ index + sin(2*pi*index/p_month) +
                 cos(2*pi*index/p_month), data = df_train)
lm_basic_test <- lm(co2 ~ index + sin(2*pi*index/p_month) +
                       cos(2*pi*index/p_month), data = df_test)
xreg_basic_train <- model.matrix(lm_basic_train)
xreg_basic_test <- model.matrix(lm_basic_test)

# WLS
arima_basic <- arima(df_train$co2, xreg=xreg_basic_train[, seq(2, 4)], order=c(1,0,0))
wls_arima <- arima_basic$coef[c(2, 3, 4, 5)]

par(mfrow=c(1,1))
acf(arima_basic$residuals)
plot(arima_basic$residuals)

predict_wls <- xreg_basic_test %*% wls_arima
fitted_wls <- xreg_basic_train %*% wls_arima

################ 1.2 OLS and WLS ####################
#### 1. Estimate the parameters in OLS - BY HAND
X <- cbind(rep(1, n_train), df_train$time,
           sin(2*pi*(df_train$month - 1)/12),
           cos(2*pi*(df_train$month - 1)/12))
p <- dim(X)[2]
betas <- solve(t(X) %*% X) %*% t(X) %*% df_train$co2
betas_r <- solve(t(X) %*% X, t(X) %*% df_train$co2)

X_train <- cbind(rep(1, n_train), df_train$time,
                sin(2*pi*df_train$time),
                cos(2*pi*df_train$time))
y_train <- df_train$co2

X_test <- cbind(rep(1, (n - n_train)), df_test$time,
                sin(2*pi*df_test$time),
                cos(2*pi*df_test$time))

betas_time <- solve(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
betas_time_r <- solve(t(X_train) %*% X_train, t(X_train) %*% y_train)

#### 2. present the parameters with the uncertantiy
sigma_hat2 <- sum((y_train - X_train %*% betas_time_r)^2)/(dim(X_train)[1] - dim(X_train)[2])
Information_matrix <- t(X) %*% X / sigma_hat2
sd_beta <- sqrt(diag(solve(Information_matrix)))
summary(lm)

####### 4. Corelation structure, relaxation algorithm ######
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
ar1_cor(3, 0.7)

number <- 5
betas_matrix <- matrix(, nrow = p, ncol = number)
Sigma_matrix <- sigma_hat2*diag(n_train)
corelations <- rep(1, number)
for (i in (1:number)) {
  betas <- solve(t(X_train) %*% solve(Sigma_matrix) %*% X_train) %*% t(X_train) %*%
    solve(Sigma_matrix) %*% df_train$co2
  print(solve(t(X_train) %*% solve(Sigma_matrix) %*% X_train))
  print("")
  betas_matrix[,i] <- betas
  residuals <- df_train$co2 - X %*% betas
  df_residuals <- data.frame(e1 = residuals[-length(residuals)],
                             e2 = residuals[-1])
  corelation_iteration <- cor(df_residuals)[1, 2]
  corelations[i] <- corelation_iteration
  sigma_hat2_relaxation <- as.numeric(t(df_train$co2 - X_train %*% betas) %*%
                                 solve(Sigma_wls) %*% (df_train$co2 - X_train %*% betas) /
                                 (dim(X_train)[1] - dim(X_train)[2] - 1))
  print(sigma_hat2_relaxation)
  Sigma_matrix <- sigma_hat2_relaxation*ar1_cor(n_train, corelation_iteration)
}
betas_wls <- betas_matrix[,number]
Sigma_wls <- ar1_cor(n_train, corelations[number])

sigma_hat2_wls <- as.numeric(t(df_train$co2 - X_train %*% betas_wls) %*%
  solve(Sigma_wls) %*% (df_train$co2 - X_train %*% betas_wls) /
  (dim(X_train)[1] - dim(X_train)[2] - 1))
sd_beta_wls <- sqrt(diag(solve(t(X_train) %*% solve(Sigma_wls) %*% X_train /
                                 sigma_hat2_wls)))

### 5. Plotting
plot(df$co2, type = "l")
lines(c(X_train %*% betas_matrix[,1],
        X_test %*% betas_matrix[,1]), col = "red")
lines(c(X_train %*% betas_matrix[,number],
        X_test %*% betas_matrix[,number]), col = "blue")
abline(v = n_train)
legend("topleft",
       c("Data","OLS",
         "OLS AR"),
       col=c("black","red","blue"),lty=1,bty='n',lwd=2)

# Residuals - bad ... 
plot(df_train$co2 - X_train %*% betas_matrix[,1], type = "l")
lines(df_train$co2 - X_train %*% betas_matrix[,number], col = "blue")

################ 1.3 LOCAL TREND MODELS ####################
L <- rbind(
  c(1, 0, 0, 0),
  c(1, 1, 0, 0),
  c(0, 0, cos(2*pi/p_month), sin(2*pi/p_month)),
  c(0, 0, -sin(2*pi/p_month), cos(2*pi/p_month))
)
f <- function(j){
  return(cbind(1, j, sin(2*pi*j/p_month), cos(2*pi*j/p_month)))
}
lambda <- 0.9

X_N <- matrix(0, nrow =dim(df_train)[1], ncol = dim(f(0))[2])
for (i in 1:dim(df_train)[1]){
  X_N[i,] <- f(-dim(df_train)[1] + i)
}
betas_trend <- solve(t(X_N) %*% X_N) %*% t(X_N) %*% df_train$co2
predict_trend <- X_N %*% betas_trend

X_N_test <- matrix(0, nrow = dim(df_test)[1], ncol = dim(f(0))[2])
for (i in 1:dim(df_test)[1]){
  X_N_test[i,] <- f(i)
}

# the predictions of the linear trend model are the same as lm R
# we can define the inverse directly since it's just a diagonal matrix
Sigma_local_inverse <- diag(lambda^seq(dim(df_train)[1] - 1, 0))
# checking that it's correct
Sigma_local_inverse[(dim(df_train)[1]- 5):dim(df_train)[1],
                    (dim(df_train)[1]- 5):dim(df_train)[1]]

# epsilon <- 1e-15
# Epsilon_matrix <- diag(epsilon*seq(dim(df_train)[1], 1)/dim(df_train)[1]))
# # range(svd(solve(Sigma_local + Epsilon_matrix))$d)
# F_N_local <- t(X_N) %*% solve(Sigma_local + Epsilon_matrix) %*% X_N
# h_N_local <- t(X_N) %*% solve(Sigma_local + Epsilon_matrix) %*% df_train$co2

F_N_local <- t(X_N) %*% Sigma_local_inverse %*% X_N
h_N_local <- t(X_N) %*% Sigma_local_inverse %*% df_train$co2

beta_local <- solve(F_N_local) %*% h_N_local
fitted_local <- X_N %*% beta_local
predict_local <- X_N_test %*% beta_local

## PLOTTING
plot(df$co2, type = "l")
lines(c(fitted_wls, predict_wls), col = "red")
lines(c(fitted_local, predict_local), col = "blue")
abline(v = dim(df_train)[1])
legend("topleft",
       c("Data","OLS",
         "Local"),
       col=c("black","red","blue"),lty=1,bty='n',lwd=2)

###### 1.5 EXTENSIONS ##########
### INCLUDING TIME SQUARED #######
lm_test <- lm(co2 ~ index + I(index^2) + sin(2*pi*index/p_month) +
                cos(2*pi*index/p_month), data = df_train)
lm_orthogonal <- lm(co2 ~ poly(index, 2) + sin(2*pi*index/p_month) +
                      cos(2*pi*index/p_month), data = df_train)
cor(model.matrix(lm_test)[, 2], model.matrix(lm_test)[, 3])
cor(model.matrix(lm_orthogonal)[, 2],
    model.matrix(lm_orthogonal)[, 3])
summary(lm_orthogonal)
summary(lm_test)

# there is a difference in correlation of columns
plot(studres(lm_orthogonal))
abline(0, 0, col = "red")
acf(studres(lm_orthogonal))
Box.test(lm_orthogonal$residuals)

### Trying out the auto.arima ######
fit_auto_arima <- auto.arima(lm_basic_train$residuals)
Acf(fit_auto_arima$residuals)
arima_co2 <- auto.arima(df_train$co2)
plot(df$co2, type = "l")
lines(arima_co2$fitted, col = "red")

plot(lm_basic_train$residuals, type = "l")
lines(fit_auto_arima$fitted, col = "red")



