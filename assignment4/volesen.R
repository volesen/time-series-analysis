library(dplyr)
library(readr)

data <- read_csv("./A4_Kulhuse.csv")

print(data[rowSums(is.na(data)) > 0, ], n = 111)

# Are there missing values?
data %>%
  filter(if_any(everything(), is.na))

# So, there is missing rows
#datetimes.missing <-
#  data %>%
#  filter(if_any(everything(), is.na)) %>%
#  select(DateTime)

# The missing times are
# 2017-09-26 12:00:00 - 2017-09-28 16:00:00 (including)

# And spurious ones
# 2017-10-01 08:00:00
# 2017-10-02 17:30:00
# 2017-10-11 14:30:00
# 2017-10-11 21:30:00
# 2017-10-12 05:30:00

# Further, there is a NaN timestamp. Lets look into that

which(is.na(data$DateTime)) # index 4608

# That is between 2017-11-28 09:30:55 and 2017-11-28 10:00:00
# We can safely omit it


# Q1

# I used this tutroial to do the inset plots
# https://www.r-bloggers.com/2016/10/create-an-inset-plot/

# Lets start with choosing the

# Lets start with an inset plot
par(fig = c(0, 1, 0, 1))
plot(
  data$DateTime,
  data$Sal,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Salinity [PSU]'
)

margin <- 0.025
par(fig = c(0.5, 1 - margin, 0.2, 0.8), new = TRUE)
plot(
  data$DateTime,
  data$Sal,
  type = 'p',
  cex  = 0.25,
  xlim = as.POSIXct(c("2017-11-30", "2017-12-06")),
  ylim = c(17, 22),
  axes = FALSE,
  frame.plot = TRUE,
  ann = FALSE
)

# Plot 2
par(fig = c(0, 1, 0, 1))
plot(
  data$DateTime,
  data$ODO,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Dissolved oxygen [mg/L]'
)

margin <- 0.025
par(fig = c(0.5, 1 - margin, 0.05, 0.7), new = TRUE)
plot(
  data$DateTime,
  data$ODO,
  type = 'p',
  cex  = 0.25,
  xlim = as.POSIXct(c("2017-11-30", "2017-12-06")),
  ylim = c(10.5, 11.5),
  frame.plot = TRUE,
  ann = FALSE,
  axis = FALSE
)


# Q2
dev.off()
plot(
  data$DateTime,
  data$ODO,
  type = 'p',
  cex  = 0.25,
  xlim = as.POSIXct(c("2017-11-30", "2017-12-06")),
  ylim = c(10.5, 11.5),
  frame.plot = TRUE,
  ann = FALSE,
  axis = FALSE
)


# Q3
A <- matrix(c(1), nrow = 1)
C <- matrix(c(1), nrow = 1)
Sigma.1 <- matrix(c(0.01), nrow = 1) # System variance
Sigma.2 <- matrix(c(0.005), nrow = 1) # observation variance

source("./kalman.R")

X0 <- data$Sal[1]

k1 <-
  kalman(
    data$Sal,
    A = A,
    C = C,
    Sigma.1 = Sigma.1,
    Sigma.2 = Sigma.2,
    Xhat0 = X0,
    verbose = TRUE
  )

pred.sigma <- as.vector(k1$Sigma.yy.pred)

plot(
  data$DateTime,
  data$Sal,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Salinity [PSU]'
)
lines(data$DateTime, k1$pred[1:5000], col = "red")
lower = k1$pred - 1.96 * sqrt(pred.sigma)
upper = k1$pred + 1.96 * sqrt(pred.sigma)
lines(data$DateTime, lower[1:5000], col = "red", lty = 2)
lines(data$DateTime, upper[1:5000], col = "red", lty = 2)

# Standardized prediction error
residuals = (data$Sal - k1$pred[1:length(data$Sal), 1])
#residuals.zero = residuals - mean(residuals, na.rm=TRUE)
residuals.std = residuals / sqrt(k1$Sigma.xx.pred[1, 1, 1:length(data$Sal)])
plot(residuals.std)



source("./kalman.modified.R")

k2 <-
  kalman.modified(
    data$Sal,
    A = A,
    C = C,
    Sigma.1 = Sigma.1,
    Sigma.2 = Sigma.2,
    Xhat0 = X0,
    verbose = TRUE
  )


k2$X.outliers

pred.sigma <- as.vector(k1$Sigma.yy.pred)

plot(
  data$DateTime,
  data$Sal,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Salinity [PSU]'
)
lines(data$DateTime, k2$pred[1:5000], col = "red")
lower = k2$pred - 1.96 * sqrt(pred.sigma)
upper = k2$pred + 1.96 * sqrt(pred.sigma)
lines(data$DateTime, lower[1:5000], col = "red", lty = 2)
lines(data$DateTime, upper[1:5000], col = "red", lty = 2)

# Standardized prediction error
residuals = (data$Sal - k2$pred[1:length(data$Sal), 1])
#residuals.zero = residuals - mean(residuals, na.rm=TRUE)
residuals.std = residuals / sqrt(k2$Sigma.xx.pred[1, 1, 1:length(data$Sal)])
plot(residuals.std)

#  PART 5

log.lik <- function(par) {
  A <- matrix(1)
  C <- matrix(1)
  
  Sigma.1 <- matrix(exp(par[1])) # System variance
  Sigma.2 <- matrix(exp(par[2])) # observation variance
  
  X0 <- data$Sal[1]
  
  N <- 800
  k <-
    kalman(
      data$Sal[1:N],
      A = A,
      C = C,
      Sigma.1 = Sigma.1,
      Sigma.2 = Sigma.2,
      Xhat0 = X0,
      verbose = TRUE
    )
  
  vars <- k$Sigma.yy.pred[1, 1, 1:N]
  residuals <- data$Sal[1:N] - k$pred[1:N]
  
  return(0.5 * (sum(log(vars) + residuals ^ 2 / vars, na.rm = TRUE)))
  
}

low <- (0.005) ^ 2

vars.opt <-
  optim(
    par = log(c(0.01, 0.005)),
    fn = log.lik,
    method = "L-BFGS-B",
    lower = log(c(0, low)),
    hessian = TRUE
  )

exp(vars.opt$par)
