library(dplyr)
library(readr)
library(grid)
library(ggplot2)

data <- read_csv("./A4_Kulhuse.csv")

# Are there missing values?
data %>%
  filter(if_any(everything(), is.na))

# So, there is missing rows
datetimes.missing <-
  data %>%
  filter(if_any(everything(), is.na)) %>%
  select(DateTime)

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


# Q1 ----
# Lets start with choosing the last week of data
last_week <- data %>% slice_tail(n = 2 * 24 * 7)

par(mfrow=c(1, 2))

ggplot(data,
       aes(DateTime, Sal)) +
  geom_point(size = 0.5) +
  xlab("Time [30 min.]") +
  ylab("Salinity [PSU]") +
  theme_light()


ggplot(last_week,
       aes(DateTime, Sal)) +
  geom_point(size = 0.5) +
  xlab("Time [30 min.]") +
  ylab("Salinity [PSU]") +
  theme_bw()


# Lets start with an inset plot
png(filename="sal.png", units="mm", width=200, height=100, res=300)
par(mfrow = c(1, 2))
plot(
  data$DateTime,
  data$Sal,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Salinity [PSU]'
)
plot(
  last_week$DateTime,
  last_week$Sal,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Salinity [PSU]'
)
dev.off()

# Plot 2
png(filename="odo.png", units="mm", width=200, height=100, res=300)
par(mfrow = c(1, 2))
plot(
  data$DateTime,
  data$ODO,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Dissolved oxygen [mg/L]'
)
plot(
  last_week$DateTime,
  last_week$ODO,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Dissolved oxygen [mg/L]'
)
dev.off()

# Q2 ----
# See latex

# Q3 ----
source("./kalman.R")

A <- matrix(1)
C <- matrix(1)

Sigma.1 <- matrix(0.01) # System variance
Sigma.2 <- matrix(0.005) # observation variance

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

## One-step predictions and residuals ----
pred = k1$pred[1:length(data$Sal)]
pred_sigma = sqrt(k1$Sigma.yy.pred[1:length(data$Sal)])

lower = pred - 1.96 * pred_sigma
upper = pred + 1.96 * pred_sigma

par(mfrow=c(1,1))
plot(
  data$DateTime,
  data$Sal,
  type = 'p',
  cex  = 0.25,
  xlab = 'Time [30 min.]',
  ylab = 'Salinity [PSU]'
)
lines(data$DateTime, pred, col = "red")
lines(data$DateTime, lower, col = "red", lty = 2)
lines(data$DateTime, upper, col = "red", lty = 2)



##

# Standardized prediction error
residuals = (data$Sal - k1$pred[1:length(data$Sal), 1])
#residuals.zero = residuals - mean(residuals, na.rm=TRUE)
residuals.std = residuals / sqrt(k1$Sigma.xx.pred[1, 1, 1:length(data$Sal)])
plot(residuals.std)


# Q4 ----
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

opt <-
  optim(
    par = log(c(0.00, 0.005)),
    fn = log.lik,
    method = "L-BFGS-B",
    lower = log(c(0, low)),
    hessian = TRUE
  )

opt.par <- exp(opt$par)
opt.par

opt.std <- sqrt(diag(solve(opt$hessian)))
exp(opt.std)
exp(opt$par[1] + c(-1.96, 1.96) * opt.std[1])
exp(opt$par[2] + c(-1.96, 1.96) * opt.std[2])
