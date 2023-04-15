library(dplyr)
library(readr)

data <- read_csv("./A4_Kulhuse.csv")

print(dataset[rowSums(is.na(dataset)) > 0, ], n = 111)

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
par(fig = c(0.5, 1-margin, 0.2, 0.8), new = TRUE)
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
par(fig = c(0.5, 1-margin, 0.05, 0.7), new = TRUE)
plot(
  data$DateTime,
  data$ODO,
  type = 'p',
  cex  = 0.25,
  xlim = as.POSIXct(c("2017-11-30", "2017-12-06")),
  #ylim = c(10.5, 11.5),
  #frame.plot = TRUE,
  #ann = FALSE
)


axis(1, at = seq(as.POSIXct("2017-11-30"), as.POSIXct("2017-12-06"), by = "day"), format = "%d-%m-%Y")



# Q2
seq(as.POSIXct("2017-11-30"), as.POSIXct("2017-12-06"), by = "day")


