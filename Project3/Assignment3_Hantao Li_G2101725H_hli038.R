# @Author: Hantao Li hli038@e.ntu.edu.sg
# @Matric No.: G2101725H
# 
# This project is to analyze financial data. 
# The data are from the daily historical Apple stock prices(open, high, low, close and adjusted prices) 
# from February 1, 2002 to January 31, 2017 extracted from the Yahoo Finance website. 
# The data has logged the prices of the Apple stock everyday and comprises of the open, 
# close, low, high and the adjusted close prices of the stock for the span of 15 years.
# The goal of the project is to discover an interesting trend in the apple stock prices 
# over the past 15 years (3775 attributes) and to design and develop the best model for forecasting.


###### Load Data and Import the Libraries #########################################

setwd("C:/Users/Hanta/OneDrive - Nanyang Technological University/AI 6123 Time Series Analysis/Assignment/3")

library(TSA)
library(astsa)
library(zoo)
library(xts)
library(quantmod)
library(fBasics)
library(forecast)
library(ggplot2)
library(fGarch)
library(rugarch)
library(tseries)
library(zoo)

global.xlab <- 'Date'
global.ylab <- 'Adjusted Closing Price (USD)'
global.stockname <- 'AAPL'

data <- getSymbols(global.stockname, from='2002-02-01', to='2017-02-01', src='yahoo', auto.assign = F) 
data <- na.omit(data)
data_AC <- data[,4]

plot(data_AC, 
     main="Adjusted Closing Price of AAPL (2002-2017)", 
     xlab=global.xlab, 
     ylab=global.ylab)

### Basic attributes of the data
min(data_AC)  # 0.234286
max(data_AC)  # 33.25
mean(data_AC) # 10.77983

### ACF/PACF
ggtsdisplay(data_AC)

### Augmented Dickey-Fuller Test
adf.test(data_AC) # 0.4114

### Seasonal decomposition using stl()
ts_AC <- ts(Ad(to.monthly(data)), frequency = 12)
fit.stl <- stl(ts_AC[,1], s.window = "period")
autoplot(fit.stl, main="STL Decomposition")


###### Data transformation #########################################
lambda = 0 # log transform
percentage = 1 # change to 1 to remove multiplication
data_R <- diff(BoxCox(data_AC, lambda))*percentage # returns
data_R <- data_R[!is.na(data_R)]

### ACF/PACF
ggtsdisplay(data_R)
ggtsdisplay(abs(data_R))
ggtsdisplay(data_R^2)
adf.test(data_R) # 0.01

### QQ Plot
qqnorm(data_R)
qqline(data_R, col = 2)
skewness(data_R) # -0.190054
kurtosis(data_R) # 5.435621

### EACF
eacf(data_R)
eacf(abs(data_R))
eacf(data_R^2)

### Split dataset
train_num <- (length(data_R) - 30)
data_train <- head(data_R, train_num) # 3745
data_test <- tail(data_R, round(length(data_R) - train_num)) # 30


###### GARCH Fitting #########################################
Garch_11=garch(data_R, order=c(1,1))
summary(Garch_11)
AIC(Garch_11) # -18554.32

### GARCH Diagnostic Checking
checkresiduals(Garch_11)
qqnorm(residuals(Garch_11)); qqline(residuals(Garch_11), col = 2)
ggtsdisplay(abs(residuals(Garch_11)))
ggtsdisplay(residuals(Garch_11)^2)
gBox(Garch_11,method='squared')



###### RUGARCH Fitting #########################################

###### Distribution

### sGARCH(1,1), Normal Distribution (Innovations)
# 9296.009 -4.9229
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0))),
  data = data_R)

### sGARCH(1,1), Skew Normal Distribution
# 9296.242 -4.9225
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "snorm"),
  data = data_R)

### sGARCH(1,1), T-Distribution
# 9472.314 -5.0158
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### sGARCH(1,1), Skew T-Distribution
# 9472.777 -5.0155
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "sstd"),
  data = data_R)

### sGARCH(1,1), Generalized Error Distribution
# 9449.263 -5.0036
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "ged"),
  data = data_R)

### sGARCH(1,1), Skew Generalized Error Distribution
# 9450.581 -5.0038
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "sged"),
  data = data_R)

### sGARCH(1,1), Normal Inverse Gaussian Distribution
# 9468.302 -5.0131
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "nig"),
  data = data_R)

### sGARCH(1,1), Generalized Hyperbolic Distribution
# 9472.504 -5.0148
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "ghyp"),
  data = data_R)

### sGARCH(1,1), Johnson's S_U Distribution
# 9471.499 -5.0148
ugarchfit(spec = ugarchspec(
  variance.model=list(garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "jsu"),
  data = data_R)


######## Sub-model

### fGARCH(1,1), GARCH, T-Distribution
# 9472.314 -5.0158
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "GARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### fGARCH(1,1), TGARCH, T-Distribution
# 9494.4 -5.0270
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "TGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### fGARCH(1,1), AVGARCH, T-Distribution
# 9494.357 -5.0264
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "AVGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### fGARCH(1,1), NGARCH, T-Distribution
# 9482.028 -5.0204
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "NGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### fGARCH(1,1), NAGARCH, T-Distribution
# 9487.081 -5.0231
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "NAGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### fGARCH(1,1), APARCH, T-Distribution
# 9494.421 -5.0264
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "APARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### fGARCH(1,1), GJRGARCH, T-Distribution
# 9482.289 -5.0206
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "GJRGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### fGARCH(1,1), ALLGARCH, T-Distribution
# 9495.035 -5.0262
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "fGARCH",
                      submodel = "ALLGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### eGARCH(1,1), T-Distribution
# 9496.097 -5.0279
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "eGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### gjrGARCH(1,1), T-Distribution
# 9482.289 -5.0206
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "gjrGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### apARCH(1,1), T-Distribution
# 9494.421 -5.0264
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "apARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### iGARCH(1,1), T-Distribution
# 9471.993 -5.0162
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "iGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)

### csGARCH(1,1), T-Distribution
# 9486.041 -5.0220
ugarchfit(spec = ugarchspec(
  variance.model=list(model = "csGARCH",
                      garchOrder=c(1,1)),
  mean.model=list(armaOrder = c(0,0)),
  distribution.model = "std"),
  data = data_R)


###### Forecast #############################################

### N-roll
garchspec <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                        variance.model=list(model = "eGARCH", garchOrder=c(1,1)), 
                        distribution.model = "std")
fta <- ugarchfit(garchspec, data_R, out.sample=length(data_test))
fwdCast = ugarchforecast(fta, n.ahead=length(data_test), n.roll=length(data_test))
plot(fwdCast, which="all")


### Simulations
garchspec <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                        variance.model=list(model = "eGARCH", garchOrder=c(1,1),
                                            variance.targeting = FALSE), 
                        distribution.model = "std")
garchfit <- ugarchfit(data = data_train, spec = garchspec)
simgarchspec <- garchspec
setfixed(simgarchspec) <- as.list(coef(garchfit))
simgarch <- ugarchpath(spec = simgarchspec, m.sim = 5,
                       n.sim = 8 * 252, rseed = 123) 
simret <- fitted(simgarch)
plot.zoo(simret)
plot.zoo(sigma(simgarch))
simprices <- exp(apply(simret, 2, "cumsum"))
matplot(simprices, type = "l", lwd = 3)

