# @Author: Hantao Li hli038@e.ntu.edu.sg
# @Matric No.: G2101725H
# 
# The wwwusage  time series data consist of  the number of users connected to the internet through a server.
# The data are collected at a time interval of one minute and there are 100 ob?ervations.
# Please fit an appropriate ARIMA model for it and submit a short report including R codes, the fitted model, the diagnostic checking, AIC, etc.


###### Load Data and Import the Libraries #########################################

getwd()
setwd?"C:/Users/Hanta/OneDrive - Nanyang Technological University/AI 6123 Time Series Analysis/Assignment/1")
y = scan("./wwwusage.txt", skip = 1)
library(forecast)
library(tseries)
library(urca)
library(stats4)
library(sarima)
library(ggplot2)


## Check the st?tionary of TS ###########################################################

# Analysis the attribute of the original TS
min(y)
max(y)
mean(y)
plot(1:100, y, 
     xlim=c(0,100),
     ylim=c(50,250), 
     main="Number of users connected to internet over Tim?", 
     xlab="Time Point", 
     ylab="Number of users connected")
lines(1:100,y,type="l")

acf(y, lag.max = 100)  # does not cut off until lag 32 (dies down quickly)
pacf(y, lag.max = 100) # cut off after lag 2


# Calculate the properate diif to convert?TS
ndiffs(y) # 1
d1_y <- diff(y)
d2_y <- diff(d1_y)

plot(1:99, d1_y, 
     xlim=c(0,100),
     ylim=c(-15,15), 
     main="Time Plot after First-order Differencing", 
     xlab="Time Point", 
     ylab="First-order Difference")
lines(1:99,d1_y,type="l")

?lot(1:98, d2_y, 
     xlim=c(0,100),
     ylim=c(-15,15), 
     main="Time Plot after Second-order Differencing", 
     xlab="Time Point", 
     ylab="Second-order Difference")
lines(1:98,d2_y,type="l")

# ADF Augmented Dickey-Fuller
y %>% ur.kpss() %>% su?mary()
d1_y %>% ur.kpss() %>% summary()
d2_y %>% ur.kpss() %>% summary()

# KPSS Kwiatkowski-Phillips-Schmidt-Shin
y %>% ur.df() %>% summary()
d1_y %>% ur.df() %>% summary()
d2_y %>% ur.df() %>% summary()

# PP Philipps-Perron
PP.test(y)
PP.test(d1_y)
PP.t?st(d2_y)


## D1 #######################################################################

acf(d1_y, lag.max = 50)  # lag 24
pacf(d1_y, lag.max = 50) # lag 3

# Yule Walker to estimate AR coefficient
ar.yw(d1_y, order.max = 10) # 3

##### ARIMA(3,1,0) #####?
# Try arima(3,1,0) from yule walker est on differenced data
fit310 = arima(y, c(3,1,0)) 
tsdiag(fit310)
checkresiduals(fit310)
AIC(fit310) # 511.994
BIC(fit310) # 522.3745

## Q-Q Plot
qqnorm(fit310$residuals)
qqline(fit310$residuals)

Box.test(fit310$res?duals, type="Ljung-Box") # 0.9749

# Forecast for arima(3,1,0)
plot(forecast(fit310,h=20), 
     xlim=c(0,120), 
     ylim=c(-100,500), 
     main="Forecast of ARIMA(3, 1, 0)",
     xlab="Time Point", 
     ylab="Number of users connected")

# Validation f?r arima(3,1,0)
# 10-ahead forecasting and One-step forecasting with only historical data up to 90 mins
fit310_train = arima(y[1:90], c(3,1,0))
pred_310_1 = forecast(fit310_train)
fit310_test = Arima(y[91:100], model = fit310_train)
pred_310_2 = forecast(fi?310_test)
plot(c(y[1:90]),
     main = "Validation of ARIMA(3, 1, 0)",
     xlab = "Time Point",
     ylab = "Number of users connected",
     type = "l",
     xaxt = "n",
     xlim = c(0,100),
     ylim = c(50,250),
     panel.first = abline(h = seq(0,250?50), v = seq(0,120,10), col = "gray95"))
axis(1, at=seq(0,120,10), labels=seq(0,120,10)); lines(90:100, y[90:100], col="blue")
lines(90:100, c(y[90], pred_310_1$mean), col="red"); lines(90:100, c(y[90],pred_310_2$fitted), col="gold")


## D2 ##############?########################################################

acf(d2_y, lag.max = 50)  # lag 27
pacf(d2_y, lag.max = 50) # lag 2

# Yule Walker to estimate AR coefficient
ar.yw(d2_y, order.max = 10) # 2

##### ARIMA(2,2,0) #####

# Try arima(2,2,0) from yule w?lker est on differenced data
fit220 = arima(y, c(2,2,0)) 
tsdiag(fit220)
checkresiduals(fit220)
AIC(fit220) # 511.4645
BIC(fit220) # 519.2194

## Q-Q Plot
qqnorm(fit220$residuals)
qqline(fit220$residuals)

Box.test(fit220$residuals, type="Ljung-Box") # 0.8?11

# Forecasting for arima(2,2,0)
plot(forecast(fit220,h=20), 
     xlim=c(0,120), 
     ylim=c(-100,500), 
     main="Forecast of ARIMA(2, 2, 0)",
     xlab="Time Point", 
     ylab="Number of users connected")

# Validation for arima(2,2,0)
fit220_train?= arima(y[1:90], c(2,2,0))
pred_220_1 = forecast(fit220_train)
fit220_test = Arima(y[91:100], model = fit220_train)
pred_220_2 = forecast(fit220_test)
plot(c(y[1:90]),
     main = "Validation of ARIMA(2, 2, 0)",
     xlab = "Time Point",
     ylab = "Numbe? of users connected",
     type = "l",
     xaxt = "n",
     xlim = c(0,100),
     ylim = c(50,250),
     panel.first = abline(h = seq(0,250,50), v = seq(0,120,10), col = "gray95"))
axis(1, at=seq(0,120,10), labels=seq(0,120,10)); lines(90:100, y[90:100], ?ol="blue")
lines(90:100, c(y[90], pred_220_1$mean), col="red"); lines(90:100, c(y[90],pred_220_2$fitted), col="gold")


## AUTO #######################################################################

fitauto <- auto.arima(y)
fitauto

##### ARIMA(1,1,1) ##?##

# Try arima(1,1,1) from yule walker est on differenced data
tsdiag(fitauto)
checkresiduals(fitauto)
AIC(fitauto) # 514.2995
BIC(fitauto) # 522.0848

## Q-Q Plot
qqnorm(fitauto$residuals)
qqline(fitauto$residuals)

Box.test(fitauto$residuals, type="Ljun?-Box") # 0.8618

# Forecasting for arima(1,1,1)
plot(forecast(fitauto,h=20), 
     xlim=c(0,120), 
     ylim=c(-100,500), 
     main="Forecast of ARIMA(1, 1, 1)",
     xlab="Time Point", 
     ylab="Number of users connected")

# Validation for arima(1,1,1?
fit111_train = arima(y[1:90], c(1,1,1))
pred_111_1 = forecast(fit111_train)
fit111_test = Arima(y[91:100], model = fit111_train)
pred_111_2 = forecast(fit111_test)
plot(c(y[1:90]),
     main = "Validation of ARIMA(1, 1, 1)",
     xlab = "Time Point",
    ?ylab = "Number of users connected",
     type = "l",
     xaxt = "n",
     xlim = c(0,100),
     ylim = c(50,250),
     panel.first = abline(h = seq(0,250,50), v = seq(0,120,10), col = "gray95"))
axis(1, at=seq(0,120,10), labels=seq(0,120,10)); lines(90:10?, y[90:100], col="blue")
lines(90:100, c(y[90], pred_111_1$mean), col="red"); lines(90:100, c(y[90],pred_111_2$fitted), col="gold")


## Brute-Force #######################################################################

# 6¡Á6 matrix to store the AIC, BI?, and MASE of models
aics_1 <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6)) 
bics_1 <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6)) 
aics_2 <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6)) 
bics_2 <- matrix(NA, 7, 7, dimnames = list(p =?0:6, q = 0:6)) 

for(p in 0:6)
  for(q in 0:6)
    if(!(p == 0 & q == 0)){
      aics_1[p+1, q+1] <- AIC(arima(y, c(p, 1, q)))
      bics_1[p+1, q+1] <- BIC(arima(y, c(p, 1, q)))
      aics_2[p+1, q+1] <- AIC(arima(y, c(p, 2, q)))
      bics_2[p+1, q+1] <-?BIC(arima(y, c(p, 2, q)))
    }

# Find the location of minimum value
aics_1 - min(aics_1, na.rm = TRUE) # (5,1,4) (3,1,0) (5,1,5)
bics_1 - min(bics_1, na.rm = TRUE) # (1,1,1) (3,1,0) (1,1,2)
aics_2 - min(aics_2, na.rm = TRUE) # (5,2,5) (3,2,1) (2,2,0) 
bi?s_2 - min(bics_2, na.rm = TRUE) # (2,2,0) (0,2,3) (1,2,2) 


##### ARIMA(5,2,5) #####

# Try arima(5,2,5) from yule walker est on differenced data
fit525 = arima(y, c(5,2,5)) 
tsdiag(fit525)
checkresiduals(fit525)
AIC(fit525) # 509.8135
BIC(fit525) # 538.2?81

## Q-Q Plot
qqnorm(fit525$residuals)
qqline(fit525$residuals)

Box.test(fit525$residuals, type="Ljung-Box") # 0.6234

# Forecasting for arima(5,2,5)
plot(forecast(fit525,h=20), 
     xlim=c(0,120), 
     ylim=c(-100,500), 
     main="Forecast of ARIMA(?, 2, 5)",
     xlab="Time Point", 
     ylab="Number of users connected")

# Validation for arima(5,2,5)
fit525_train = arima(y[1:90], c(5,2,5))
pred_525_1 = forecast(fit525_train)
fit525_test = Arima(y[91:100], model = fit525_train)
pred_525_2 = forecast(?it525_test)
plot(c(y[1:90]),
     main = "Validation of ARIMA(5, 2, 5)",
     xlab = "Time Point",
     ylab = "Number of users connected",
     type = "l",
     xaxt = "n",
     xlim = c(0,100),
     ylim = c(50,250),
     panel.first = abline(h = seq(0,2?0,50), v = seq(0,120,10), col = "gray95"))
axis(1, at=seq(0,120,10), labels=seq(0,120,10)); lines(90:100, y[90:100], col="blue")
lines(90:100, c(y[90], pred_525_1$mean), col="red"); lines(90:100, c(y[90],pred_525_2$fitted), col="gold")


## Accuracy ######?################################################################

acctest <- window(y, start=91, end=100)
accuracy(forecast(fit310), acctest)
accuracy(forecast(fitauto), acctest)
accuracy(forecast(fit220), acctest)
accuracy(forecast(fit525), acctest)

#  (?,1,0)            ME      RMSE       MAE        MPE     MAPE      MASE         ACF1 Theil's U
# Training set  0.230588  3.044632  2.367157  0.2748377 1.890528 0.5230995 -0.003095066        NA
# Test set     -2.168862 12.071916 10.086919 -1.2910728 4.816833 ?.2290289  0.659781410  1.616196

#  (1,1,1)            ME      RMSE       MAE        MPE     MAPE      MASE         ACF1 Theil's U
# Training set  0.3035616  3.113754 2.405275  0.2805566 1.917463 0.5315228 -0.01715517        NA
# Test set     -2.5855588 11?351388 9.265086 -1.4666403 4.437330 2.0474185  0.64113200  1.477598

#  (2,2,0)            ME      RMSE       MAE        MPE     MAPE      MASE         ACF1 Theil's U
# Training set  0.02797758  3.150308  2.511921 0.2062350 1.994727 0.5550897 -0.0235521   ?    NA
# Test set      2.39255695 14.750798 13.118737 0.7787425 6.155561 2.8990065  0.6808541  2.186538

#  (5,2,5)            ME      RMSE       MAE        MPE     MAPE      MASE         ACF1 Theil's U
# Training set  0.1330466  2.713488 2.098567  0.21695?1 1.642098 0.4637459 0.04838289        NA
# Test set     -8.3673760 12.528193 9.415103 -4.1207123 4.581967 2.0805697 0.64399383  1.752434
