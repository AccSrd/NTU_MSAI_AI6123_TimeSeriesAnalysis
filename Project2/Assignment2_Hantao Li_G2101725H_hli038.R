# @Author: Hantao Li hli038@e.ntu.edu.sg
# @Matric No.: G2101725H
# 
# This data is for monthly anti-diabetic drug sales in Australia from 1992 to 2008. 
# Total monthly scripts for pharmaceutical products falling under ATC code A10, as recorded by the Australian Health Insurance Commission. 
# Please build a good model to predict the drug sales.


###### Load Data and Import the Libraries #########################################

getwd()
setwd("C:/Users/Hanta/OneDrive - Nanyang Technological University/AI 6123 Time Series Analysis/Assignment/2")
library(forecast)
library(tseries)
library(urca)
library(stats4)
library(sarima)
library(ggplot2)
library(dplyr)
library(lubridate)
library(astsa)
library(gridExtra)

data = read.delim("./drug.txt", header=TRUE, sep = ",")

data <- data %>% mutate(date = ymd(date))
data.start_time <<- ymd(as.Date(data$date[[1]])) %>% print()
data.end_time <<- ymd(as.Date(data$date[[nrow(data)]])) %>% print()
data.freq <<- 12

### Convert to time series data
data_all <- ts(data$value, start=c(year(data.start_time), month(data.start_time)),
               end=c(year(data.end_time), month(data.end_time)), frequency=data.freq)
plot(data_all, 
     main="Monthly Anti-diabetic Drug Sales in Australia", 
     xlab="Year", 
     ylab="Drug Sales")

### Basic attributes of the data
min(data_all)  # 2.81452
max(data_all)  # 29.66536
mean(data_all) # 10.69443

data_all %>% ur.kpss() %>% summary()
v %>% ur.df() %>% summary()

### Seasonal decomposition using stl()
autoplot(stl(data_all, s.window="period"))

### ACF/PACF
ggtsdisplay(data_all)



###### Split dataset #########################################

data_train <- head(data_all, round(0.8 * length(data_all)))
data_test <- tail(data_all, round(length(data_all) - round(0.8 * length(data_all))))
length(data_train) # 163
length(data_test)  # 41


###### Find suitable SARIMA model #########################################

## Find lambda then transform data using BoxCox to reduce variance
lambda_zero <- 0
lambda_auto <- BoxCox.lambda(data_all, 'guerrero') %>% print() # 0.1313326

data_train_t_zero <- BoxCox(data_train,lambda_zero)
data_train_t_auto <- BoxCox(data_train,lambda_auto)

autoplot(stl(data_train_t_zero, s.window="period"))
autoplot(stl(data_train_t_auto, s.window="period"))

# which is better?
data_train_t <- data_train_t_zero
# data_train_t <- data_train_t_auto

# diff
data_train_lag = diff(data_train_t, lag = data.freq)
ggtsdisplay(data_train_lag,
            main  = "Lag-12 Seasonal Differencing",
            ylab  = "Drug Sales",
            xlab  = "Year")

data_train_diff_lag = diff(data_train_diff, lag = data.freq)
ggtsdisplay(data_train_diff_lag,
            main  = "Lag-12 + Diff-1 Seasonal Differencing",
            ylab  = "Drug Sales",
            xlab  = "Year")

# So, we choose D=1, d=1
# p, q = 5(PACF), 4(ACF) (7?)
# P, Q = 0, 2


############ ARMA model ############
model.514012 = stats::arima(data_train_t, 
                            order=c(5,1,4), 
                            seasonal=list(order=c(0,1,2), 
                                          period=data.freq))
checkresiduals(model.514012)
tsdiag(model.514012)
qqnorm(model.514012$residuals)
qqline(model.514012$residuals)
AIC(model.514012) # -441.2474
BIC(model.514012) # -405.1198
summary(model.514012)
#    RMSE        MAE
# 0.04691566 0.03377543

pred.514012 <- InvBoxCox(forecast(model.514012, h = length(data_test)+20)$mean, lambda = lambda_zero)
plot(pred.514012,
     main="Prediction of Drug Sales of SARIMA(5,1,4)(0,1,2)", 
     xlab="Year", 
     ylab="Drug Sales")
lines(data_test, col="blue");
accuracy(pred.514012, data_test)
#    RMSE        MAE
#  1.500039    1.217112


############ AR model ############
model.510010 = stats::arima(data_train_t, 
                            order=c(5,1,0), 
                            seasonal=list(order=c(0,1,0), 
                                          period=data.freq))
checkresiduals(model.510010)
tsdiag(model.510010)
qqnorm(model.510010$residuals)
qqline(model.510010$residuals)
AIC(model.510010) # -416.2912
BIC(model.510010) # -398.2274
summary(model.510010)
#    RMSE        MAE
# 0.05547478 0.04048461

pred.510010 <- InvBoxCox(forecast(model.510010, h = length(data_test)+20)$mean, lambda = lambda_zero)
plot(pred.510010,
     main="Prediction of Drug Sales of SARIMA(5,1,0)(0,1,0)", 
     xlab="Year", 
     ylab="Drug Sales")
lines(data_test, col="blue");
accuracy(pred.510010, data_test)
#    RMSE        MAE
#  3.400101    2.956103


############ MA model ############
model.014012 = stats::arima(data_train_t, 
                            order=c(0,1,4), 
                            seasonal=list(order=c(0,1,2), 
                                          period=data.freq))
checkresiduals(model.014012)
tsdiag(model.014012)
qqnorm(model.014012$residuals)
qqline(model.014012$residuals)
AIC(model.014012) # -441.8861
BIC(model.014012) # -420.8116
summary(model.014012)
#    RMSE        MAE
# 0.04889734 0.03598504

pred.014012 <- InvBoxCox(forecast(model.014012, h = length(data_test)+20)$mean, lambda = lambda_zero)
plot(pred.014012,
     main="Prediction of Drug Sales of SARIMA(0,1,4)(0,1,2)", 
     xlab="Year", 
     ylab="Drug Sales")
lines(data_test, col="blue");
accuracy(pred.014012, data_test)
#    RMSE        MAE
#  1.545632   1.257072



############ ARMA model with auto lambda ############
model.514012_auto = stats::arima(data_train_t_auto, 
                                 order=c(5,1,4), 
                                 seasonal=list(order=c(0,1,2), 
                                               period=data.freq))
checkresiduals(model.514012_auto)
tsdiag(model.514012_auto)
qqnorm(model.514012_auto$residuals)
qqline(model.514012_auto$residuals)
AIC(model.514012_auto) # -366.7601
BIC(model.514012_auto) # -330.6324
summary(model.514012_auto)
#    RMSE        MAE
# 0.05841829 0.04262249

pred.514012_auto <- InvBoxCox(forecast(model.514012_auto, h = length(data_test)+20)$mean, lambda = lambda_auto)
plot(pred.514012_auto,
     main="Prediction of Drug Sales of SARIMA(5,1,4)(0,1,2) with auto-lambda", 
     xlab="Year", 
     ylab="Drug Sales")
lines(data_test, col="blue");
accuracy(pred.514012_auto, data_test)
#    RMSE        MAE
#  1.664577    1.37995


## Brute-Force ################################################

# 6×6 matrix to store the AIC, BIC of models
aic_mat <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6)) 
bic_mat <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6))
rsme_mat <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6))
mae_mat <- matrix(NA, 7, 7, dimnames = list(p = 0:6, q = 0:6))

for(p in 0:6)
  for(q in 0:6)
    if(!(p == 0 & q == 0)){
      model.cur = stats::arima(data_train_t,
                               order=c(p,1,q),
                               seasonal=list(order=c(0,1,2), 
                                             period=data.freq))
      pred.cur <- InvBoxCox(forecast(model.cur, h = length(data_test)+20)$mean, lambda = lambda_zero)
      aic_mat[p+1, q+1] <- AIC(model.cur)
      bic_mat[p+1, q+1] <- BIC(model.cur)
      rsme_mat[p+1, q+1] <- accuracy(pred.cur, data_test)[2]
      mae_mat[p+1, q+1] <- accuracy(pred.cur, data_test)[3]
    }

# Find the location of minimum value
aic_mat - min(aic_mat, na.rm = TRUE)
bic_mat - min(bic_mat, na.rm = TRUE)
rsme_mat - min(rsme_mat, na.rm = TRUE)
mae_mat - min(mae_mat, na.rm = TRUE)


############ best AIC&BIC model ############
model.213012 = stats::arima(data_train_t, 
                            order=c(2,1,3), 
                            seasonal=list(order=c(0,1,2), 
                                          period=data.freq))

pred.213012 <- InvBoxCox(forecast(model.213012, h = length(data_test)+20)$mean, lambda = lambda_zero)
plot(pred.213012,
     main="Prediction of Drug Sales of SARIMA(2,1,3)(0,1,2)", 
     xlab="Year", 
     ylab="Drug Sales")
lines(data_test, col="blue");


## Holt-Winter ################################################

model.HW_add <- hw(data_train_t, seasonal = "additive")
model.HW_mul <- hw(data_train_t, seasonal = "multiplicative")

pred.HW_add <- InvBoxCox(forecast(model.HW_add$mean, h=length(data_test))$mean, lambda = lambda_zero)
pred.HW_mul <- InvBoxCox(forecast(model.HW_mul$mean, h=length(data_test))$mean, lambda = lambda_zero)

plot(pred.HW_add,
     main="Prediction of Drug Sales of Holt-Winter Model", 
     xlab="Year", 
     ylab="Drug Sales")
lines(pred.HW_mul, col='gold')
lines(data_test, col="blue");

accuracy(pred.HW_add, data_test) %>% print()
#    RMSE        MAE
#  1.857903   1.593702

accuracy(pred.HW_mul, data_test) %>% print()
#    RMSE        MAE
#  3.258563 2.777955