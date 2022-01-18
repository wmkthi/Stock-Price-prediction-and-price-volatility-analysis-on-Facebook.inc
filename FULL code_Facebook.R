setwd("C:\\Users\\Wimukthi\\OneDrive - University of Moratuwa\\Semester 7\\DA4420 - Financial Econometrics\\Take Home Assignment")

#Libraries
library("zoo")
library("strucchange")
library("fpp2")
library("tseries")
library("hydroTSM")
library("xts")
library("MuMIn")
library("aTSA")
library("fGarch")
library("rugarch")

#Import data
fb <- read.csv('fb.csv')

#Handling missing values
dates = data.frame(Date = seq(as.Date('2015-01-02'), as.Date('2021-01-11'), by = 'days'))
fb$Date = as.Date(fb$Date, "%Y-%m-%d")
fb_full <- merge(fb,dates,by="Date", all = T)
FB <-na.locf(fb_full, fromLast = TRUE)

#ARIMA model for predicting Stock prices
close=FB[,"Adj.Close"]
close_log = log(close)#log transformation
close_vals <- ts(data = close_log, frequency = 356, start = c(2015,1))
plot(close_vals, col = "Blue4")

#Check for structural breaks
close_brk <- breakpoints(close_vals~1, h=0.1)
summary(close_brk)
plot(close_brk)

##consider 7 structural breaks
breakdates(close_brk, breaks =7)
plot(close_vals)
lines(fitted(close_brk, breaks=7), col=4)
lines(confint(close_brk, breaks=7))

FB$break1 <- 1
FB$break1[1:257] <- 0
FB$break2 <- 1
FB$break2[1:482] <- 0
FB$break3 <- 1
FB$break3[1:787] <- 0
FB$break4 <- 1
FB$break4[1:1007] <- 0
FB$break5 <- 1
FB$break5[1:1338] <- 0
FB$break6 <- 1
FB$break6[1:1558] <- 0
FB$break7 <- 1
FB$break7[1:1965] <- 0
FB$close_log <- close_log

lm1<-lm(close_log~break1+break2+break3+break4+break5+break6+break7,data = FB)
summary(lm1)

##consider single structural break
plot(close_vals)
lines(fitted(close_brk, breaks=1), col=4)
lines(confint(close_brk, breaks=1))

resids1 <- lm1$residuals
adf.test(resids1)

#Create new dataset after controlling for structural breaks
close_vals_latest <-close_vals[798:2202]
close_vals_latest_ts <- ts(data = close_vals_latest, frequency = 356, start = c(2017,86))
plot(close_vals_latest_ts, col = "Blue4", ylab = "Adjusted close price")

#Check for monthly seasonal effects
date = seq(as.Date('2017-03-9'), as.Date('2021-01-11'), by = 'days')
means<- monthlyfunction(data.matrix(close_vals_latest_ts), mean, na.rm=TRUE, dates = date,out.type = 'data.frame')
M_means <- c(means[1],means[2],means[3],means[4],means[5],means[6],means[7],means[8],means[9],means[10],means[11],means[12])

plot(M_means, ylim = c(4,6), type = "b", col = "Red",xaxt = 'n', ylab='Monthly averages',xlab = 'Month') #plot monthly averages
axis(1, at=1:12, labels=c("March","April","May","June","July","August","September","October","November","December","January","February"), cex.axis=1)

#Box-Jenkins methology for estimating ARIMA model
##Identification
###initial guess for d
adf.test(close_vals_latest_ts) #Check for stationarity
d.close_vals=diff(close_vals_latest_ts,differences = 1) #Differening
adf.test(d.close_vals) #check for stationarity of the difference data
plot.ts(d.close_vals, col = "bLUE4", ylab = "Close price returns (differenced close price)")

###Initial guesses for p and q
acf(d.close_vals, lag=30, ylim = c(-0.1,0.2))
pacf(d.close_vals, lag=30)

##Estimation
arima_1 <- Arima(close_vals_latest_ts, order=c(1,1,1))
arima_1
fit <- auto.arima(close_vals_latest_ts, seasonal=FALSE) #auto.arima function to ensure the model selection
fit 

###Information criterion 
AIC(arima_1,fit)
BIC(arima_1,fit)
AICc(arima_1,fit)

##Diagnostic checking
checkresiduals(arima_1)

##Forecast
###train test split
close_vals.train <- window(close_vals_latest_ts, end=c(2020,8))#train set
close_vals.test <- window(close_vals, start=c(2020,9))#test set

###model training
arima.train <- Arima(close_vals.train, order=c(1,1,1), include.drift = F) #ARIMA (1,1,1)
arima.train1 <- Arima(close_vals.train, order=c(0,1,1), include.drift = T) #ARIMA (0,1,1) with a drift
###evaluation metrics
accuracy(forecast(arima.train, h=414), close_vals.test) #ARIMA (1,1,1)
accuracy(forecast(arima.train1, h=414), close_vals.test) #ARIMA (0,1,1) with a drift

###Predicted results
forecast <-forecast(arima.train, h=414) #ARIMA (1,1,1)
autoplot(close_vals) + autolayer(forecast, series = "ARIMA(1,1,1)\nforecasts", alpha = 0.5) 
forecast1 <-forecast(arima.train1, h=414) #ARIMA (0,1,1) with a drift
autoplot(close_vals) + autolayer(forecast1, series = "ARIMA(0,1,1)\nforecasts with a drift", alpha = 0.5) 


#Multivariate model for predicting stock prices

#Variables
##Adjusted close
adj.close=FB[,"Adj.Close"]
adj.close_log = log(adj.close)#log transformation
adj.close_log_vals <- ts(data = adj.close_log, frequency = 356, start = c(2015,1))
##Volume
volume=FB[,"Volume"]
volume_log = log(volume)#log transformation
volume_vals <- ts(data = volume_log, frequency = 356, start = c(2015,1))


plot(adj.close_log_vals, type="l", xlab="January 02nd, 2015 to January 11th, 2021", ylab="Stock prices", col="blue",lwd = 2) #adjustedclose
par(new=T)
plot(volume_vals, type="l",axes=F, xlab="", ylab="", col = "gray19") #volume
legend("topleft",c("Adjusted close","Volume"),fill=c("blue","black"))
par(new=F)

#check for stationary
adf.test(adj.close_log_vals) 
adf.test(volume_vals)

#differencing Adj close
d.adj.close=diff(adj.close_log_vals,differences = 1)
adf.test(d.adj.close)

#Set the VAR dataset
v1 <- cbind(d.adj.close, volume_vals)
v1 <- v1[-1,]
colnames(v1) <- cbind("Adj.close","volume")

#Lag selection
lagselect <- VARselect(v1, lag.max = 15, type = "const")
lagselect$selection 

#Model estimation
Model1 <- VAR(v1, p = 6, type = "const", season = NULL, exog = NULL) 
summary(Model1)

#Granger causality
cause.Adj.close <- causality(Model1, cause = "Adj.close")
cause.Adj.close
cause.volume <- causality(Model1, cause = "volume")
cause.volume

#Impulse response functions
##adjusted close->adjusted close
ADJ.C.irf <- irf(Model1, impulse = "Adj.close", response = "Adj.close", n.ahead = 40, boot = TRUE)
plot(ADJ.C.irf, ylab = "Adj.close", main = "Adj.close shock to Adj.close")
##volume->adjusted close
v.irf <- irf(Model1, impulse = "volume", response = "Adj.close", n.ahead = 40, boot = TRUE)
plot(v.irf, ylab = "volume", main = "volume shock to Adj.close")

#Variance decomposition
bvs.vardec <- fevd(Model1, n.ahead = 10)
plot(bvs.vardec)

#Adding contemporaneous effects using a structural VAR
# Estimate structural coefficients
a <- diag(1, 2)
a[lower.tri(a)] <- NA

#Contemporaneous effects
svar_est <- SVAR(Model1, Amat = a, max.iter = 10000)
svar_est #coefficients
svar_est$Ase #standard errors

#IRFs and variance decomposition for SVAR
##IRF
ADJ.C.irf <- irf(svar_est, impulse = "Adj.close", response = "Adj.close", n.ahead = 40, boot = TRUE)
plot(ADJ.C.irf, ylab = "Adj.close", main = "Adj.close shock to Adj.close")
v.irf <- irf(svar_est, impulse = "volume", response = "Adj.close", n.ahead = 40, boot = TRUE)
plot(v.irf, ylab = "volume", main = "volume shock to Adj.close")
##Variance decomposition
bvs.vardec <- fevd(svar_est, n.ahead = 10)
plot(bvs.vardec)

#Modeling volatility in data
##plot of log returns of stock prices
plot.ts(d.close_vals, col = "bLUE4", ylab = "Close price returns (differenced close price)")

##Testing for the ARCH effects formally
##ARIMA (0,1,1) with a drift
m1 <- arima(d.close_vals, order=c(0,0,1), include.mean = T)

##ACF and PACF of Squared returns
sqreturns <- d.close_vals^2
acf(sqreturns, ylim = c(-0.01,0.15))
pacf(sqreturns, ylim = c(-0.02,0.15))

##ARCH heteroskedasticity test for residuals
arch.test(m1, output=TRUE)

#Fit an ARCH model
model1 <- garchFit( ~arma(0,1)+garch(1,0), d.close_vals, include.mean=T, trace=FALSE)
model1 
summary(model1)

#Fit a GARCH model
model2 <- garchFit( ~arma(0,1)+garch(1,1), d.close_vals, include.mean=T, trace=FALSE)
model2 
summary(model2) 

#predictions for GARCH model
prediction <- predict(model2, n.ahead=50, plot=TRUE)

##Extension for the GARCH model
#GJR GARCH model
gjrgarch1 <- garchFit(~arma(0,1)+garch(1,1), leverage=T,d.close_vals,trace=F,include.mean=T)
summary(gjrgarch1)

#Exponential GARCH model
egarch1 <- ugarchfit(ugarchspec(mean.model=list(armaOrder=c(0,1),include.mean=T),variance.model=list(model="eGARCH",garchOrder=c(1,1))),d.close_vals)
egarch1 

#predictions for Exponenttial GARCH model
forecast <- ugarchforecast(egarch1,n.ahead=50)
plot(forecast)
1