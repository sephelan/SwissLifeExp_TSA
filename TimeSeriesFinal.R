
# setwd("~/Desktop/Time Series")
#Dependencies
library(MASS)
library(TSA)
library(tseries)
library(ggResidpanel)
library(lmtest)

#read in data & data cleaning
par(mfrow=c(1,1))
life = read.csv("SwissLifeExpectancy.csv",header = TRUE)
life.ts.trunc = ts((life$life.expectancy[28:99]),start=life$year[28],end=life$year[99])
plot ( life.ts.trunc,ylab = 'Life Expectancy (Yrs.)')
#idenfy possible transformation
BoxCox.ar(life.ts.trunc)
#test for stationarity
adf.test(life.ts.trunc)
life.diff.trunc = diff(life.ts.trunc)
plot(life.diff.trunc,ylab='Difference in Life Expectancy (Yrs.)')
adf.test(life.diff.trunc)

#Model Specification & diagnostics
par(mfrow=c(1,2))
acf(life.diff.trunc)
pacf(life.diff.trunc)
eacf(life.diff.trunc)
res=armasubsets(y=life.diff.trunc,nar=7,nma=7,
                y.name='test', ar.method='ols')
par(mfrow=c(1,1))

plot(res)
#AR(4), AR(5), and ARMA(5,1) are all suggested above

#Parameter estimation
par(mfrow=c(1,2))

#ar4 model only has phi1&4 significant
ar4.diff.trunc = arima(life.diff.trunc,order = c(4,0,0),method='ML')
ar4.diff.trunc #check for significant coefficients, std err, & AIC
BIC(ar4.diff.trunc) 
plot(ar4.diff.trunc$residuals,ylab='AR(4) Residuals')
qqnorm(rstandard(ar4.diff.trunc)) #normal
qqline(rstandard(ar4.diff.trunc)) #normal
runs(rstandard(ar4.diff.trunc)) #independence
shapiro.test(rstandard(ar4.diff.trunc)) #normality
par(mfrow=c(1,1))
tsdiag(ar4.diff.trunc,gof=20,omit.initial=F) #Box-ljung & resid. acf
par(mfrow=c(1,2))

#ar5 component is not significant -> no overfitting
ar5.diff.trunc = arima(life.diff.trunc,order = c(5,0,0),method='ML')
ar5.diff.trunc



par(mfrow=c(1,2))
#ARMA41 has an insignificant rho1, meaning ma component not valid.
arma41.diff.trunc = arima(life.diff.trunc,order = c(4,0,1),method='ML')
arma41.diff.trunc
BIC(ar4.diff.trunc)
plot(ar4.diff.trunc$residuals,ylab='ARMA(4,1) Residuals')
qqnorm(rstandard(ar4.diff.trunc))
qqline(rstandard(ar4.diff.trunc))
runs(rstandard(ar4.diff.trunc))
shapiro.test(rstandard(ar4.diff.trunc))
par(mfrow=c(1,1))
tsdiag(ar4.diff.trunc,gof=20,omit.initial=F)
par(mfrow=c(1,2))

#AR4 with only phi1&4 outperforms ar4 full model
par(mfrow=c(1,2))
ar4.reduced = arima(life.diff.trunc,order=c(4,0,0),include.mean = TRUE,fixed=c(NA,0,0,NA,NA),method='ML')
ar4.reduced
BIC(ar4.reduced)
plot(ar4.reduced$residuals,ylab='Reduced AR(4) Residuals')
qqnorm(rstandard(ar4.reduced))
qqline(rstandard(ar4.reduced))
runs(rstandard(ar4.reduced))
shapiro.test(rstandard(ar4.reduced))
par(mfrow=c(1,1))

tsdiag(ar4.reduced,gof=20,omit.initial=F)
par(mfrow=c(1,2))
#Becuase onlt the ar1&4 coefficients are significant, we test for overfitting with an ar1 model.
#ar1 model performs similar to the reduced ar4 model, but with fatter tails in the residuals.
ar1.diff.trunc = arima(life.diff.trunc,order = c(1,0,0),method='ML')
ar1.diff.trunc
BIC(ar1.diff.trunc)
plot(ar1.diff.trunc$residuals,ylab='AR(1) Residuals')
qqnorm(rstandard(ar1.diff.trunc))
qqline(rstandard(ar1.diff.trunc))
runs(rstandard(ar1.diff.trunc))
shapiro.test(rstandard(ar1.diff.trunc))
par(mfrow=c(1,1))
  tsdiag(ar1.diff.trunc,gof=20,omit.initial=F)

#We see the ar1 and reduced ar4 models are most appropriate for swiss life expectancy

#Forecasting
ar4.pred = predict(ar4.reduced,n.ahead=8)
ar4.real.pred = cumsum(c(life.ts.trunc[72],ar4.pred$pred[1:7]))

ar1.invert = c(life.ts.trunc[1:72])
plot(ts(ar1.invert,start=life$year[28],end=life$year[99]),ylab="Life Expectancy",xlim=c(1950,2028),ylim=c(68,90),main="ARIMA(4,1,0) predictions")
points(c(2022:2028),ar4.real.pred[2:8],pch=1)

lines(y=ar4.real.pred[2:8]+1.96*cumsum(ar4.pred$se[1:7]),x=c(2022:2028),lwd=2,col="red",lty="dashed")
lines(y=ar4.real.pred[2:8]-1.96*cumsum(ar4.pred$se[1:7]),x=c(2022:2028),lwd=2,col="red",lty="dashed")


ar4.pred.interval=data.frame("lower"=c(ar4.real.pred[2:8]-1.96*cumsum(ar4.pred$se[1:7])),"upper"=c(ar4.real.pred[2:8]+1.96*cumsum(ar4.pred$se[1:7])))
ar4.pred.interval <- cbind(ar4.pred.interval,data.frame('predicted'=ar4.real.pred[2:8]))
ar4.pred.interval <- cbind(ar4.pred.interval,data.frame('actual'=c(83.5,83.9,NA,NA,NA,NA,NA))) 
# 2022&2023 actual comes from OECD data explore (public health -> life expectancy -> switzerland -> all genders)
ar4.pred.interval

#we see the actual fall within our prediction interval
# we see the ar1 seems like a less noisy ar4 model. it also outperforms ar4 on the predictions which are closer to our observed data.

ar1.pred = predict(ar1.diff.trunc,n.ahead = 8)
ar1.pred
ar1.invert = c(life.ts.trunc[1:72])
plot(ts(ar1.invert,start=life$year[28],end=life$year[99]),ylab="Life Expectancy",xlim=c(1950,2028),ylim=c(68,90),main="ARIMA(1,1,0) predictions")
ar1.real.pred = cumsum(c(life.ts.trunc[72],ar1.pred$pred[1:7]))
points(c(2022:2028),ar1.real.pred[2:8],pch=1)

lines(y=ar1.real.pred[2:8]+1.96*cumsum(ar1.pred$se[1:7]),x=c(2022:2028),lwd=2,col="red",lty="dashed")
lines(y=ar1.real.pred[2:8]-1.96*cumsum(ar1.pred$se[1:7]),x=c(2022:2028),lwd=2,col="red",lty="dashed")

plot(ar1.diff.trunc,main='AR(1) Difference Predictions')

ar1.pred.interval=data.frame("lower"=c(ar1.real.pred[2:8]-1.96*cumsum(ar1.pred$se[1:7])),"upper"=c(ar1.real.pred[2:8]+1.96*cumsum(ar1.pred$se[1:7])))
ar1.pred.interval <- cbind(ar1.pred.interval,data.frame('predicted'=ar1.real.pred[2:8]))
ar1.pred.interval <- cbind(ar1.pred.interval,data.frame('actual'=c(83.5,83.9,NA,NA,NA,NA,NA))) 
ar1.pred.interval
