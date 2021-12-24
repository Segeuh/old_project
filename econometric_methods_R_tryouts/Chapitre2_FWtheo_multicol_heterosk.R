#Different topics of the chapter 2 in action : 

#Frisch-Waugh Theorem
#Multicolinearity
#Heteroscedasticity : correction using white robust SE + Breuch Pagan

### Example: Frisch-Waugh Theorem ###################################

rm(list=ls())

#install.packages("wooldridge")
library(wooldridge)


data("wage2")

wages <- wage2

lwage <- wages$lwage
wage <- wages$wage
educ <- wages$educ
exper <- wages$exper


wage_reg <- lm(lwage ~ educ + exper)
summary(wage_reg)

wage_educ <- lm(lwage ~ educ)
summary(wage_educ)
resid1 <- residuals(wage_educ)

exper_educ <- lm(exper ~ educ)
summary(exper_educ)
resid2 <- residuals(exper_educ)

resid_reg <- lm(resid1 ~ resid2)
summary(resid_reg)


### Example: Multicollinearity ########################################

set.seed(3)

N <- 150

x1 <- runif(N)
x2 <- x1 + 0.1*runif(N)
x3 <- runif(N)

y <- -1 + 0.3*x1 + 0.2*x3 + rnorm(150, 0, 0.04)


plot(x1,x2)
plot(x1,x3)

X <- cbind(x1,x2,x3)

cor(X)


fit <- lm(y ~ x1 +x2 +x3)
summary(fit)


fit_x3 <- lm(y ~  x3)
summary(fit_x3)


fit_true <- lm(y ~ x1 +x3)
summary(fit_true)

#install.packages("faraway")
library(faraway)

vif(fit)


fit_x2_x3 <- lm(y ~  x2 + x3)
summary(fit_x2_x3)

AIC(fit)
AIC(fit_x3)
AIC(fit_true)
AIC(fit_x2_x3)

### Example heteroskedasticity  ###########################################


rm(list=ls())

## obtain data from package PoEdata, downloaded from github
install.packages("devtools")
devtools::install_git("https://github.com/ccolonescu/PoEdata")

library(PoEdata) 
library(ggplot2)

data("food",package="PoEdata")
N <- 40

## simple linear regression model foodexp=beta1+ beta2*income+ eps

mod1 <- lm(food_exp ~ income, data=food)
summary(mod1)

ggplot(data = food, aes(y = food_exp, x = income)) + geom_point(col = 'blue') + 
  geom_abline(slope= mod1$coefficients[2], intercept=mod1$coefficients[1])

## Regression residual and plot
food$residuals <- mod1$residuals
food$fitted <- mod1$fitted.values
ggplot(data = food, aes(y = residuals, x = income)) + geom_point(col = 'blue') + geom_abline(slope = 0)
ggplot(data = food, aes(y = residuals, x = fitted)) + geom_point(col = 'blue') + geom_abline(slope = 0)

## Weighted least squares, known form of heteroskedasticity
foodexp_gls <- food$food_exp/sqrt(food$income)
inter_gls <- rep(1,40)/sqrt(food$income)
mod_gls <- lm(foodexp_gls ~ inter_gls + sqrt(food$income) -1)
summary(mod_gls)

## Residuals from GLS/WLS model
food$resid_gls <- mod_gls$residuals
food$fitted_gls <- mod_gls$fitted.values
ggplot(data = food, aes(y = resid_gls, x = fitted_gls)) + geom_point(col = 'blue') + geom_abline(slope = 0)

## Breusch Pagan test for heteroskedasticity
#install.packages("lmtest")
library(lmtest)
bptest(mod1)

## White test for heteroskedasticity
#install.packages("broom")
library(broom)
ressq <- resid(mod1)^2
modres <- lm(ressq~income+I(income^2), data=food)
gmodres <- glance(modres)
Rsq <- gmodres$r.squared
S <- gmodres$df #Number of regression coefficients in the model
chisq <- N*Rsq  #Statistic of white test
pval <- 1-pchisq(chisq, S-1)

## OLS estimation with robust standard errors
library(sandwich)
coeftest(mod1, vcov = vcovHC(mod1, "HC1"))

