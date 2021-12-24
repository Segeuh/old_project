#Different projects over the topics of the chapter 3

# Example OLS classic
# Example regression sur equation 2.3 log
# Autocorrelation
# Test de Chow


#### Example: OLS estimator######################################################################
rm(list=ls())

#install.packages("ggplot2")
library(ggplot2)


# generate variables y and x

set.seed(20)


n <- 100
k <- 2

eps <- rnorm(n,0,2)   # eps~N(0,2) normal distribution with mean zero and standard deviation two

x1 <- rep(1,n)
#x2 <- c(1:n)
x2 <- rnorm(n, 0, 1)  # alternative of stochastic regressor
beta <- c(1,0.5)

# define matrix X and vector y
X <- cbind(x1,x2)

y <- X%*% beta + eps

ggplot(data.frame(x=x2,y=y), aes(x=x2, y=y)) + geom_point(size=4, colour = "red", alpha=0.5) + 
  geom_abline(intercept=beta[1], slope=beta[2]) 

## Calculate OLS estimator

beta_hat <- solve(t(X) %*% X) %*% t(X) %*% y
fit1 <- lm(y ~ x2)      # function includes intercept by default, alternative lm(y ~ x1 + x2 -1)
summary(fit1)

ggplot(data.frame(x=x2,y=y), aes(x=x2, y=y)) + geom_point(size=4, colour = "red", alpha=0.5) + 
  geom_abline(intercept=beta[1], slope=beta[2]) + geom_abline(intercept=beta_hat[1], slope=beta_hat[2], colour = "blue") 


e <- y - X%*%beta_hat      # estimated residuals
# e <- residuals(fit1)

sigma_hat <- (t(e) %*% e) / (n-k)   # estimated variance
sd_hat <- sqrt(sigma_hat)           # estimated standard deviation
var_betahat <- drop(sigma_hat) * solve(t(X) %*% X)    # estimated variance of beta_hat, drop() for making sigma_hat a scalar
sd_beta1 <- sqrt(var_betahat[1,1])
sd_beta2 <- sqrt(var_betahat[2,2])


R2 <- 1 - (t(e) %*% e) / (t(y - mean(y)) %*% (y - mean(y))) 
R2_adj <-  1 - ((n-1)* t(e) %*% e) / ((n-k) * t(y - mean(y)) %*% (y - mean(y))) 



### Example regression equation 2.3 ##############################################"
rm(list=ls())

library(readxl)

lifeexpec <- read_excel("lifeexpec.xls")

# define variables
lifeexp <- lifeexpec$`Life expectancy at birth, total (years)`
gdp <- lifeexpec$`GDP per capita (current US$)`
lgdp <- log(gdp)
educ <- lifeexpec$`Public spending on education, total (% of GDP)`

# fit linear regression model: lifeexp = beta1 + beta2*lgdp + beta3*educ + eps
fit <- lm(lifeexp ~ lgdp + educ)
summary(fit)

SSE <- sum(residuals(fit)^2)
SST <- sum((lifeexp-mean(lifeexp))^2)

# fit linear regression on constant only: lifeexp = beta1 + eps
fit_red <- lm(lifeexp ~ 1)
summary(fit_red)

SSE_red <- sum(residuals(fit_red)^2)


### Split sample

fit_full <- lm(lifeexp ~ lgdp)

plot(lgdp, lifeexp)

lgdp_high <- lgdp[lgdp>9]
lgdp_low <- lgdp[lgdp<9]

lifeexp_high <- lifeexp[lgdp>9]
lifeexp_low <- lifeexp[lgdp<9]

fit_high <- lm(lifeexp_high ~ lgdp_high)
fit_low <- lm(lifeexp_low ~ lgdp_low)



abline(a=fit_high$coefficients[1], b=fit_high$coefficients[2], col="blue" )
abline(a=fit_low$coefficients[1], b=fit_low$coefficients[2], col="red" )
abline(a=fit_full$coefficients[1], b=fit_full$coefficients[2] )


SSE_high <- sum(residuals(fit_high)^2)
SSE_low <- sum(residuals(fit_low)^2)

GQ <-(SSE_low/(67-2))/(SSE_high/(42-2))

crit5 <- qf(0.95,(67-2),(42-2))
p_GQ <- 1 -pf(GQ,(67-2),(42-2))


### Mean and plot example :

rm(list=ls())


mu <- 2
n <- 200
eps <- rnorm(n,0,1)

Y <- 2 + eps


library(ggplot2)
ggplot(data.frame(x=1:n,y=Y), aes(x=1:n, y=y)) + geom_point(size=4, colour = "red", alpha=0.5) + 
  geom_abline(intercept=mu, slope=0) 

mu_hat1 <- sum(Y)/n


### Example autocorrelation #######################################"

rm(list=ls())

#install.packages("devtools")
devtools::install_git("https://github.com/ccolonescu/PoEdata")

data("fred", package="PoEdata")
fred <- ts(fred, start=c(1960,1),end=c(2009,4),frequency=4)

ts.plot(fred[,"c"],fred[,"y"], type="l", 
        lty=c(1,2), col=c(1,2))
legend("topleft", border=NULL, legend=c("c","y"), 
       lty=c(1,2), col=c(1,2))

c <- fred[,"c"]
y <- fred[,"y"]
T <- length(y)

plot( y, c)

mod1 <- lm(c ~ y)
summary(mod1)
abline(mod1)

library(lmtest)
dwtest(c ~ y)

resid <- mod1$residuals

plot(resid, type="l")
#plot(rnorm(200), type="l")

library(Hmisc)
resid_reg <- lm(resid ~ Lag(resid, +1) - 1) 
summary(resid_reg)
rho <- resid_reg$coefficients[1]


y_mod <- y -rho*Lag(y, +1)
c_mod <- c -rho*Lag(c, +1)

y0 <- rep(1-rho,T)
c_mod[1] <- sqrt(1-rho^2)*c[1]
y_mod[1] <- sqrt(1-rho^2)*y[1]
y0[1] <- sqrt(1-rho^2)

reg_mod <- lm( c_mod ~ y0 + y_mod -1)
summary(reg_mod)

plot(residuals(reg_mod), type="l")


library(sandwich)
coeftest(mod1, vcov = vcovHC(mod1, "HC1"))


### Example: Chow-test #############################################
rm(list=ls())

#install.packages("strucchange")
library(strucchange) 
data(USIncExp)

USIncExp2 <- window(USIncExp, start = c(1985,12))

fit1 <- lm(expenditure ~ income, data = USIncExp2)
summary(fit1)


## Potentially spurious effects
plot(USIncExp2[,"income"], USIncExp2[,"expenditure"])

## Regressing percentage changes instead
USIncExp2 <- cbind(USIncExp2,diff((USIncExp2[,"income"])),diff((USIncExp2[,"expenditure"])))
colnames(USIncExp2) <- c("income", "expenditure", "diff.income", "diff.expenditure")
plot(USIncExp2[,"diff.income"], USIncExp2[,"diff.expenditure"])
USIncExp2 <- window(USIncExp2, start = c(1986,1), end = c(2001,2))

fit2 <- lm(diff.expenditure ~ diff.income, data=USIncExp2)
summary(fit2)


fit3 <- lm(diff.expenditure[-100<diff.income & diff.income<100] ~ diff.income[-100<diff.income & diff.income<100], data=USIncExp2)
summary(fit3)

abline(fit2$coefficients[1], fit2$coefficients[2], col="blue")

plot(USIncExp2[,"diff.income"], USIncExp2[,"diff.expenditure"], xlim=c(-100,100))
abline(fit3$coefficients[1], fit3$coefficients[2], col="red")


fs <- Fstats(fit2, from = c(1990, 1), to = c(1999,6), data = USIncExp2)
plot(fs)



### Improvement: account for spurious relation by error correction model 

USIncExp2 <- window(USIncExp, start = c(1985,12))

coint.res <- residuals(lm(expenditure ~ income, data = USIncExp2))
coint.res <- lag(ts(coint.res, start = c(1985,12), freq = 12), k = -1)
USIncExp2 <- cbind(USIncExp2, diff(USIncExp2), coint.res)
USIncExp2 <- window(USIncExp2, start = c(1986,1), end = c(2001,2))
colnames(USIncExp2) <- c("income", "expenditure", "diff.income","diff.expenditure", "coint.res")
ecm.model <- diff.expenditure ~ coint.res + diff.income
fs <- Fstats(ecm.model, from = c(1990, 1), to = c(1999,6), data = USIncExp2)
summary(fs)

plot(fs)
