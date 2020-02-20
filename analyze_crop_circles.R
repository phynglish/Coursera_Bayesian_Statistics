setwd("D:/Users/phyng/Coursera/Bayesian Statistics")
df_init = read.csv(file = "crop_circle.csv", header = TRUE)
df_cropcircles = as.data.frame(table(df_init$Country, df_init$Year))
colnames(df_cropcircles) <- c("Country","Year","Number")

levels(df_cropcircles$Country)
summary(df_cropcircles)
pairs(df_cropcircles)
hist(df_cropcircles$Number)
boxplot(df_cropcircles$Number ~ df_cropcircles$Country, ylab = "Number of crop circles")
boxplot(df_cropcircles$Number ~ df_cropcircles$Year, ylab = "Number of crop circles")

library(rjags)
set.seed(113)

##### Model 1

mod_string1 = " model {
for (i in 1:length(Number)) {
  Number[i] ~ dpois(lam)
}

lam ~ dgamma(alpha,beta)

alpha = mu^2 / sig^2
beta = mu / sig^2

mu ~ dgamma(20.0, 1.0/1.0e6)
sig ~ dexp(1.0/10.0)

} "

data_jags1 = as.list(df_cropcircles)

params1 = c("lam", "mu", "sig")

mod1 = jags.model(textConnection(mod_string1), data=data_jags1, n.chains=3)
update(mod1, 5e3)

mod_sim1 = coda.samples(model=mod1,
                       variable.names=params1,
                       n.iter=1e4)
mod_csim1 = as.mcmc(do.call(rbind, mod_sim1))

## convergence diagnostics
plot(mod_sim1)

gelman.diag(mod_sim1)
autocorr.diag(mod_csim1)
autocorr.plot(mod_csim1)
effectiveSize(mod_csim1)

## compute DIC
(dic1 = dic.samples(mod1, n.iter=5e3))
summary(mod_csim1)


##### Model 2

mod_string2 = " model {
for (i in 1:length(Number)) {
  Number[i] ~ dpois(lam[Country[i]])
}

for (j in 1:max(Country)) {
  lam[j] ~ dgamma(alpha, beta)
}

alpha = mu^2 / sig^2
beta = mu / sig^2

mu ~ dgamma(20.0, 1.0/1.0e6)
sig ~ dexp(1.0/10.0)

} "

data_jags2 = as.list(df_cropcircles)

params2 = c("lam", "mu", "sig")

mod2 = jags.model(textConnection(mod_string2), data=data_jags2, n.chains=3)
update(mod2, 5e3)

mod_sim2 = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=1e4)
mod_csim2 = as.mcmc(do.call(rbind, mod_sim2))

## convergence diagnostics
plot(mod_sim2, ask = TRUE)

gelman.diag(mod_sim2)
autocorr.diag(mod_csim2)
autocorr.plot(mod_csim2)
effectiveSize(mod_csim2)

## compute DIC
(dic2 = dic.samples(mod2, n.iter=5e3))
summary(mod_csim2)


##### Model 3

mod_string3 = " model {
for (i in 1:length(Number)) {
  Number[i] ~ dpois(lam[i])
  log(lam[i]) = const[Country[i]] + dt[Country[i]]*Year[i]
}

for (j in 1:max(Country)) {
  const[j] ~ dnorm(0.0, 1.0/1.0e6)
  dt[j] ~ dnorm(0.0, 1.0/1.0e6)
}

} "

data_jags3 = as.list(df_cropcircles)
data_jags3$Year = as.numeric(data_jags3$Year)
head(data_jags3)
params3 = c("const", "dt")

mod3 = jags.model(textConnection(mod_string3), data=data_jags3, n.chains=3)
update(mod3, 5e3)

mod_sim3 = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=1e4)
mod_csim3 = as.mcmc(do.call(rbind, mod_sim3))

## convergence diagnostics
plot(mod_sim3, ask = TRUE)

gelman.diag(mod_sim3)
autocorr.diag(mod_csim3)
autocorr.plot(mod_csim3)
effectiveSize(mod_csim3)

## compute DIC
(dic1 = dic.samples(mod1, n.iter=5e3))
(dic2 = dic.samples(mod2, n.iter=5e3))
(dic3 = dic.samples(mod3, n.iter=5e3))
summary(mod_csim3)

library(plot3D)
scatter3D(as.numeric(df_cropcircles$Country),as.numeric(df_cropcircles$Year),df_cropcircles$Number)
as.numeric(df_cropcircles$Country)
library(ggplot2)
ggplot(df_cropcircles, aes(x = Year, y = Number, shape = factor(Country))) + geom_point()
with(df_cropcircles, plot(x = Year, y = Number, pch = Country))

summary(mod_csim3)

pm_params = colMeans(mod_csim3)
loglam = rep(pm_params[1:9], each = 7) + rep(pm_params[10:18], each = 7) * as.numeric(df_cropcircles[order(df_cropcircles$Country),]$Year)
yhat = exp(loglam)
resid = df_cropcircles[order(df_cropcircles$Country),]$Number - yhat

par(mfrow = c(1,2))
plot(resid)
plot(yhat, resid)
var(resid[yhat>30])


par(mfrow = c(3,3))
hist(mod_csim3[,10], xlab = expression(paste(Delta, "t"[1])), prob = TRUE, main = "Belgium")
hist(mod_csim3[,11], xlab = expression(paste(Delta, "t"[2])), prob = TRUE, main = "Canada")
hist(mod_csim3[,12], xlab = expression(paste(Delta, "t"[3])), prob = TRUE, main = "Czech")
hist(mod_csim3[,13], xlab = expression(paste(Delta, "t"[4])), prob = TRUE, main = "England")
hist(mod_csim3[,14], xlab = expression(paste(Delta, "t"[5])), prob = TRUE, main = "Germany")
hist(mod_csim3[,15], xlab = expression(paste(Delta, "t"[6])), prob = TRUE, main = "Netherlands")
hist(mod_csim3[,16], xlab = expression(paste(Delta, "t"[7])), prob = TRUE, main = "Italy")
hist(mod_csim3[,17], xlab = expression(paste(Delta, "t"[8])), prob = TRUE, main = "Switzerland")
hist(mod_csim3[,18], xlab = expression(paste(Delta, "t"[9])), prob = TRUE, main = "USA")
