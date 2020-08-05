data = read.csv("slash_survsex.csv")
head(data)
is.na(data)
n=length(data[,1])
mean(data$Female)
mean(data$SexualActivity)
table(data$Female, data$SexualActivity, dnn=c("female", "SexualActivity"))
mean(data[data[,1] == 1,]$SexualActivity)
mean(data[data[,1] == 0,]$SexualActivity)
mean(data[data[,2]== 1,]$Female)
mean(data$Survival)
mean(data[data[,1] == 1,]$Survival)
mean(data[data[,1] == 0,]$Survival)
colMeans(data)

pairs(sapply(data, jitter))

library("rjags")

mod1_string = " model {
    for (i in 1:length(Survival)) {
        Survival[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*Female[i] + b[2]*SexualActivity[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:2) {
        b[j] ~ dnorm(0.0, 1.0/25.0)
    }
} "

set.seed(42)

data_jags = as.list(data)
params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string), data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)
dev.print(pdf, "mod1.pdf")

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)

mod2_string = " model {
    for (i in 1:length(Survival)) {
        Survival[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*Female[i] + b[2]*SexualActivity[i] + b[3]*Female[i]*SexualActivity[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:3) {
        b[j] ~ dnorm(0.0, 1.0/25.0)
    }
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)
update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

## convergence diagnostics
plot(mod2_sim, ask=TRUE)
dev.print(pdf, "mod2.pdf")

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

## calculate DIC
dic2 = dic.samples(mod2, n.iter=1e3)

summary(mod2_csim)
HPDinterval(mod2_csim)
head(mod2_csim)
mean(mod2_csim[,1] > 0)
mean(mod2_csim[,2] < 0)
mean(mod2_csim[,3] < 0)

summary(mod1_csim)
HPDinterval(mod1_csim)
head(mod1_csim)
mean(mod1_csim[,1] > 0)
mean(mod1_csim[,2] < 0)

(pm_coef = colMeans(mod2_csim))
X = as.matrix(data[,-3])
X = cbind(X, with(data, Female*SexualActivity))
pm_Xb = pm_coef["int"] + X %*% pm_coef[1:3]
phat = 1.0 / (1.0 + exp(-pm_Xb))
head(phat)
plot(phat, jitter(data$Survival))
hist(phat)
dev.print(pdf, "predXlab.pdf")
(tab0.5 = table(phat > 0.5, data_jags$Survival))
sum(diag(tab0.5)) / sum(tab0.5)



(pm_coef1 = colMeans(mod1_csim))
pm_Xb1 = pm_coef1["int"] + X %*% pm_coef1[1:3]
phat1 = 1.0 / (1.0 + exp(-pm_Xb1))
head(phat1)
plot(phat1, jitter(data$Survival))
hist(phat1)
dev.print(pdf, "predXlab.pdf")
(tab0.5_1 = table(phat1 > 0.5, data_jags$Survival))
sum(diag(tab0.5_1)) / sum(tab0.5_1)

male_nosex = c(0, 0, 0)
male_sex = c(0,1,0)
female_nosex = c(1, 0,0)
female_sex = c(1,1,1)
X2 = rbind(male_nosex, male_sex, female_nosex, female_sex)
X2
1.0 / (1.0 + exp(-(pm_coef["int"] + X2 %*% pm_coef[1:3])))