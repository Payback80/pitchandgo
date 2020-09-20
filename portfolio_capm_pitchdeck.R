library(quantmod)
library(timeSeries)
library(stockPortfolio) # Base package for retrieving returns
library(ggplot2) # Used to graph efficient frontier
library(reshape2) # Used to melt the data
library(quadprog) #Needed for solve.QP
library(rjags)


#retrive data 
stocks <- c(
  "NSRGY" = .0,
  "RHHBY" = .0,
  "ASML" = .0,
  "SAP" = .0,
  "BTI" = .0,
  "SNY" = .0,
  "NVS" = .0,
  "AZN.L" = .0,
  
  "LVMUY" = .0,
  "ENIA" =.0,
  "E"    =.0
  )

pippo <- getSymbols(names(stocks), src = "yahoo", from = "2010-01-01", to = "2015-01-01")
#filling missing values
AZN.L <- na.approx(AZN.L)

Nestle <- NSRGY$NSRGY.Adjusted
Roche <- RHHBY$RHHBY.Adjusted
Asmlsv <- ASML$ASML.Adjusted
SapSA <- SAP$SAP.Adjusted
British <- BTI$BTI.Adjusted
Sanofi <- SNY$SNY.Adjusted
Novartis <- NVS$NVS.Adjusted
Astrazeneca <- AZN.L$AZN.L.Adjusted
Lvmh <- LVMUY$LVMUY.Adjusted
Enel <- ENIA$ENIA.Adjusted
Eni <-  E$E.Adjusted


assets <- merge(Nestle, Roche  )
assets <- merge(assets, Asmlsv)
assets <- merge(assets, SapSA)
assets <- merge(assets, British)
assets <- merge(assets, Sanofi)
assets <- merge(assets, Novartis)
assets <- merge(assets, Astrazeneca)
assets <- merge(assets, Lvmh)
assets <- merge(assets, Enel)
assets <- merge(assets, Eni)

#remove NA
assets <- na.omit(assets)



#pass the datafram as input

y = assets
n = dim(y)[1]
m = dim(y)[2]  
# m = dim(y)[2] - 1
return <- returns(y, percentage = TRUE) #calculate  percentage of daily returns  
return <- return[-c(1), ]
r = return
k1 = 250
k2 = k1 + 250
rtrain = r[1:k1, 1:m]
mkt_train = r[1:k1, 11]     #secund argument is the number of the assets involved in the calculation
rtest = r[(k1+1):k2, 1:m]
lambda = 2      #lambda range optimal between [2, 8] 10

data = list(R = rtrain, N = k1, mkt = mkt_train, m = m)
inits.Capm = function(){list(beta = rep(1,m))}
Capm.jags <- jags.model("BayesCapm.bug", data = data, inits = inits.Capm,
                        n.chains = 3, n.adapt = 1000, quiet = FALSE)
nthin = 10
N = 500
Capm.coda = coda.samples(Capm.jags, 
                         c("beta", "tauepsilon", "taubeta"), N * nthin, thin = nthin)
MCMC_out = Capm.coda[[1]]

summ = as.matrix(summary(Capm.coda)[[1]][,1])
summ
beta = summ[1:11] #10
taubeta = summ[12]
tauepsilon = summ[13:23]
sigmaepsilon = tauepsilon^(-.5)

pdf("two_betas.pdf", width = 6, height = 5)    ##  Figure 20.11
lmFit = lm(as.matrix(rtrain) ~ mkt_train)
beta_lm = lmFit$coeff[2,]
par(mfrow =  c(1,1))
plot(beta_lm, beta, xlim = c(0.7, 1.55), ylim = c(0.7, 1.55), 
     xlab = "beta - least squares", ylab = "beta - Bayes")
grid(lwd=2)
abline(0,1)
graphics.off()

ExUtil = function(w){
  -1 + exp(-lambda * (1 + t(w) %*% mu) + lambda^2 * t(w) %*% Omega %*% w /2 )
}

lambda = 2   #check
mu_model_free = colMeans(rtrain)
Omega_model_free = cov(rtrain)
mu_Capm = beta * mean(mkt_train)
Omega_Capm = beta %o% beta * var(mkt_train) + diag(sigmaepsilon^2)

mu = mu_model_free
Omega = Omega_model_free

library(quadprog)
opt1 = solve.QP(Dmat = as.matrix(lambda^2 * Omega), dvec = lambda * mu, 
                Amat = as.matrix(rep(1,11)), bvec = 1, meq = 1)
w_model_free = opt1$solution 

mu = mu_Capm
Omega = Omega_Capm
opt2 = solve.QP(Dmat = as.matrix(lambda^2 * Omega), dvec = lambda * mu, 
                Amat = as.matrix(rep(1,11)), bvec = 1, meq = 1)
w_Capm = opt2$solution

pdf("two_w_lambda3.pdf", width = 6, height = 5)    ##  Figure 20.10
plot(w_model_free, type = "b", col = "blue", main = paste("lambda = ", lambda), 
     ylab = "w", ylim = c(-2.25, 2.25), lwd = 2 )
points(w_Capm, type = "b", col = "red", pch = "*", lwd = 2)
legend("bottomleft", c("model-free", "CAPM"), lty = 1, 
       pch = c("o", "*"), col = c("blue", "red"), lwd = 2)
abline(h=0)
graphics.off()

return_model_free = as.matrix(rtest) %*% w_model_free
ExUt_model_free = mean(1 - exp(-lambda * return_model_free))

return_Capm = as.matrix(rtest) %*% w_Capm
ExUt_Capm = mean(1 - exp(-lambda * return_Capm))

print(c(ExUt_model_free, ExUt_Capm), digits = 2)
print(c(mean(return_model_free), mean(return_Capm)), digits = 2)
print(c(sd(return_model_free), sd(return_Capm)), digits = 2)
