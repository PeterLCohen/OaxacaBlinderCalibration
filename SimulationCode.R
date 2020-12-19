#logistic
library(mvtnorm)
library(sandwich)
options(warn=2) # Warnings as errors

set.seed(123)

nsim = 100
B = 100

N=200
nt = round(0.3*N)
nc = N - nt
treatind = c(rep(T, nt), rep(F, nc))

muhatOLS0_Record = data.frame(Intercept = numeric(0), fittedT = numeric(0), fittedC = numeric(0))
muhatOLS1_Record = data.frame(Intercept = numeric(0), fittedT = numeric(0), fittedC = numeric(0))

# Response functions
p_c = Vectorize(function(V)
{
  m = -.85
  s = 1.3
  return(1* dnorm(V, mean = m, sd = s) / max(dnorm(seq(from = -5, to = 5, length.out = 100), mean = m, sd = s)))
})

p_t = Vectorize(function(V){
  m = 1
  s = 1
  return(dnorm(V, mean = m, sd = s) / max(dnorm(seq(from = -5, to = 5, length.out = 100), mean = m, sd = s)))  
})

# Visualize the probabilities as a function of the covariates
x = seq(from = -5, to = 5, length.out = 100)
plot(x, p_c(x), ylim = c(0, 1), type = 'l', col ='red', main = 'Response Functions')
lines(x, p_t(x), ylim = c(0, 1), col = 'blue')
legend(x = 'topleft', col = c('red', 'blue'), legend = c('pc', 'pt'), lwd = c(2,2))

# Simulate over many data sets
Vdim <- Vgob <- Vcal1 <- Vcal2 <- rep(0, B)
for(b in 1:B)
{
  # Generate the covariates
  V = runif(N, -5, 5)
  X = cbind(V)
  
  # Generate underlying probabilities for treated and control groups
  pt = p_t(V)
  pc = p_c(V)
  
  
  # Generate potential outcomes
  U = runif(N)
  Yc = 1*(U<=pc)
  Yt = 1*(U<=pt)
  
  # The true difference in means
  taubar = mean(Yt-Yc)
  
  # Simulate many realizations of treatment allocation on a single data set
  modelFlag = rep(FALSE, nsim)
  taudim <- taugob <- taucal1 <-taucal2 <- rep(0, nsim)
  i = 1
  while(i <= nsim)
  {
    try({
      # Sample a CRE treatment allocation
      treat = sample(treatind)
      
      # What the experimenter gets to see
      Ytreat = Yt[treat]
      Ycontrol = Yc[!treat]
      Xtreat = X[treat,,drop = F]
      Xcontrol = X[!treat,, drop = F]
      colnames(Xcontrol) = c("X1")
      colnames(Xtreat) = colnames(Xcontrol)
      colnames(X) = colnames(Xcontrol)
      X = data.frame(X)
      Xtreat = data.frame(Xtreat)
      Xcontrol = data.frame(Xcontrol)
      Y = Yt*treat + Yc*(1-treat)
      
      # The standard difference in means (without any adjustment or imputation)
      taudim[i] = mean(Ytreat) - mean(Ycontrol)
      
      # "In-sample" logistic regressions
      lmT = glm(Ytreat~., data = Xtreat, family = "binomial")
      lmC = glm(Ycontrol~., data = Xcontrol, family = "binomial")
      fittedT = predict(lmT, newdata = X, type = "response")
      fittedC = predict(lmC, newdata = X, type = "response")
      
      # The generalized Oaxaca-Blinder estimator without any calibration
      yti = fittedT*(!treat) + Y*treat
      yci = fittedC*(treat) + Y*(!treat)
      taugob[i] =  mean(yti-yci)
      
      # Cohen & Fogarty calibration (regress on *both* imputed values)
      Xtilde = data.frame(fittedT, fittedC)
      XtildeT = Xtilde[treat,]
      XtildeC = Xtilde[!treat,]
      lmT2 = lm(Ytreat~., data = XtildeT)
      lmC2 = lm(Ycontrol~., data = XtildeC)
      fittedT2 = predict(lmT2, newdata = Xtilde)
      fittedC2 = predict(lmC2, newdata = Xtilde)
      taucal2[i] = mean(fittedT2-fittedC2)
      # yti = fittedT2*(!treat) + Y*treat
      # yci = fittedC2*(treat) + Y*(!treat)
      # taucal2[i] = mean(yti-yci)
      # Record the regressions
      muhatOLS0_Record = rbind(muhatOLS0_Record, lmC2$coefficients)
      muhatOLS1_Record = rbind(muhatOLS1_Record, lmT2$coefficients)
      
      # Calibration with only a single imputed value (i.e., Guo & Basse's Equation 8)
      # See https://arxiv.org/pdf/2004.11615.pdf (Version 1) for reference
      fittedT = data.frame(fittedT)
      fittedTT = fittedT[treat,,drop = F]
      fittedC = data.frame(fittedC)
      fittedCC = fittedC[!treat,,drop = F]
      lmT3 = lm(Ytreat~., data = fittedTT)
      lmC3 = lm(Ycontrol~., data = fittedCC)
      fittedT3 = predict(lmT3, newdata = fittedT)
      fittedC3 = predict(lmC3, newdata = fittedC)
      taucal1[i] = mean(fittedT3-fittedC3)
      i = i + 1
    })
  }
  
  # Record sample variances
  Vdim[b] = var(taudim - taubar)
  Vgob[b] = var(taugob - taubar)
  Vcal1[b]= var(taucal1 - taubar)
  Vcal2[b]= var(taucal2 - taubar)
  
  # print(b)
}

#cat("N = ", N , "\n")
#cat("Vgob / Vdim \t", mean(Vgob / Vdim), "\n")
#cat("Vcal1 / Vdim \t", mean(Vcal1 / Vdim), "\n")
#cat("Vcal2 / Vdim \t", mean(Vcal2 / Vdim), "\n")

colnames(muhatOLS0_Record) = c("Intercept", "FittedT", "FittedC")
colnames(muhatOLS1_Record) = c("Intercept", "FittedT", "FittedC")

cat("\nmuhatOLS0\nmeans\n")
print(round(colMeans(head(muhatOLS0_Record, n=nsim)), digits = 3))
cat("SDs\n")
print(round(apply(X = head(muhatOLS0_Record, n=nsim), 2, FUN = sd), digits = 3))

cat("\nmuhatOLS1\nmeans\n")
print(round(colMeans(head(muhatOLS1_Record, n=nsim)), digits = 3))
cat("SDs\n")
print(round(apply(X = head(muhatOLS1_Record, n=nsim), 2, FUN = sd), digits = 3))

vars = data.frame(Vcal2, Vcal1, Vgob, Vdim)
coeffs = data.frame(muhatOLS0_Record, muhatOLS1_Record)

varsName = paste("SamplingVariances_N", N, ".csv", sep = '')
coeffsName = paste("coefficients_N", N, ".csv", sep = '')

write.csv(vars, file = varsName, row.names = FALSE)
write.csv(coeffs, file = coeffsName, row.names = FALSE)



