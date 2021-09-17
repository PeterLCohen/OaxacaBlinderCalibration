N=200
p = .8
nt = round(p*N)
nc = N - nt
treatind = c(rep(T, nt), rep(F, nc))
S = 100
B = 100
taudimCentered <- taugobCentered <- taucal1Centered <-taucal2Centered <- taubarCentered <- matrix(0, nrow = S, ncol = B) # Each row corresponds to a single realization of covariates
SATE <- rep(0, S)
vtaucal1 = rep(0, S)
vtaucal2 = rep(0, S)
vtaudim = rep(0, S)
vtaugob = rep(0, S)

mult = 1*(p-1)/p

f1 = function(x)
{
  return(exp(-3 + 2*x)/(1+exp(-3 + 2*x)))
}
f = Vectorize(f1)


PATE = integrate(f, -8, 8)$value/16 - (mult*1.6*(integrate(f,-8, 8)$value/16 -1))

for(i in 1:S)
{
  X = runif(N, -8, 8)
  Xcenter = cbind(X)
  pc= mult*(1.6)*(f(X)-1)
  pt= f(X)
  
  Yc = rbinom(N, 1, pc)
  Yt = rbinom(N, 1, pt)
  SATE[i] = mean(Yt-Yc)
  for (j in 1:B)
  {
    treat = sample(treatind)
    
    Ytreat = Yt[treat]
    Ycontrol = Yc[!treat]
    Xtreat = Xcenter[treat,,drop = F]
    Xcontrol = Xcenter[!treat,, drop = F]
    colnames(Xcontrol) = c("X1")
    colnames(Xtreat) = colnames(Xcontrol)
    colnames(Xcenter) = colnames(Xcontrol)
    Xcenter = data.frame(Xcenter)
    Xtreat = data.frame(Xtreat)
    Xcontrol = data.frame(Xcontrol)
    Y = Yt*treat + Yc*(1-treat)
    
    # DiM
    taudimCentered[i, j] = mean(Ytreat) - mean(Ycontrol) - SATE[i]
    
    data = data.frame(Y = Y, X = X, treat = treat)
    mu1 = glm(Y ~ .-treat, family=binomial, subset(data, treat==1))
    mu0 = glm(Y ~ .-treat, family=binomial, subset(data, treat==0))
    fittedT = predict(mu1, data, "r")
    fittedC = predict(mu0, data, "r")
    # taugob[i] = mean(fittedT - fittedC) # The uncalibrated difference in means
    
    # No calibration
    yti = fittedT*(!treat) + Y*treat
    yci = fittedC*(treat) + Y*(!treat)
    taugobCentered[i, j] =  mean(yti-yci) - SATE[i]
    
    # Calibration with both predicted outcomes
    Xtilde = data.frame(fittedT, fittedC)
    XtildeT = Xtilde[treat,]
    XtildeC = Xtilde[!treat,]
    lmT2 = lm(Ytreat~., data = XtildeT)
    lmC2 = lm(Ycontrol~., data = XtildeC)
    fittedT2 = predict(lmT2, newdata = Xtilde)
    fittedC2 = predict(lmC2, newdata = Xtilde)
    
    yti = fittedT2*(!treat) + Y*treat
    yci = fittedC2*(treat) + Y*(!treat)
    taucal2Centered[i, j] = mean(yti-yci) - SATE[i]
    
    # Calibration using only one predicted value (without \hat{\mu}_{1-treat})
    fittedT = data.frame(fittedT)
    fittedTT = fittedT[treat,,drop = F]
    fittedC = data.frame(fittedC)
    fittedCC = fittedC[!treat,,drop = F]
    lmT3 = lm(Ytreat~., data = fittedTT)
    lmC3 = lm(Ycontrol~., data = fittedCC)
    fittedT3 = predict(lmT3, newdata = fittedT)
    fittedC3 = predict(lmC3, newdata = fittedC)
    taucal1Centered[i, j] = mean(fittedT3-fittedC3) - SATE[i]
  }
  vtaucal1[i] = var(taucal1Centered[i, ]) #Guo and Basse's singly calibrated
  vtaucal2[i] = var(taucal2Centered[i, ]) #our fully calibrated estimator
  vtaudim[i] = var(taudimCentered[i, ]) #unadjusted difference in means
  vtaugob[i] = var(taugobCentered[i, ]) #Guo and Basse's uncalibrated Oaxaca-Blinder estimator
}


#SATE Inference
goB_vs_dim = mean(vtaugob / vtaudim)
cal1_vs_dim = mean(vtaucal1 / vtaudim)
cal2_vs_dim = mean(vtaucal2 / vtaudim)

cat("SATE\n")

#print(goB_vs_dim)
#print(cal1_vs_dim)
#print(cal2_vs_dim)

vars = data.frame(goB_vs_dim, cal1_vs_dim, cal2_vs_dim)

varsName = paste("Logistic_SATE_VarRatios_N", N, ".csv", sep = '')

write.csv(vars, file = varsName, row.names = FALSE)
