N=200
p = .8
nt = ceiling(p*N)
nc = N - nt
treatind = c(rep(T, nt), rep(F, nc))
S = 100
B = 100
rothe_calLOO_Centered <- rothe_LOO_Centered <- taudimCentered <- taugobCentered <- taucal1Centered <-taucal2Centered <- taubarCentered <- matrix(0, nrow = S, ncol = B) # Each row corresponds to a single realization of covariates and outcomes

PATE_Estimated = matrix(0, nrow = S, ncol = B) # Each row corresponds to a single realization of covariates and outcomes

MSEtaucal1 = rep(0, S)
MSEtaucal2 = rep(0, S)
MSEtaudim = rep(0, S)
MSEtaugob = rep(0, S)
MSErothe_LOO = rep(0, S)
MSErothe_calLOO = rep(0, S)
# vLOOP = rep(0, S)

f = function(x)
{
  return(exp(x))
}
g = function(x)
{
  return(1.8*(40 - .25*exp(x)))
}

# Numerically integrate for the population average treatment effect
PATE = integrate(f, -5, 5)$value/10 - ((integrate(g,-5, 5)$value))/10 # implicitly marginalizes over all realizations of potential outcomes (handy property of Poisson random variables)

# Computes:
#   Rothe's leave-one-out Poisson regression estimator (see http://www.christophrothe.net/papers/fca_apr2020.pdf) - Based on the paper "Flexible Covariate Adjustments in Randomized Experiments" by Christoph Rothe
#   A calibrated version of  Rothe's leave-one-out Poisson regression estimator
#   The LOOP Estimator of Gagnon-Barch & Wu (https://arxiv.org/pdf/1708.01229.pdf)
LeaveOneOut = function(observations, covariates, allocation)
{
  # U = (N/nt)*allocation - (N/nc)*(1 - allocation) # See https://arxiv.org/pdf/1708.01229.pdf equation 7 for definition
  # LOOPVec = numeric(N)

  psiCalVec  <- psiVec <- numeric(N)
  for(i in 1:N)
  {
    # Hold out the ith data point
    Obs_holdOuti = observations[-i]
    Cov_holdOuti = covariates[-i]
    treat_holdOuti = allocation[-i]
    data = data.frame(Y = Obs_holdOuti, X = Cov_holdOuti, treat = treat_holdOuti)

    # Fit on the data set without the ith data point
    mu1 = glm(Y ~ .-treat, family=poisson, subset(data, treat==1))
    mu0 = glm(Y ~ .-treat, family=poisson, subset(data, treat==0))
    tiHat = predict(mu1, data.frame(Y = observations[i],
                                    X = covariates[i],
                                    treat = allocation[i]), "r")
    ciHat = predict(mu0, data.frame(Y = observations[i],
                                    X = covariates[i],
                                    treat = allocation[i]), "r")

    # # Gagnon-Barch & Wu's LOOP Estimator
    # mihat = (1 - (nt/N))*tiHat + (nt/N)*ciHat
    # LOOPVec[i] = (observations[i] - mihat)*U[i]

    # The uncalibrated Poisson Rothe Leave-One-Out estimator
    # See the equation on the bottom of pg. 5 of http://www.christophrothe.net/papers/fca_apr2020.pdf
    psiVec[i] = (tiHat - ciHat) +
      allocation[i]*(N/nt)*(observations[i] - tiHat) -
      (1 - allocation[i])*(N/nc)*(observations[i] - ciHat)

    # Impute for all but the ith data point
    fittedT = predict(mu1, data, "r")
    fittedC = predict(mu0, data, "r")

    # Calibration with both predicted outcomes
    Xtilde = data.frame(fittedT, fittedC)
    XtildeT = Xtilde[treat_holdOuti,]
    XtildeC = Xtilde[!treat_holdOuti,]
    Ytreat_minusi = Obs_holdOuti[treat_holdOuti]
    Ycontrol_minusi = Obs_holdOuti[!treat_holdOuti]
    lmT2 = lm(Ytreat_minusi~., data = XtildeT)
    lmC2 = lm(Ycontrol_minusi~., data = XtildeC)


    ithPoint = data.frame(fittedT = tiHat, fittedC = ciHat)
    fittedT2 = predict(lmT2, newdata = ithPoint)
    fittedC2 = predict(lmC2, newdata = ithPoint)

    # See section 3.4 of http://www.christophrothe.net/papers/fca_apr2020.pdf
    psiCalVec[i] = (fittedT2 - fittedC2) +
      allocation[i]*(N/nt)*(observations[i] - fittedT2) -
      (1 - allocation[i])*(N/nc)*(observations[i] - fittedC2)
  }
  # LOOPEstimate = as.numeric(mean(LOOPVec))
  rothLOOEstimate = mean(psiVec)
  rothLOOEstimate_Calibrated = mean(psiCalVec)

  # return(c(rothLOOEstimate, rothLOOEstimate_Calibrated, LOOPEstimate))
  return(c(rothLOOEstimate, rothLOOEstimate_Calibrated))

}

for(i in 1:S)
{
  for (j in 1:B)
  {
    X = runif(N, -5, 5)
    Xcenter = cbind(X)

    YcTemp= rpois(N, g(X))
    YtTemp = rpois(N, f(X))
    Yt = YtTemp
    Yc = YcTemp

    PATE_Estimated[i,j] = mean(Yt - Yc)

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
    taudimCentered[i, j] = mean(Ytreat) - mean(Ycontrol) - PATE

    data = data.frame(Y = Y, X = X, treat = treat)
    mu1 = glm(Y ~ .-treat, family=poisson, subset(data, treat==1))
    mu0 = glm(Y ~ .-treat, family=poisson, subset(data, treat==0))
    fittedT = predict(mu1, data, "r")
    fittedC = predict(mu0, data, "r")
    # taugob[i, j] = mean(fittedT - fittedC) # The uncalibrated difference in means


    # No calibration
    yti = fittedT*(!treat) + Y*treat
    yci = fittedC*(treat) + Y*(!treat)
    taugobCentered[i, j] =  mean(yti-yci) - PATE
    # vhatgob[i] = var((fittedT-Y)[treat])/nt + var((fittedC-Y)[!treat])/nc + var(fittedT-fittedC)/(nc+nt) # This is an asymptotic approximation (since really this is a random variable based upon the sample mu0 and mu1)

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
    taucal2Centered[i, j] = mean(yti-yci) - PATE
    # vhatcal2[i] = var((fittedT2-Y)[treat])/nt + var((fittedC2-Y)[!treat])/nc+var(fittedT2-fittedC2)/(nc+nt) # This is an asymptotic approximation (since really this is a random variable based upon the sample mu0 and mu1)

    # Calibration using only one predicted value (without \hat{\mu}_{1-treat})
    fittedT = data.frame(fittedT)
    fittedTT = fittedT[treat,,drop = F]
    fittedC = data.frame(fittedC)
    fittedCC = fittedC[!treat,,drop = F]
    lmT3 = lm(Ytreat~., data = fittedTT)
    lmC3 = lm(Ycontrol~., data = fittedCC)
    fittedT3 = predict(lmT3, newdata = fittedT)
    fittedC3 = predict(lmC3, newdata = fittedC)
    taucal1Centered[i, j] = mean(fittedT3-fittedC3) - PATE
    # vhatcal1[i] = var((fittedT3-Y)[treat])/nt + var((fittedC3-Y)[!treat])/nc + var(fittedT3-fittedC3)/(nc+nt)# This is an asymptotic approximation (since really this is a random variable based upon the sample mu0 and mu1)

    # Compute Rothe LOO, and Rothe LOO Calibrated
    leaveOneOutResults = LeaveOneOut(observations = Y, covariates = X, allocation = treat)

    # Rothe's LOO estimator (see https://arxiv.org/pdf/1708.01229.pdf for details)
    rothe_LOO_Centered[i, j] = leaveOneOutResults[1] - PATE

    # Calibrating Rothe's LOO estimator
    rothe_calLOO_Centered[i, j] = leaveOneOutResults[2] - PATE

    # # LOOP estimator (see https://arxiv.org/pdf/1708.01229.pdf for details)
    # LOOPCentered[i, j] = leaveOneOutResults[3] - PATE
  }
  
  # Variances
#  vtaucal1[i] = var(taucal1Centered[i, ]) #Guo and Basse's singly calibrated
#  vtaucal2[i] = var(taucal2Centered[i, ]) #our fully calibrated estimator
#  vtaudim[i] = var(taudimCentered[i, ]) #unadjusted difference in means
#  vtaugob[i] = var(taugobCentered[i, ]) #Guo and Basse's uncalibrated Oaxaca-Blinder estimator
#  vrothe_LOO[i] = var(rothe_LOO_Centered[i, ]) #Rothe's uncalibrated LOO estimator
#  vrothe_calLOO[i] = var(rothe_calLOO_Centered[i, ]) #Rothe's calibrated LOO estimator
#  vLOOP[i] = var(LOOPCentered[i, ]) # LOOP estimator

  # MSEs
  MSEtaucal1[i] = mean(taucal1Centered[i, ]^2) #Guo and Basse's singly calibrated
  MSEtaucal2[i] = mean(taucal2Centered[i, ]^2) #our fully calibrated estimator
  MSEtaudim[i] = mean(taudimCentered[i, ]^2) #unadjusted difference in means
  MSEtaugob[i] = mean(taugobCentered[i, ]^2) #Guo and Basse's uncalibrated Oaxaca-Blinder estimator
  MSErothe_LOO[i] = mean(rothe_LOO_Centered[i, ]^2) #Rothe's uncalibrated LOO estimator
  MSErothe_calLOO[i] = mean(rothe_calLOO_Centered[i, ]^2) #Rothe's calibrated LOO estimator

}

# Notice that LOOP & Rothe (without calibration) are the same.
# all(abs(LOOPCentered - rothe_LOO_Centered) < 1E-10)

#PATE Inference
goB_vs_dim = mean(MSEtaugob / MSEtaudim)
cal1_vs_dim = mean(MSEtaucal1 / MSEtaudim)
cal2_vs_dim = mean(MSEtaucal2 / MSEtaudim)
LOO_vs_dim = mean(MSErothe_LOO / MSEtaudim)
calLOO_vs_dim = mean(MSErothe_calLOO / MSEtaudim)
# LOOP_vs_dim = mean(MSELOOP / MSEtaudim)

cat("PATE\n")

cat('gOB vs DiM\t\t', goB_vs_dim, '\n')
cat('cal1 vs DiM\t\t', cal1_vs_dim, '\n')
cat('cal2 vs DiM\t\t', cal2_vs_dim, '\n')
cat('Rothe LOO vs DiM\t', LOO_vs_dim, '\n')
cat('Cal. Rothe LOO vs DiM\t', calLOO_vs_dim, '\n')
# cat('LOOP vs DiM\t\t', LOOP_vs_dim, '\n')

vars = data.frame(goB_vs_dim, cal1_vs_dim, cal2_vs_dim, LOO_vs_dim, calLOO_vs_dim) #, LOOP_vs_dim)

varsName = paste("Poisson_PATE_MSERatios_withLeaveOneOuts_N", N, ".csv", sep = '')

write.csv(vars, file = varsName, row.names = FALSE)
