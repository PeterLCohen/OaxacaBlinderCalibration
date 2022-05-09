library(devtools)
devtools::install_github("swager/randomForest")
library(randomForest)

# Based on code from swager/randomForest GitHub repo
ate.randomForest_calibrated = function(X, Y, W, nodesize = 20, conf.level=.9) {
  
  if (prod(W %in% c(0, 1)) != 1) {
    stop("Treatment assignment W must be encoded as 0-1 vector.")
  }
  
  nobs = nrow(X)
  pobs = ncol(X)
  
  if(length(unique(Y)) > 2) {
    # Sample-Split (two folds)
    tauVec = numeric(2)
    tauVec_cal = numeric(2)
    for(i in 1:2)
    {
      # Hold out half of the data points
      inds = (seq(1, N) %% 2 == (i - 1)) # select based on parity of the index 
      Y_leaveouti = Y[!inds]
      X_leaveouti = X[!inds,]
      W_leaveouti = W[!inds]
      
      initialFit0 = rep(NA, length(W_leaveouti))
      initialFit1 = rep(NA, length(W_leaveouti))
      
      yhat.0 = rep(NA, length(W_leaveouti))
      yhat.1 = rep(NA, length(W_leaveouti))
      
      rf_leaveouti.0 = randomForest::randomForest(X_leaveouti[W_leaveouti==0,], Y_leaveouti[W_leaveouti==0], nodesize = nodesize)
      rf_leaveouti.1 = randomForest::randomForest(X_leaveouti[W_leaveouti==1,], Y_leaveouti[W_leaveouti==1], nodesize = nodesize)
      
      # Form the initial predictions based upon the random forests (on the data points that aren't withheld)
      initialFit0[W_leaveouti==0] = stats::predict(rf_leaveouti.0)
      initialFit0[W_leaveouti==1] = stats::predict(rf_leaveouti.0, newdata = X_leaveouti[W_leaveouti==1,])
      initialFit1[W_leaveouti==1] = stats::predict(rf_leaveouti.1)
      initialFit1[W_leaveouti==0] = stats::predict(rf_leaveouti.1, newdata = X_leaveouti[W_leaveouti==0,])
      
      # Initial fit on the held out data
      tiHat = stats::predict(rf_leaveouti.1, newdata = X[inds,])
      ciHat = stats::predict(rf_leaveouti.1, newdata = X[inds,])
      
      # Form the calibrated predictions using the initial fits as covariates
      Xtilde = data.frame(initialFit1, initialFit0)
      Xtilde1 = Xtilde[W_leaveouti==1,]
      Xtilde0 = Xtilde[W_leaveouti==0,]
      lm1 = lm(Y_leaveouti[W_leaveouti==1]~., data = Xtilde1)
      lm0 = lm(Y_leaveouti[W_leaveouti==0]~., data = Xtilde0)
      calFit1 = predict(lm1, newdata = data.frame(initialFit1 = tiHat, initialFit0 = ciHat))
      calFit0 = predict(lm0, newdata = data.frame(initialFit1 = tiHat, initialFit0 = ciHat))
      
      
      tauVec[i] = mean((tiHat - ciHat) + W[inds]*(N/n1)*(Y[inds] - tiHat) - (1 - W[inds])*(N/n0)*(Y[inds] - ciHat))
      tauVec_cal[i] = mean((calFit1 - calFit0) + W[inds]*(N/n1)*(Y[inds] - calFit1) - (1 - W[inds])*(N/n0)*(Y[inds] - calFit0))
    }
  } else {
    cat("Error: Not enough unique values of Y")
  }
  
  list(uncalibrated = mean(tauVec), calibrated = mean(tauVec_cal))
}

set.seed(123)

library(MASS) # Required for the Boston data set
N = 2000
p = .6 # Hypothetical treatment portion
n1 = round(N*p)
n0 = N - n1
treatind = c(rep(T, n1), rep(F, n0))

nsims = 1000
uncalibrated = numeric(nsims)
calibrated = numeric(nsims)

for(i in 1:nsims)
{
  # Form a hypothetical realization of a CRE for which the outcome of interest is crime and the features are the remaining columns of the Boston data set
  # IID sample from the rows of the Boston data set to form a realization of a superpopulation experiment
  rowInds = sample(1:nrow(Boston), size = N, replace = TRUE)
  realization = Boston[rowInds,]
  crime = realization$crim
  Yt = crime + rexp(n = N)
  Yc = crime - rexp(n = N)
  features = realization[, -1]
  Z = sample(treatind) # Form a realization of the hypothetical treatment allocation
  outcomes = Yt*Z + Yc*(1 - Z)
  
  ate = ate.randomForest_calibrated(X = features, Y = outcomes, W = Z)
  uncalibrated[i] = ate$uncalibrated
  calibrated[i] = ate$calibrated
  cat("Iteration", i, "\n")
}

hist(sqrt(N)*uncalibrated, xlab = '', freq = F)
hist(sqrt(N)*calibrated, xlab = '', freq = F)

cat("N-scaled variance of uncalibrated estimator\t", N*var(uncalibrated), "\n")
cat("N-scaled variance of calibrated estimator\t", N*var(calibrated), "\n")

save(uncalibrated, file = paste("wagerUncalibrated_NeymanNull_SampleSplit_N", N, ".RData", sep=''))
save(calibrated, file = paste("wagerCalibrated_NeymanNull_SampleSplit_N", N, ".RData", sep=''))

dtf = data.frame(uncal = N*var(uncalibrated),cal = N*var(calibrated))
write.csv(dtf, file = paste("wager_NeymanNull_SampleSplit_N", N, ".csv", sep=''), row.names = F)
