#' Calibration of Oaxaca blinder estimators
#'
#' Implements calibration method of Cohen & Fogarty (2020) to build non-inferior
#'  Oaxaca-Blinder imputation estimators.
#'
#' @param muhat model to use in the treated and control groups (must be a family of GLMs)
#' @param Y the observed outcomes
#' @param Z the observed treatment allocations
#' @param X covariate matrix (each row is a single unit's covariates)
#' @param combineWithLinear augment calibration with linear covariates (defaults to TRUE; cannot be combined with linear regression for muhat)
#'
#' @return estimate of sample average treatment effect & asymptotically valid confidence interval
#'
#' @export

calibrateOB = function(muhat,
                       Y,
                       Z,
                       X,
                       combineWithLinear = TRUE)
{
  # Number of total units
  N = length(Z)
  
  # Make sure that Z is binary
  assertthat::assert_that(sort(unique(Z)) == c(0, 1))
  
  # Make sure that dimensions of Y, X, and Z are compatible
  assertthat::assert_that(length(Y) == N & dim(X)[1] == N)
  
  # Make sure that the model is a valid type of glm
  cat("TO DO\n")
  
  # Format dataframe
  data = data.frame(X)
  data$Y = Y
  data$Z = Z
  
  # Train the prediction functions muhatTreated and muhatControl
  mu1 = glm(Y ~ .-Z, family=muhat, subset(data, Z==1))
  mu0 = glm(Y ~ .-Z, family=muhat, subset(data, Z==0))
  
  if(combineWithLinear == FALSE)
  {
    # Using the regressions, form imputed treatment & control outcomes for each individual
    imputed = data.frame(TreatmentOutcomes = predict(mu1, data, "r"),
                         ControlOutcomes = predict(mu0, data, "r"),
                         Z = data$Z,
                         Y = data$Y)
  }else{
    # Add the original covariates back to the dataset
    imputed = cbind(imputed, data[, !(names(data) %in% c('Y', 'Z'))])
  }
  
  # Linearly calibrate 
  lin1 = lm(Y ~ .-Z, subset(imputed, Z==1))
  lin0 = lm(Y ~ .-Z, subset(imputed, Z==0))
  
  # Form Oaxaca-Blinder estimator using the predictions from the linear regression
  tau.hat = mean(predict(lin1, imputed) - predict(lin0, imputed))
  CI = tau.hat + t.test(residuals(lin1), residuals(lin0))$conf.int
  
  # Return the point estimat of sample average treatment effect & asymptotically valid confidence interval
  return(list(point = tau.hat, lower = CI[1], upper = CI[2]))
}