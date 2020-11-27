#' Calibration of Oaxaca blinder estimators
#'
#' Implements calibration method of Cohen & Fogarty (2020) to build non-inferior
#'  Oaxaca-Blinder imputation estimators.
#'
#' @param muhat model to use in the treated and control groups (must be a family of GLMs)
#' @param Yobs the observed outcomes
#' @param Z the observed treatment allocations
#' @param X covariate matrix (each row is a single unit's covariates)
#'
#' @return estimate of sample average treatment effect
#'
#' @export

calibrateOB = function(muhatTreated,
                       muhatControl,
                       Yobs,
                       Z,
                       X)
{
  # Number of total units
  N = length(Z)
  
  # Make sure that the treated and control models inherit from the lm class
  # assertthat::assert_that('lm' %in% class(muhatTreated))
  # assertthat::assert_that('lm' %in% class(muhatControl))
  
  # Make sure that Z is binary
  assertthat::assert_that(sort(unique(Z)) == c(0, 1))
  
  # Make sure that dimensions of Yobs, X, and Z are compatible
  assertthat::assert_that(length(Yobs) == N & dim(X)[1] == N)
  
  # Train the prediction functions muhatTreated and muhatControl
  mu1 = glm(Yobs ~ .-Z, family=muhatTreated, subset(data, Z==1))
  
}