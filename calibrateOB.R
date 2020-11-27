#' Calibration of Oaxaca blinder estimators 
#'
#' Implements calibration method of Cohen & Fogarty (2020) to build non-inferior
#'  Oaxaca-Blinder imputation estimators.
#'
#' @param muhatTreated imputation algorithm to use in the treated group
#' @param muhatControl imputation algorithm to use in the treated group
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
  
}