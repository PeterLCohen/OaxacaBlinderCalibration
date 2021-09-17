# OaxacaBlinderCalibration
Calibration procedure for asymptotically non-inferior imputation estimators.  See "No-harm calibration for generalized Oaxaca-Blinder estimators" (Preprint available at: https://arxiv.org/abs/2012.09246v2)

Authors: Peter Cohen & Colin Fogarty

# calibrateOB.R:
A short script that allows a user to input:
  *) observed outcomes (Y)
  *) covariates (X)
  *) binary vector of treatment allocation (Z)
Returns an estimate of the treatment effect by regression-adjusting according to the user's modeling choice (muhat)
Note: The model muhat must be a valid family of glm

# OaxacaBlinderExample.R:
An example script to demonstrate using calibrateOB.R on data.

# SimulationCode.R:
Reproduces the simulation of "No-harm calibration for generalized Oaxaca-Blinder estimators"
