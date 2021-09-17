# OaxacaBlinderCalibration
Calibration procedure for asymptotically non-inferior imputation estimators.  See "No-harm calibration for generalized Oaxaca-Blinder estimators" (Preprint available at: https://arxiv.org/abs/2012.09246v2)

Authors: Peter Cohen & Colin Fogarty

# calibrateOB.R:
A short script that allows a user to input:
  * observed outcomes (Y)
  * covariates (X)
  * binary vector of treatment allocation (Z)
Returns an estimate of the treatment effect by regression-adjusting according to the user's modeling choice (muhat)
Note: The model muhat must be a valid family of glm

# OaxacaBlinderExample.R:
An example script to demonstrate using calibrateOB.R on data.

# PoissonSimulations
A subdirectory containing the simulations based upon Poisson regression in the ArXiv v2 preprint.

# Logistic_twoPhaseSampling_SATE.R
Simulation code to recreate the logistic regression simulations of the ArXiv v2 preprint.

# Paper_Thiotepa_Analysis.R
Reproducibility code to generate the analysis of thiotepa for prevention of reccurent bladder tumors.  relies upon *Bladder_Cancer_2.csv* which is data extracted from *Data* by Andrews and Herzberg (Springer 1985).
