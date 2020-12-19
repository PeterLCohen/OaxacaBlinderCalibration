# Example to demonstrate using calibration with non-linear predictions for 
#   non-inferior imputation estimators 

# install.packages('AER')
library(AER) # For the Fatalities data
data("Fatalities")

N = dim(Fatalities)[1] # Total number of 'units'
p = .5 # Proportion of treated units
n1 = floor(N*p) # Number of treated units
n0 = N - n1 # Number of control units
treatmentVec = c(rep(TRUE, n1), rep(FALSE, n0))

X = Fatalities[,c("pop", "miles", "income")] # Control for the same covariates as Guo & Basse (See Fig. 2 of their paper)
logX = log(X)
Y = Fatalities$fatal # Number of vehicle fatalities.

# Compute the non-inferior Oaxaca-blinder estimator
OB = calibrateOB(muhat = 'poisson',
                 Y = Y,
                 Z = sample(treatmentVec),
                 X = logX,
                 combineWithLinear = FALSE)
