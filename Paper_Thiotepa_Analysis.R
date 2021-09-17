# Read in the data
dat = read.csv("Bladder_Cancer_2.csv")
colnames(dat) <- c('Treatment', 'Number_Recurrences','Months_Follow_Up', 'Tumors_Baseline', 'Largest_Tumor_Baseline')
dat = dat[dat$Treatment!="Pyridoxine",]

# Treatment allocation
Treatment = dat$Treatment

# Number of recurrences per patient
Number_Recurrences = dat$Number_Recurrences

# The number of patients in the control group
nPlacebo = table(dat$Treatment)[1]
cat("# control", nPlacebo, '\n')

# The number of patients in the thiotepa treatment group
nThiotepa = table(dat$Treatment)[2]
cat("# Thiotepa", nThiotepa, '\n')

cat('Study size', nPlacebo + nThiotepa, '\n')


# Observed outcomes
yPlacebo = Number_Recurrences[Treatment == "Placebo"] # control population observed outcomes
yThiotepa = Number_Recurrences[Treatment=="Thiotepa"] # thiotepa treatment population ovserved outcomes

################################################################################
###                   The unadjusted difference in means                     ###
################################################################################
tau_unadj =  mean(yThiotepa) - mean(yPlacebo) 



################################################################################
### The uncalibrated generalized Oaxaca-Blinder estimator of Guo and Basse   ###
################################################################################
# Base prediction functions of outcome given log(months), # initial tumors, size initial largest tumor
muPlacebo = glm(Number_Recurrences~log(Months_Follow_Up) + Tumors_Baseline+Largest_Tumor_Baseline, family = poisson, data = subset(dat, Treatment == "Placebo"))
muThiotepa = glm(Number_Recurrences~log(Months_Follow_Up) + Tumors_Baseline+Largest_Tumor_Baseline, family = poisson, data = subset(dat, Treatment == "Thiotepa"))

# Fitted values of muPlacebo and muThiotepa
fPlacebo = predict(muPlacebo, dat, "r")
fThiotepa = predict(muThiotepa, dat, "r")

tau_gob = mean(fThiotepa - fPlacebo)

# Residuals of muPlacebo and muThiotepa
rPlacebo = (Number_Recurrences - fPlacebo)[Treatment=="Placebo"]
# rThiotepa = (Number_Recurrences - fThiotepa)[Treatment=="Thiotepa"]
rThiotepa = (Number_Recurrences - fThiotepa)[Treatment=="Thiotepa"]


################################################################################
###               Guo & Basses's singly calibrated estimator                 ###
################################################################################
# Form pseudo-covariates using the fitted values of muPlacebo and muThiotepa
dattilde = data.frame(Treatment, Number_Recurrences, fPlacebo, fThiotepa)

# Joint linear calibration of prediction functions
mu1Placebo = lm(Number_Recurrences~fPlacebo, subset(dattilde, Treatment == "Placebo"))
mu1Thiotepa = lm(Number_Recurrences~fThiotepa, subset(dattilde, Treatment == "Thiotepa"))
f1Placebo = predict(mu1Placebo, dattilde)
f1Thiotepa = predict(mu1Thiotepa, dattilde)
tau_cal1 = mean(f1Thiotepa - f1Placebo)

# Residuals after joint linear calibration of the prediction functions
r1Placebo = (Number_Recurrences - f1Placebo)[Treatment=="Placebo"]
r1Thiotepa = (Number_Recurrences - f1Thiotepa)[Treatment=="Thiotepa"]


################################################################################
###           Cohen & Fogarty's joint linearly calibrated estimator          ###
################################################################################
# Form pseudo-covariates using the fitted values of muPlacebo and muThiotepa
dattilde = data.frame(Treatment, Number_Recurrences, fPlacebo, fThiotepa)

# Joint linear calibration of prediction functions 
mu2Placebo = lm(Number_Recurrences~fPlacebo + fThiotepa, subset(dattilde, Treatment == "Placebo"))
mu2Thiotepa = lm(Number_Recurrences~fPlacebo + fThiotepa, subset(dattilde, Treatment == "Thiotepa"))
f2Placebo = predict(mu2Placebo, dattilde)
f2Thiotepa = predict(mu2Thiotepa, dattilde)
tau_cal2 = mean(f2Thiotepa - f2Placebo)

# Residuals after joint linear calibration of the prediction functions
r2Placebo = (Number_Recurrences - f2Placebo)[Treatment=="Placebo"]
r2Thiotepa = (Number_Recurrences - f2Thiotepa)[Treatment=="Thiotepa"]

################################################################################
###                         Variance Estimators                              ###
################################################################################
var_unajd = var(yPlacebo)/nPlacebo + var(yThiotepa)/nThiotepa # Neyman variance estimator of the unadjusted difference in means
var_gob = var(rPlacebo)/nPlacebo + var(rThiotepa)/nThiotepa # Variance estimator of the uncalibrated generalized Oaxaca-Blinder estimator of Guo and Basse
var_cal1 = var(r1Placebo)/nPlacebo + var(r1Thiotepa)/nThiotepa # Variance estimator of the Guo & Basse's singly calibrated estimator
var_cal2 = var(r2Placebo)/nPlacebo + var(r2Thiotepa)/nThiotepa # Variance estimator of the linearly calibrated estimator

################################################################################
###                               Print Table                                ###
################################################################################

cat("$\\tauhat_{unadj}$ \t&  ", round(tau_unadj, digits = 3),"                &  ", round(var_unajd, digits = 3) ,"    \\\\
$\\tauhat_{gOB}$     \t&  ", round(tau_gob, digits = 3) ,"                &  ", round(var_gob, digits = 3) ,"    \\\\
$\\tauhat_{GBcal}$   \t&  ", round(tau_cal1, digits = 3),"                &  ", round(var_cal1, digits = 3) ,"    \\\\
$\\tauhat_{cal}$     \t&  ", round(tau_cal2, digits = 3),"                &  ", round(var_cal2, digits = 3)    )