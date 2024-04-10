# Code Appendix for the CIRCL Study Design Example

# Setup ------------------------------------------------------------------------

# Clear environment
rm(list = ls(all.names = TRUE))

# Sets the repository to be the CodeAppendix folder; assumes current
# directory is the main GitHub repository on your computer
setwd("./CodeAppendix")

# Source file that automatically checks and downloads required packages
source("./RequiredPackages.R")

# Source code for Method 5: Conjunctive IU test from
#     https://github.com/siyunyang/coprimary_CRT
source("./powerSampleCal_varCluster_ttest.R")

# Source helper functions
source("./HelperFunctions.R")

# Input Parameters -------------------------------------------------------------
# Creating objects for the input design parameters
K_input     <- 15 # number of clusters
m_input     <- 300 # number of participants
power_input <- 0.8 # statistical power
alpha_input <- 0.05 # family-wise alpha value

# Treatment effects for outcomes
beta1_input <- 0.1 # Y1, controlled BP
beta2_input <- 0.1 # Y2, reach of Kaiser bundle

# Correlations
rho01_input <- 0.025 # ICC(Y1)
rho02_input <- 0.025 # ICC(Y2)
rho1_input  <- 0.01 # Inter-subject between-endpoint ICC
rho2_input  <- 0.05 # Intra-subject between-endpoint ICC
varY1_input <- 0.23 # total variance of Y1
varY2_input <- 0.25 # total variance of Y2
r_input     <- 1 # treatment ratio

# Find non-centrality parameter that corresponds with inputted power and alpha
# Based on Chi2 Distribution or F
ncp         <- calc_ncp_chi2(alpha_input, power_input, df = 1)
ncp2        <- calc_ncp_chi2(alpha_input, power_input, df = 2) # for 2 DF
ncpF        <- calc_ncp_F(alpha_input, K_input, power_input)

# Method 1: P-Value Adjustments ------------------------------------------------

# Adjusted p-values for three methods based on inputted alpha level
alpha_B <- alpha_input/2 # Bonferroni
alpha_S <- 1 - (1 - alpha_input)^(1/2) # Sidak
alpha_D <- 1 - (1 - alpha_input)^(1/(2^(1 - rho2_input))) # D/AP

# Critical values for three methods
cv_B <- qchisq(1 - alpha_B, df = 1, ncp = 0) # Bonferroni
cv_S <- qchisq(1 - alpha_S, df = 1, ncp = 0) # Sidak
cv_D <- qchisq(1 - alpha_D, df = 1, ncp = 0) # D/AP

# Lambda values for calculating K and m for three methods
lambda_B <- calc_ncp_chi2(alpha_B, power_input, df = 1)
lambda_S <- calc_ncp_chi2(alpha_S, power_input, df = 1)
lambda_D <- calc_ncp_chi2(alpha_D, power_input, df = 1)

# Power for Method 1
method1_lambda1 <- (beta1_input^2)/(2*(varY1_input/(K_input*m_input))*
                                      (1 + (m_input-1)*rho01_input))
method1_lambda2 <- (beta2_input^2)/(2*(varY2_input/(K_input*m_input))*
                                      (1 + (m_input-1)*rho02_input))

method1_power_bonf  <- c(1 - pchisq(cv_B, 1, ncp = method1_lambda1,
                                    lower.tail = TRUE), # Power for Y1
                         1 - pchisq(cv_B, 1, ncp = method1_lambda2,
                                    lower.tail = TRUE)) # Power for Y2
method1_power_sidak <- c(1 - pchisq(cv_S, 1, ncp = method1_lambda1,
                                    lower.tail = TRUE), # Power for Y1
                         1 - pchisq(cv_S, 1, ncp = method1_lambda2,
                                    lower.tail = TRUE)) # Power for Y2
method1_power_dap   <- c(1 - pchisq(cv_D, 1, ncp = method1_lambda1,
                                    lower.tail = TRUE), # Power for Y1
                         1 - pchisq(cv_D, 1, ncp = method1_lambda2,
                                    lower.tail = TRUE)) # Power for Y2

# K for Method 1
method1_K_bonf  <- c((2*lambda_B*varY1_input*(1 + (m_input - 1)*rho01_input))/
                       (m_input*(beta1_input^2)),
                     (2*lambda_B*varY2_input*(1 + (m_input - 1)*rho02_input))/
                       (m_input*(beta2_input^2)))
method1_K_sidak <- c((2*lambda_S*varY1_input*(1 + (m_input - 1)*rho01_input))/
                       (m_input*(beta1_input^2)),
                     (2*lambda_S*varY2_input*(1 + (m_input - 1)*rho02_input))/
                       (m_input*(beta2_input^2)))
method1_K_dap   <- c((2*lambda_D*varY1_input*(1 + (m_input - 1)*rho01_input))/
                       (m_input*(beta1_input^2)),
                     (2*lambda_D*varY2_input*(1 + (m_input - 1)*rho02_input))/
                       (m_input*(beta2_input^2)))

# m for Method 1
method1_m_bonf  <- c((2*lambda_B*varY1_input*(1 - rho01_input))/
                       ((beta1_input^2)*K_input - 2*lambda_B*varY1_input*rho01_input),
                     (2*lambda_B*varY2_input*(1 - rho02_input))/
                       ((beta2_input^2)*K_input - 2*lambda_B*varY2_input*rho02_input))
method1_m_sidak <- c((2*lambda_S*varY1_input*(1 - rho01_input))/
                       ((beta1_input^2)*K_input - 2*lambda_S*varY1_input*rho01_input),
                     (2*lambda_S*varY2_input*(1 - rho02_input))/
                       ((beta2_input^2)*K_input - 2*lambda_S*varY2_input*rho02_input))
method1_m_dap   <- c((2*lambda_D*varY1_input*(1 - rho01_input))/
                       ((beta1_input^2)*K_input - 2*lambda_D*varY1_input*rho01_input),
                     (2*lambda_D*varY2_input*(1 - rho02_input))/
                       ((beta2_input^2)*K_input - 2*lambda_D*varY2_input*rho02_input))

# Method 2: Combined Outcomes  -------------------------------------------------

# Combined outcome effect size
betac_input <- beta1_input + beta2_input

# Calculate variance for the combined outcome
varYc_input <- round(varY1_input + varY2_input + 2*rho2_input*
                       sqrt(varY1_input)*sqrt(varY2_input), 2)

# Endpoint specific ICC for combined outcome
rhoC_input  <- (rho01_input*varY1_input + rho02_input*varY2_input +
                  2*rho1_input*sqrt(varY1_input*varY2_input))/
  (varY1_input + varY2_input + 2*rho2_input*sqrt(varY1_input*varY2_input))

# Critical value for Method 2
method2_crit    <- qchisq(p = alpha_input, df = 1, lower.tail = FALSE)

# Power for Method 2
method2_lambda  <- (betac_input^2)/(2*(varYc_input/(K_input*m_input))*
                                     (1 + (m_input - 1)*rhoC_input))
method2_power   <- 1 - pchisq(method2_crit, 1,
                            ncp = method2_lambda, lower.tail = TRUE)

# K for Method 2
method2_K       <- (2*ncp*varYc_input*(1 + (m_input - 1)*
                                  rhoC_input))/(m_input*(betac_input^2))

# m for Method 2
method2_m       <- (2*ncp*varYc_input*
                      (1 - rhoC_input))/((betac_input^2)*
                                     K_input - (2*ncp*varYc_input*rhoC_input))

# Method 3: Single 1-DF Combined Test ------------------------------------------

# Critical value for Method 3
method3_crit <- qchisq(p = alpha_input, df = 1, lower.tail = FALSE)

# Power for Method 3

# Calculate test statistic for first and second outcome
Z1.sq <- (beta1_input^2)/(((2*varY1_input)/(K_input*m_input))*
                            (1 + (m_input - 1)*rho01_input)) # Z1^2
Z2.sq <- (beta2_input^2)/(((2*varY2_input)/(K_input*m_input))*
                            (1 + (m_input - 1)*rho02_input)) # Z2^2

CorrZ1Z2 <- (rho2_input + (m_input - 1)*rho1_input)/
  sqrt((1 + (m_input - 1)*rho01_input)*(1 + (m_input - 1)*rho02_input))

method3_lambda <- ((sqrt(Z1.sq) + sqrt(Z2.sq))^2)/(2*(1 + CorrZ1Z2))

method3_power   <- 1 - pchisq(method3_crit, ncp = method3_lambda,
                              df = 1, lower.tail = TRUE)

# K for Method 3
method3_K <- (2*ncp*(1 + CorrZ1Z2))/
  (sqrt((beta1_input^2)/((2*varY1_input/m_input)*(1+(m_input-1)*rho01_input))) +
     sqrt((beta2_input^2)/((2*varY2_input/m_input)*(1+(m_input-1)*rho02_input))))^2

# m for Method 3
power.eq.m <- function(m, para){ # Function for equation solve for m
  # Calculate test statistic for first and second outcome
  Z1.sq <- (para[1]^2)/((2*para[4]/(m*para[3]))*(1 + (m - 1)*para[6])) # Z1^2
  Z2.sq <- (para[2]^2)/((2*para[5]/(m*para[3]))*(1 + (m - 1)*para[7])) # Z1^2
  # Calculate Corr(Z1, Z2)
  CorrZ1Z2 <- (para[9] + (m - 1)*para[8])/
    sqrt((1+(m-1)*para[6])*(1+(m- 1)*para[7]))
  # Non-centrality parameter
  lambda <- ((sqrt(Z1.sq) + sqrt(Z2.sq))^2) / (2*(1 + CorrZ1Z2)) - para[10]
}

# Find m using multiroot function
find.para.m <- multiroot(power.eq.m, start = 10,
                         para = c(beta1_input, beta2_input, K_input,
                                  varY1_input, varY2_input,
                                  rho01_input, rho02_input,
                                  rho1_input, rho2_input,
                                  ncp),
                         positive = TRUE)
method3_m <- find.para.m$root

# Method 4: Disjunctive 2-DF Test (Chi2 Distribution) --------------------------

# Critical value for method 4, df = 2
method4_crit_chi2 <- qchisq(1 - alpha_input, df = 2,
                            ncp = 0, lower.tail = TRUE, log.p = FALSE)

# Calculate VIFs
VIF1 <- 1 + (m_input - 1)*rho01_input
VIF2 <- 1 + (m_input - 1)*rho02_input
VIF12 <- rho2_input + (m_input - 1)*rho1_input

# Power for Method 4 Chi2
method4_lambda <- K_input*m_input*((beta1_input^2)*varY2_input*VIF2 +
                                     (beta2_input^2)*varY1_input*VIF1 -
                 2*beta1_input*beta2_input*sqrt(varY1_input)*
                   sqrt(varY2_input)*VIF12)/(2*varY1_input*
                                               varY2_input*(VIF1*VIF2-VIF12^2))

method4_power_chi2 <- 1 - pchisq(method4_crit_chi2, df = 2,
                                 ncp = method4_lambda, lower.tail = TRUE)

# K for Method 4 Chi2
method4_K_chi2 <- (ncp2*2*varY1_input*varY2_input*(VIF1*VIF2-VIF12^2))/
  (m_input*((beta1_input^2)*varY2_input*VIF2 +
                             (beta2_input^2)*varY1_input*VIF1 -
                             2*beta1_input*beta2_input*sqrt(varY1_input)*
                             sqrt(varY2_input)*VIF12))

# m for Method 4 Chi2
# Temporary equation to use multiroot() on to solve for m
find.para.m <- function(m, para){
  VIF1 <- 1 + (m - 1)*para[6]
  VIF2 <- 1 + (m - 1)*para[7]
  VIF12 <- para[9] + (m - 1)*para[8]
  lambda <- para[3]*m*(para[1]^2*para[5]*VIF2+para[2]^2*para[4]*VIF1 -
                         2*para[1]*para[2]*sqrt(para[4]*para[5])*VIF12)/
    (2*para[4]*para[5]*(VIF1*VIF2 - VIF12^2)) - para[10]
}

# Finding value
find.para.m_chi2 <- multiroot(find.para.m, start = 3,
                              para = c(beta1_input, beta2_input, K_input,
                                       varY1_input, varY2_input,
                                       rho01_input, rho02_input,
                                       rho1_input, rho2_input, ncp2),
                         positive = TRUE)
method4_m_chi2 <- find.para.m_chi2$root # Cluster size required

# Method 4: Disjunctive 2-DF Test (F Distribution) -----------------------------

# Value needed for power calculations using F distribution
Fscore <- qf(1 - alpha_input, df1 = 2, df2 = K_input*2 - 2*2, ncp = 0,
             lower.tail = TRUE, log.p = FALSE)

# Calculate VIFs
VIF1 <- 1 + (m_input - 1)*rho01_input
VIF2 <- 1 + (m_input - 1)*rho02_input
VIF12 <- rho2_input + (m_input - 1)*rho1_input

# Power for Method 4
method4_lambda <- K_input*m_input*((beta1_input^2)*varY2_input*VIF2 +
                                     (beta2_input^2)*varY1_input*VIF1 -
                                     2*beta1_input*beta2_input*sqrt(varY1_input)*
                                     sqrt(varY2_input)*VIF12)/(2*varY1_input*
                                                                 varY2_input*(VIF1*VIF2-VIF12^2))

method4_power_F <- 1 - pf(Fscore, df1 = 2, df2 = K_input*2 - 2*2,
                          method4_lambda, lower.tail = TRUE, log.p = FALSE)

# K for Method 4
# Defining necessary parameters based on input values
betas <- c(beta1_input, beta2_input) # Beta vector
r_alt <- 1/(1 + 1)
Q <- 2 # Number of outcomes, could extend this to more than 2 in the future
vars <- c(varY1_input, varY2_input) # Vector of variances
rho01_mat <- matrix(c(rho01_input, rho1_input,
                      rho1_input, rho02_input), 2, 2)
rho2_mat <- matrix(c(1, rho2_input, rho2_input, 1), 2, 2)
clus_var <- 0 # cluster variation, placeholder for future extensions
pred.power <- 0 # Initializing predictive power
sigmaz.square <- r_alt*(1-r_alt) # Variance of treatment assignment
K_total_temp <- 2*Q
while(pred.power < power_input){ # Iterate while the current predictive power is lower than desired power
  K_total_temp <- K_total_temp + 1
  omega_temp <- calCovbetas(vars, rho01_mat, rho2_mat, clus_var,
                            sigmaz.square, m_input, Q)
  tau_temp <- K_total_temp*t(betas) %*% solve(omega_temp) %*% betas
  Fscore_temp <- qf(1 - alpha_input, df1 = Q, df2 = K_total_temp - 2*Q, ncp = 0,
               lower.tail = TRUE, log.p = FALSE)
  pred.power <- round(1 - pf(Fscore_temp, df1 = Q, df2 = K_total_temp - 2*Q, tau_temp,
                             lower.tail = TRUE, log.p = FALSE), 4)
}
method4_K_F <- K_total_temp/2

# m for Method 4
# Temporary equation to use multiroot() on to solve for m
find.para.m <- function(m, para){
  VIF1 <- 1 + (m - 1)*para[6]
  VIF2 <- 1 + (m - 1)*para[7]
  VIF12 <- para[9] + (m - 1)*para[8]
  lambda <- para[3]*m*(para[1]^2*para[5]*VIF2+para[2]^2*para[4]*VIF1 -
                         2*para[1]*para[2]*sqrt(para[4]*para[5])*VIF12)/
    (2*para[4]*para[5]*(VIF1*VIF2 - VIF12^2)) - para[10]
}

# Finding value
find.para.m_F <- multiroot(find.para.m, start = 3,
                           para = c(beta1_input, beta2_input, K_input,
                                    varY1_input, varY2_input,
                                    rho01_input, rho02_input,
                                    rho1_input, rho2_input, ncpF),
                           positive = TRUE)
method4_m_F <- find.para.m_F$root # Cluster size required

# Method 5: Conjunctive IU Test (T-Distribution) -------------------------------

# Power for Method 5
method5_power <- calPower_ttestIU(betas = c(beta1_input, beta2_input),
                                  K = 2, # In this function, K = # of outcomes
                                  m = m_input,
                                  deltas = c(0, 0),
                                  vars = c(varY1_input, varY2_input),
                                  rho01 = matrix(c(rho01_input, rho1_input,
                                                   rho1_input, rho02_input),
                                                 2, 2),
                                  rho2 = matrix(c(1, rho2_input,
                                                  rho2_input, 1),
                                                2, 2),
                                  r = 0.5,
                                  N = 2*K_input,
                                  alpha = alpha_input)

# K for Method 5
method5_K <- calSampleSize_ttestIU(betas = c(beta1_input, beta2_input),
                                   m = m_input,
                                   power = power_input,
                                   deltas = c(0, 0),
                                   vars = c(varY1_input, varY2_input),
                                   rho01 = matrix(c(rho01_input, rho1_input,
                                                    rho1_input, rho02_input),
                                                  2, 2),
                                   rho2 = matrix(c(1, rho2_input,
                                                   rho2_input, 1),
                                                 2, 2),
                                   r = 0.5,
                                   K = 2,
                                   alpha = alpha_input)/2 # Divide by 2 because
                                                          # function returns
                                                          # total number of clusters,
                                                          # which is 2K

# m for Method 5
# Write custom function for finding m, uses same logic as K
calculateClusterSize_IUTest <- function(betas,
                                        deltas,
                                        vars,
                                        rho01,
                                        rho2,
                                        K,
                                        r,
                                        N,
                                        alpha,
                                        power,
                                        dist = "T"){

  # Initialize m and predictive power
  m <- 0
  pred.power <- 0
  while(pred.power < power){
    m <- m + 1
    pred.power <- calPower_ttestIU(betas = betas, # Treatment effect for outcomes
                                   deltas = deltas, # Non-inferiority margins
                                   vars = vars, # Marginal variances for endpoints
                                   rho01 = rho01, # KxK matrix for correlations
                                   rho2 = rho2, # KxK matrix for correlations
                                   N = N, # Number of clusters
                                   r = r, # Proportion of clusters in intervention arm
                                   m = m, # Mean cluster size
                                   K = K, # Number of endpoints
                                   alpha = alpha, # Type I error rate upper bound
                                   dist = dist)
    if(m > 10000){
      m <- Inf
      break
    }
    #cat("m = ", m, "\n", "power = ", pred.power, "\n")
  }
  return(m)
} # End calculateCluster_IUTest()

method5_m <- calculateClusterSize_IUTest(betas = c(beta1_input, beta2_input),
                                         N = 2*K_input,
                                         power = power_input,
                                         deltas = c(0, 0),
                                         vars = c(varY1_input, varY2_input),
                                         rho01 = matrix(c(rho01_input, rho1_input,
                                                          rho1_input, rho02_input),
                                                        2, 2),
                                         rho2 = matrix(c(1, rho2_input,
                                                         rho2_input, 1),
                                                       2, 2),
                                         r = 0.5,
                                         K = 2,
                                         alpha = alpha_input)

# Method 5: Conjunctive IU Test (MVN-Distribution) -----------------------------

# Power for Method 5
method5_power_MVN <- calPower_ttestIU(betas = c(beta1_input, beta2_input),
                                      K = 2, # In this function, K = # of outcomes
                                      m = m_input,
                                      deltas = c(0, 0),
                                      vars = c(varY1_input, varY2_input),
                                      rho01 = matrix(c(rho01_input, rho1_input,
                                                       rho1_input, rho02_input),
                                                     2, 2),
                                      rho2 = matrix(c(1, rho2_input,
                                                      rho2_input, 1),
                                                    2, 2),
                                      r = 0.5,
                                      N = 2*K_input,
                                      alpha = alpha_input,
                                      dist = "MVN")

# K for Method 5
method5_K_MVN <- calSampleSize_ttestIU(betas = c(beta1_input, beta2_input),
                                       m = m_input,
                                       power = power_input,
                                       deltas = c(0, 0),
                                       vars = c(varY1_input, varY2_input),
                                       rho01 = matrix(c(rho01_input, rho1_input,
                                                        rho1_input, rho02_input),
                                                      2, 2),
                                       rho2 = matrix(c(1, rho2_input,
                                                       rho2_input, 1),
                                                     2, 2),
                                       r = 0.5,
                                       K = 2,
                                       alpha = alpha_input,
                                       dist = "MVN")/2 # Divide by 2 to get
                                                       # clusters per treatment group

method5_m_MVN <- calculateClusterSize_IUTest(betas = c(beta1_input, beta2_input),
                                             N = 2*K_input,
                                             power = power_input,
                                             deltas = c(0, 0),
                                             vars = c(varY1_input, varY2_input),
                                             rho01 = matrix(c(rho01_input, rho1_input,
                                                              rho1_input, rho02_input),
                                                            2, 2),
                                             rho2 = matrix(c(1, rho2_input,
                                                             rho2_input, 1),
                                                           2, 2),
                                             r = 0.5,
                                             K = 2,
                                             alpha = alpha_input,
                                             dist = "MVN")

# Summary of Results -----------------------------------------------------------

# Table for final power, K, and m calculations
summaryTable <- tibble(`Method` = c("1. P-Value Adjustments",
                                    "a. Bonferroni",
                                    "b. Sidak",
                                    "c. D/AP",
                                    "2. Combined Outcomes",
                                    "3. Single 1-df Combined Test",
                                    "4. Disjunctive 2-df Test",
                                    "a. Chi-2 Distribution",
                                    "b. F Distribution",
                                    "5. Conjunctive IU Test",
                                    "a. Multivariate Normal Distribution",
                                    "b. T Distribution"),
                       `Power`  = c(NA,
                                    min(method1_power_bonf),
                                    min(method1_power_sidak),
                                    min(method1_power_dap),
                                    method2_power,
                                    method3_power,
                                    NA,
                                    method4_power_chi2,
                                    method4_power_F,
                                    NA,
                                    method5_power_MVN,
                                    method5_power),
                       `K`      = c(NA,
                                    max(method1_K_bonf),
                                    max(method1_K_sidak),
                                    max(method1_K_dap),
                                    method2_K,
                                    method3_K,
                                    NA,
                                    method4_K_chi2,
                                    method4_K_F,
                                    NA,
                                    method5_K_MVN,
                                    method5_K),
                       `m`      = c(NA,
                                    max(method1_m_bonf),
                                    max(method1_m_sidak),
                                    max(method1_m_dap),
                                    method2_m,
                                    method3_m,
                                    NA,
                                    method4_m_chi2,
                                    method4_m_F,
                                    NA,
                                    method5_m_MVN,
                                    method5_m))

# Table displaying p-value adjustment method comparisons between Y1 and Y2
pValueTable <- tibble(`Method`    = c("Bonferroni",
                                      "Sidak",
                                      "D/AP"),
                      `Power(Y1)` = c(method1_power_bonf[1],
                                      method1_power_sidak[1],
                                      method1_power_dap[1]),
                      `Power(Y2)` = c(method1_power_bonf[2],
                                      method1_power_sidak[2],
                                      method1_power_dap[2]),
                      `K(Y1)`     = c(method1_K_bonf[1],
                                      method1_K_sidak[1],
                                      method1_K_dap[1]),
                      `K(Y2)`     = c(method1_K_bonf[2],
                                      method1_K_sidak[2],
                                      method1_K_dap[2]),
                      `m(Y1)`     = c(method1_m_bonf[1],
                                      method1_m_sidak[1],
                                      method1_m_dap[1]),
                      `m(Y2)`     = c(method1_m_bonf[2],
                                      method1_m_sidak[2],
                                      method1_m_dap[2])
)

View(summaryTable)
View(summaryTable %>% mutate(Power = round(100*Power, 2),
                             K = ceiling(K),
                             m = ceiling(m)))
