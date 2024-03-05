# Helper Function --------------------------------------------------------------

calc_ncp_chi2 <- function(alpha, power, df = 1){

  # First compute critical value based on alpha level and DF
  crit_val <- qchisq(alpha, df = df, lower.tail = FALSE)

  # Start with a small value and increase until it results in sufficient power
  ncp_start <- 1
  while(pchisq(crit_val, df = df, lower.tail = FALSE, ncp = ncp_start) < power){
    ncp_start <- ncp_start + 0.001
  }

  # Find exact NCP value to return
  ncp_out <- uniroot(function(ncp)
    return(pchisq(crit_val, df = df, lower.tail = FALSE, ncp = ncp) - power),
    c(0, ncp_start))$root

  # Return exact NCP value
  return(ncp_out)
} # End calc_ncp_chi2()

calc_ncp_F <- function(alpha_input, K_input, power_input){

  # First compute F-score based on K total, alpha level, and DF
  F_score <- qf(1 - alpha_input, df1 = 2, df2 = K_input*2 - 2*2, ncp = 0,
                lower.tail = TRUE, log.p = FALSE)

  # Start with a small value and increase until it results in sufficient power
  ncp_start <- 1
  while(1 - pf(F_score, df1 = 2, df2 = K_input*2 - 2*2, ncp_start,
               lower.tail = TRUE, log.p = FALSE) < power_input){
    ncp_start <- ncp_start + 0.001
  }

  # Find exact NCP value to return
  ncp_out <- uniroot(function(ncp)
   return(pf(F_score, df1 = 2, df2 = K_input*2 - 2*2,
             lower.tail = FALSE, ncp = ncp) - power_input),
   c(0, ncp_start))$root

  # Return exact value
  return(ncp_out)
} # End calc_ncp_F()


# Calculates covariance between betas
calCovbetas <- function(vars, rho01, rho2, cv, sigmaz.square, m, Q){
  sigmaE <- constrRiE(rho01, rho2, Q, vars)
  sigmaP <- constrRiP(rho01, Q, vars)
  tmp <- solve(diag(1,Q) - cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP) ))
  covMatrix <- 1/(m*sigmaz.square)*(sigmaE + m*sigmaP)%*%tmp
  covMatrix <- (covMatrix + t(covMatrix))/2  # symmerize the off-diagonal
  return(covMatrix)
}

# Constructs covariance matrix Sigma_E for Y_i
constrRiE <- function(rho01, rho2, Q, vars){
  rho0k <- diag(rho01)
  SigmaE_Matrix <- diag((1 - rho0k)*vars)
  for(row in 1:Q) {
    for(col in 1:Q) {
      if(row != col){
        SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col] - rho01[row,col])
      }
    }
  }
  return(SigmaE_Matrix)
}

# Constructs covariance matrix Sigma_phi for Y_i
constrRiP <- function(rho01, Q, vars) {
  rho0k <- diag(rho01)
  SigmaP_Matrix <- diag(rho0k*vars)
  for(row in 1:Q) {
    for(col in 1:Q) {
      if(row != col){
        SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
      }
    }
  }
  return(SigmaP_Matrix)
}
