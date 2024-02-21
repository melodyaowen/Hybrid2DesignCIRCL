rm(list = ls(all.names = TRUE))
setwd("/Users/melodyowen/Desktop/GitHub/PowerHybrid2/CodeAppendix")

# Code Appendix for the CIRCL Study Design Example

# Setup ------------------------------------------------------------------------

# Required packages

# Package names
packages <- c("tidyverse", "haven", "MASS",
              "geepack", "lme4", "rootSolve",
              "devtools", "pracma", "mvtnorm", "nlme",
              "knitr", "kableExtra")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Input Parameters -------------------------------------------------------------

# Creating objects for the input design parameters we want to vary
simCases <- data.frame(beta1_input = c(0.1, 0.5, 1, 5),
                       beta2_input = c(0.1, 0.5, 1, 5),
                       m_input     = c(5, 15, 50, 100),
                       rho01_input = c(0.1, 0.2, 0.4, 0.5),
                       rho02_input = c(0.1, 0.2, 0.4, 0.5),
                       varY1_input = c(2, 5, 8, 10),
                       varY2_input = c(2, 5, 8, 10),
                       rho1_input  = c(0.01, 0.03, 0.05, 0.1),
                       rho2_input  = c(0.05, 0.1, 0.2, 0.3))

#K_input <- 15
#alpha_input <- 0.05 # family-wise alpha value

allCases <- expand.grid(simCases) %>%
  mutate(VIF1 = 1 + (m_input-1)*rho01_input,
         VIF2 = 1 + (m_input-1)*rho02_input,
         VIF12 = rho2_input + (m_input-1)*rho1_input,
         beta_ratio = beta1_input/beta2_input,
         var_ratio = varY1_input/varY2_input
  )

# Method 1: Combined Outcomes  -------------------------------------------------

method1_fun <- function(beta1_input, beta2_input,
                        K_input, m_input,
                        rho01_input, rho02_input,
                        varY1_input, varY2_input, 
                        rho1_input, rho2_input,
                        alpha_input){
  # Combined outcome effect size
  betac_input <- beta1_input + beta2_input
  
  # Calculate variance for the combined outcome
  varYc_input <- round(varY1_input + varY2_input + 2*rho2_input*
                         sqrt(varY1_input)*sqrt(varY2_input), 2)
  
  # Endpoint specific ICC for combined outcome
  rhoC_input  <- (rho01_input*varY1_input + rho02_input*varY2_input + 
                    2*rho1_input*sqrt(varY1_input*varY2_input))/
    (varY1_input + varY2_input + 2*rho2_input*sqrt(varY1_input*varY2_input))
  
  # Critical value for Method 1
  method1_crit    <- qchisq(p = alpha_input, df = 1, lower.tail = FALSE)
  
  # Power for Method 1
  method1_lambda  <- (betac_input^2)/(2*(varYc_input/(K_input*m_input))*
                                        (1 + (m_input - 1)*rhoC_input))
  method1_power   <- 1 - pchisq(method1_crit, 1, 
                                ncp = method1_lambda, lower.tail = TRUE)
  
  # method1_power_dat <- data.frame(method1_power = method1_power,
  #                                 beta1_input = beta1_input,
  #                                 beta2_input = beta2_input,
  #                                 m_input     = m_input,
  #                                 rho01_input = rho01_input,
  #                                 rho02_input = rho02_input,
  #                                 varY1_input = varY1_input,
  #                                 varY2_input = varY2_input,
  #                                 rho1_input  = rho1_input,
  #                                 rho2_input  = rho2_input
  # )
  
  # method1_power_dat <- c(method1_power,
  #                        beta1_input, beta2_input, m_input,
  #                        rho01_input, rho02_input,
  #                        varY1_input, varY2_input,
  #                        rho1_input, rho2_input
  # )
  
  return(method1_power)
}

# Method 3: Single 1-DF Combined Test ------------------------------------------

method3_fun <- function(beta1_input, beta2_input,
                        K_input, m_input,
                        rho01_input, rho02_input,
                        varY1_input, varY2_input, 
                        rho1_input, rho2_input,
                        alpha_input){
  
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
  
  # method3_power_dat <- data.frame(method3_power = method3_power,
  #                                 beta1_input = beta1_input,
  #                                 beta2_input = beta2_input,
  #                                 m_input     = m_input,
  #                                 rho01_input = rho01_input,
  #                                 rho02_input = rho02_input,
  #                                 varY1_input = varY1_input,
  #                                 varY2_input = varY2_input,
  #                                 rho1_input  = rho1_input,
  #                                 rho2_input  = rho2_input
  # )
  
  return(method3_power)
}

# Run Sim Cases ----------------------------------------------------------------

 method1_sim_power <- mapply(method1_fun,
                             beta1_input = allCases$beta1_input,
                             beta2_input = allCases$beta2_input,
                             m_input     = allCases$m_input,
                             rho01_input = allCases$rho01_input,
                             rho02_input = allCases$rho02_input,
                             varY1_input = allCases$varY1_input,
                             varY2_input = allCases$varY2_input,
                             rho1_input  = allCases$rho1_input,
                             rho2_input  = allCases$rho2_input,
                             K_input     = rep(10, nrow(simCases)),
                             alpha_input = rep(0.05, nrow(simCases))
                             )

 method3_sim_power <- mapply(method3_fun,
                             beta1_input = allCases$beta1_input,
                             beta2_input = allCases$beta2_input,
                             m_input     = allCases$m_input,
                             rho01_input = allCases$rho01_input,
                             rho02_input = allCases$rho02_input,
                             varY1_input = allCases$varY1_input,
                             varY2_input = allCases$varY2_input,
                             rho1_input  = allCases$rho1_input,
                             rho2_input  = allCases$rho2_input,
                             K_input     = rep(10, nrow(simCases)),
                             alpha_input = rep(0.05, nrow(simCases))
 )


simResults <- cbind(allCases, method1_sim_power, method3_sim_power) %>%
  rename(method2_power = method1_sim_power, method3_power = method3_sim_power) %>%
  mutate(method2_power = round(method2_power*100, 2),
         method3_power = round(method3_power*100, 2)) %>%
  mutate(compare = ifelse(method2_power == method3_power, "Equal",
                          ifelse(method2_power > method3_power, "1 Better",
                                 ifelse(method3_power > method2_power, "3 Better",
                                        "CHECK"))),
         method3minus2 = method3_power - method2_power) %>%
  relocate(compare, method3minus2, method2_power, method3_power)

nrow(filter(simResults, method2_power == 100))
nrow(filter(simResults, method3_power == 100))
nrow(simResults)

hist(simResults$method2_power)

nrow(filter(simResults, method2_power > 75, method2_power < 95))

nrow(filter(simResults, method2_power == 0))
nrow(filter(simResults, method3_power == 0))






summaryTable <- allResults %>%
  mutate(score = beta_ratio_1 + var_ratio_1 + VIF1_equal_VIF2,
         all3 = ifelse(score == 3, 1, 0),
         no3 = ifelse(score == 0, 1, 0),
         
         beta_ratio_1_only = ifelse(score == 1 & beta_ratio_1 == 1, 1, 0),
         var_ratio_1_only = ifelse(score == 1 & var_ratio_1 == 1, 1, 0),
         VIF1_equal_VIF2_only = ifelse(score == 1 & VIF1_equal_VIF2 == 1, 1, 0)) %>%
  group_by(compare) %>%
  summarize(n = n(),
            beta1_equal_beta2 = sum(beta_ratio_1),
            var1_equal_var2 = sum(var_ratio_1),
            VIF1_equal_VIF2 = sum(VIF1_equal_VIF2),
            all3 = sum(all3),
            no3 = sum(no3)) %>%
  add_row(compare = "Total", summarise(., across(where(is.numeric), sum)))

# Overall proportions of comparisons for simulation
propTable <- tibble(`Equal` = paste0(round(nrow(subsetEqual)/nrow(simResults), 2)*100, "%"),
                    `1 Better` = paste0(round(nrow(subsetMethod1Better)/nrow(simResults), 2)*100, "%"),
                    `3 Better` = paste0(round(nrow(subsetMethod3Better)/nrow(simResults), 2)*100, "%"))
View(propTable)

# subsetEqual <- simResults %>%
#   filter(compare == "Equal")
# View(subsetEqual)
# 
# subsetMethod1Better <- simResults %>%
#   filter(compare == "1 Better")
# View(subsetMethod1Better)
# 
# subsetMethod3Better <- simResults %>%
#   filter(compare == "3 Better")
# View(subsetMethod3Better)