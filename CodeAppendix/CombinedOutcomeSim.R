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

# simCases <- data.frame(beta1_input = c(0.1, 2, 5, 10),
#                        beta2_input = c(0.1, 2, 5, 10),
#                        m_input     = c(10, 50, 100, 300),
#                        rho01_input = c(0.01, .025, 0.05, 0.1),
#                        rho02_input = c(0.01, .025, 0.05, 0.1),
#                        varY1_input = c(0.2, 0.5, 1, 5),
#                        varY2_input = c(0.2, 0.5, 1, 5),
#                        rho1_input  = c(0.01, 0.02, 0.03, 0.04),
#                        rho2_input  = c(0.05, 0.07, 0.08, 0.1))

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
  
  method1_power_dat <- data.frame(method1_power = method1_power,
                                  beta1_input = beta1_input,
                                  beta2_input = beta2_input,
                                  m_input     = m_input,
                                  rho01_input = rho01_input,
                                  rho02_input = rho02_input,
                                  varY1_input = varY1_input,
                                  varY2_input = varY2_input,
                                  rho1_input  = rho1_input,
                                  rho2_input  = rho2_input
                                    )

  # method1_power_dat <- c(method1_power,
  #                        beta1_input, beta2_input, m_input,
  #                        rho01_input, rho02_input,
  #                        varY1_input, varY2_input,
  #                        rho1_input, rho2_input
  # )
  
  return(method1_power_dat)
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
  
  method3_power_dat <- data.frame(method3_power = method3_power,
                                  beta1_input = beta1_input,
                                  beta2_input = beta2_input,
                                  m_input     = m_input,
                                  rho01_input = rho01_input,
                                  rho02_input = rho02_input,
                                  varY1_input = varY1_input,
                                  varY2_input = varY2_input,
                                  rho1_input  = rho1_input,
                                  rho2_input  = rho2_input
  )
  
  return(method3_power_dat)
}

# Run Sim Cases ----------------------------------------------------------------

# method1_results <- setNames(data.frame(matrix(ncol = 10, nrow = 0)),
#                  c("method1_power", "beta1_input", "beta2_input",
#                    "m_input", "rho01_input", "rho02_input",
#                    "varY1_input", "varY2_input",
#                    "rho1_input", "rho2_input"))
# 
# for(i in 1:nrow(allCases)){
#   row_i <- method1_fun(beta1_input = allCases$beta1_input[i],
#                        beta2_input = allCases$beta2_input[i],
#                        m_input     = allCases$m_input[i],
#                        rho01_input = allCases$rho01_input[i],
#                        rho02_input = allCases$rho02_input[i],
#                        varY1_input = allCases$varY1_input[i],
#                        varY2_input = allCases$varY2_input[i],
#                        rho1_input  = allCases$rho1_input[i],
#                        rho2_input  = allCases$rho2_input[i],
#                        K_input     = 10,
#                        alpha_input = 0.05
#                        )
#   method1_results <- bind_rows(method1_results, row_i)
# }
# 
# method3_results <- setNames(data.frame(matrix(ncol = 10, nrow = 0)),
#                             c("method3_power", "beta1_input", "beta2_input",
#                               "m_input", "rho01_input", "rho02_input",
#                               "varY1_input", "varY2_input",
#                               "rho1_input", "rho2_input"))
# 
# for(i in 1:nrow(allCases)){
#   row_i <- method3_fun(beta1_input = allCases$beta1_input[i],
#                        beta2_input = allCases$beta2_input[i],
#                        m_input     = allCases$m_input[i],
#                        rho01_input = allCases$rho01_input[i],
#                        rho02_input = allCases$rho02_input[i],
#                        varY1_input = allCases$varY1_input[i],
#                        varY2_input = allCases$varY2_input[i],
#                        rho1_input  = allCases$rho1_input[i],
#                        rho2_input  = allCases$rho2_input[i],
#                        K_input     = 10,
#                        alpha_input = 0.05
#   )
#   method3_results <- bind_rows(method3_results, row_i)
# }
# 
# LoopResults <- full_join(method1_results, method3_results,
#                         by = c("beta1_input", "beta2_input",
#                                "m_input", "rho01_input", "rho02_input",
#                                "varY1_input", "varY2_input",
#                                "rho1_input", "rho2_input"))
# 
# write.csv(LoopResults, "LoopResults.csv", row.names = FALSE)

# Summary of Results -----------------------------------------------------------

allResults <- read.csv("./LoopResults.csv") %>%
  rename(method2_power = method1_power) %>%
  mutate(method2_power = round(method2_power*100, 4),
         method3_power = round(method3_power*100, 4),
         method3minus2 = method3_power - method2_power,
         compare = ifelse(method2_power == method3_power, "Equal", 
                          ifelse(method2_power > method3_power, "2 Better",
                                 ifelse(method3_power > method2_power, "3 Better", 
                                        "CHECK"))),
         VIF1 = 1 + (m_input-1)*rho01_input,
         VIF2 = 1 + (m_input-1)*rho02_input,
         VIF12 = rho2_input + (m_input-1)*rho1_input,
         beta_ratio = beta1_input/beta2_input,
         var_ratio = varY1_input/varY2_input,
         rho_ratio = rho01_input/rho02_input,
         VIF_ratio = VIF1/VIF2,
         
         beta_ratio_1 = ifelse(beta_ratio == 1, 1, 0),
         var_ratio_1 = ifelse(var_ratio == 1, 1, 0),
         rho_ratio_1 = ifelse(rho01_input == rho02_input, 1, 0),
         VIF1_equal_VIF2 = ifelse(VIF1 == VIF2, 1, 0),
         
         rho_var_same = ifelse(varY1_input == varY2_input & rho01_input == rho02_input, 1, 0)
         ) %>%
  mutate(compare = ifelse(rho_var_same == 1 & abs(method3minus2) < 0.3, "Equal", compare)) %>%
  relocate(compare, rho_var_same, method3minus2, method2_power, method3_power)

scenarioSummary <- allResults %>%
  group_by(var_ratio_1, rho_ratio_1) %>%
  dplyr::summarise(n = n())

summaryTable <- allResults %>%
  group_by(compare, var_ratio_1, rho_ratio_1) %>%
  dplyr::summarise(n = n()) %>%
  mutate(Var = ifelse(var_ratio_1 == 1, "VarE", "VarNE"),
         Rho = ifelse(rho_ratio_1 == 1, "RhoE", "RhoNE")
         ) %>%
  unite("case", c("Var", "Rho"), sep = "_") %>%
  ungroup() %>%
  dplyr::select(compare, n, case) %>%
  spread(case, n) %>%
  replace(is.na(.), 0) %>%
  add_row(compare = "Total", summarise(., across(where(is.numeric), sum))) %>%
  mutate(Total = rowSums(across(where(is.numeric))))

histDat <- allResults %>%
  dplyr::select(method2_power, method3_power) %>%
  gather("Method", "Power") %>%
  mutate(Method = ifelse(Method == "method2_power", "Method 2", "Method 3"))

ggplot(histDat, aes(Power)) + facet_wrap(~Method) + geom_histogram()

nrow(filter(allResults, method2_power > 99.9999)) / nrow(allResults)
nrow(filter(allResults, method3_power > 99.9999)) / nrow(allResults)


allResultsLess100 <- allResults %>%
  filter(method2_power < 99.999, method3_power < 99.999) %>%
  mutate(Scenario = ifelse(var_ratio_1 == 1 & rho_ratio_1 == 1, "Equal Variances and ICC's",
                           ifelse(var_ratio_1 == 1 & rho_ratio_1 == 0, "Equal Variances Only",
                                  ifelse(var_ratio_1 == 0 & rho_ratio_1 == 1, "Equal ICC's Only",
                                         ifelse(var_ratio_1 == 0 & rho_ratio_1 == 0, "Unequal Variances and ICC's", "CHECK")))))
  
summaryTableLess100 <- allResultsLess100 %>%
  group_by(compare, var_ratio_1, rho_ratio_1) %>%
  dplyr::summarise(n = n()) %>%
  mutate(Var = ifelse(var_ratio_1 == 1, "VarE", "VarNE"),
         Rho = ifelse(rho_ratio_1 == 1, "RhoE", "RhoNE")
  ) %>%
  unite("case", c("Var", "Rho"), sep = "_") %>%
  ungroup() %>%
  dplyr::select(compare, n, case) %>%
  spread(case, n) %>%
  replace(is.na(.), 0) %>%
  add_row(compare = "Total", summarise(., across(where(is.numeric), sum))) %>%
  mutate(Total = rowSums(across(where(is.numeric))))

histDatLess100 <- allResultsLess100 %>%
  dplyr::select(method2_power, method3_power) %>%
  gather("Method", "Power") %>%
  mutate(Method = ifelse(Method == "method2_power", "Method 2", "Method 3"))

ggplot(histDatLess100, aes(Power)) + facet_wrap(~Method) + geom_histogram()

# Percent Equivalent
filter(summaryTableLess100, compare == "Equal")$Total/filter(summaryTableLess100, compare == "Total")$Total

# Percent Method 2 better than Method 3
filter(summaryTableLess100, compare == "2 Better")$Total/filter(summaryTableLess100, compare == "Total")$Total

# Percent Method 3 better than Method 2
filter(summaryTableLess100, compare == "3 Better")$Total/filter(summaryTableLess100, compare == "Total")$Total


# How many input scenarios have VIF1 = VIF2, Var(Y1) = Var(Y2), and B1 = B2
nrow(filter(allResults, VIF1 == VIF2,
            varY1_input == varY2_input,
            beta1_input == beta2_input)); nrow(allCases)

# Does beta ratio matter?
beta_rat <- allResultsLess100 %>%
  group_by(beta_ratio_1, compare) %>%
  dplyr::summarize(n = n()) %>%
  ungroup() %>%
  group_by(compare) %>%
  spread(beta_ratio_1, n) %>%
  dplyr::select(compare, 'BetaE' = `1`, 'BetaNE' = `0`) %>%
  ungroup() %>%
  add_row(compare = "Total", summarise(., across(where(is.numeric), sum))) %>%
  mutate(Total = rowSums(across(where(is.numeric))))


# Running a regression ---------------------------------------------------------

# View response variable
ggplot(data = allResultsLess100, aes(method3minus2)) + geom_histogram() + 
  xlab("Method 3 Power - Method 2 Power")

m1 <- lm(method3minus2 ~ var_ratio + rho_ratio, 
         data = allResultsLess100)
summary(m1)

m2 <- lm(method3minus2 ~ beta_ratio + var_ratio + rho_ratio, 
             data = allResultsLess100)
summary(m2)

m3 <- lm(method3minus2 ~ beta_ratio + var_ratio + rho_ratio + m_input, 
         data = allResultsLess100)
summary(m3)

m4 <- lm(method3minus2 ~ beta_ratio + var_ratio + VIF_ratio, 
         data = allResultsLess100)
summary(m4)

m5 <- lm(method3minus2 ~ beta_ratio + var_ratio + rho_ratio + m_input + 
           VIF_ratio + rho1_input + rho2_input, 
         data = allResultsLess100)
summary(m5)

predDat <- allResultsLess100 %>%
  mutate(pred1 = predict(m1, newdata = allResultsLess100),
         pred2 = predict(m2, newdata = allResultsLess100),
         pred3 = predict(m3, newdata = allResultsLess100),
         pred4 = predict(m4, newdata = allResultsLess100),
         pred5 = predict(m5, newdata = allResultsLess100)
         ) %>%
  relocate(pred1, pred2, pred3, pred4, pred5)

backward <- step(m5, direction = 'backward')
backward$coefficients
summary(backward)

#ggplot(dat = allResultsLess100, aes(method3minus2)) + geom_histogram() + 
#  facet_wrap(~Scenario)

# Table of difference statistics per scenario
View(mosaic::favstats(data = allResultsLess100, method3minus2 ~ Scenario))
#ggplot(dat = allResultsLess100, aes(var_ratio, method3minus2)) + geom_point()


filter(allResultsLess100, 
       
       varY1_input == 10,
       varY2_input == 10,
       
       rho01_input == 0.1,
       rho02_input == 0.5,
       
       beta1_input == 5,
       beta2_input == 0.1,
       
       m_input     == 15,
       
       rho1_input  == 0.05,
       rho2_input  == 0.2

       )$method3minus2


m3 <- lm(method3minus2 ~ var_ratio + VIF_ratio, 
         data = allResultsLess100)
summary(m3)

m4 <- lm(method3minus2 ~ var_ratio + VIF_ratio, 
         data = allResultsLess100)
summary(m4)


backward <- step(m.full, direction = 'backward')
backward$coefficients
summary(backward)



m6 <- lm(method3minus2 ~ beta_ratio*var_ratio + beta_ratio*rho_ratio + 
           beta_ratio + var_ratio + rho_ratio + m_input + rho1_input + rho2_input,
         data = allResultsLess100)
summary(m6)

m.red <- lm(diff_3minus1 ~ beta_ratio + var_ratio + 
              VIF1 + VIF2 + VIF12 + beta_ratio*var_ratio, 
         data = modelData)
summary(m.red)




# Separating cases and looking at distributions

method3better100 <- allResultsLess100 %>%
  filter(compare == "3 Better")
round(nrow(method3better100)/nrow(allResultsLess100), 4)*100
mosaic::tally(method3better100$Scenario)

method2better100 <- allResultsLess100 %>%
  filter(compare == "2 Better")
round(nrow(method2better100)/nrow(allResultsLess100), 4)*100
mosaic::tally(method2better100$Scenario)

equal100 <- allResultsLess100 %>%
  filter(compare == "Equal")
round(nrow(equal100)/nrow(allResultsLess100), 4)*100
mosaic::tally(equal100$Scenario)


outliers100 <- allResultsLess100 %>%
  filter(compare == "Equal", rho_var_same == 0)

mosaic::tally(method3better100$rho_var_same)
mosaic::tally(method2better100$rho_var_same)
mosaic::tally(equal100$rho_var_same)



# Graphs and Tables ------------------------------------------------------------

# Creating graphing dataset
graphData <- modelData %>%
  mutate(beta_ratio_1 = ifelse(beta_ratio == 1, "Beta_Ratio_1", "Beta_Ratio_Not1"),
         var_ratio_1 = ifelse(var_ratio == 1, "Var_Ratio_1", "Var_Ratio_Not1"),
         VIF_ratio = VIF1/VIF2,
         VIF_ratio_1 = ifelse(VIF_ratio == 1, "VIF_Ratio_1", "VIF_Ratio_Not1"))

# Betas ratio
ggplot(data = filter(graphData), 
       aes(y = diff_3minus1, x = beta_ratio)) + facet_wrap(~var_ratio_1*VIF_ratio_1) + 
  geom_point()

# Vars
ggplot(data = filter(graphData), 
       aes(y = diff_3minus1, x = var_ratio)) + facet_wrap(~beta_ratio_1*VIF_ratio_1) +
  geom_point()

# VIFs
ggplot(data = filter(graphData), 
       aes(y = diff_3minus1, x = VIF_ratio)) + facet_wrap(~beta_ratio_1*var_ratio_1) +
  geom_point()


mosaic::tally(data = simResults, ~compare)

# Histogram of beta proportions
ggplot(data = simResults, aes(x = beta_ratio)) + 
  geom_histogram(bins = 20) + facet_wrap(~compare)

# Histogram of variance proportions
ggplot(data = simResults, aes(x = var_ratio)) + 
  geom_histogram(bins = 20) + facet_wrap(~compare)

# Histogram of VIF1
ggplot(data = simResults, aes(x = VIF1)) + 
  geom_histogram(bins = 20) + facet_wrap(~compare)
