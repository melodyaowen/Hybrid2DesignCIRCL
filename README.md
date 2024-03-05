# Hybrid2DesignCIRCL

## Overview

This respitory stores the code necessary for conducting power calculations for a hybrid type 2 study. These calculations use input parameters motivated by a study conducted by the Chicago Implementation Research Center. The code is supplementary material for the paper _The Design of Hybrid Type 2 Studies_<sup>1</sup>, published in [Insert Journal and DOI once publisehd].

The Chicago Implementation Research Center (CIRCL) is conducting a study to assess the impact of adapting a community-drive Kaiser implementation strategy bundle in Chicago to improve hypertension control.<sup>2,3</sup> The goal of this study is to evaluate the use of practice facilitation (PF), which refers to supportive and collaborative approaches used to help healthcare practices, clinics, and providers improve patient care quality. The use of PF in this study is to help community health centers better utilize the community-adapted Kaiser bundle. This bundle includes supportive care with the key compnent being "drop-in" blood pressure checks. The first primary outcome is the proportion of patients with controlled blood pressure (yes/no), and the second primary outcome is reach, which is defined as the proportion of patients, among those who were eligible, who received the Kaiser bundle (yes/no). This is a parallel 2-armed cluster-randomized trial (CRT), where each cluster is a clinic. 

## Required Input Parameters

_Table 1. Description and Values of Study Design Input Parameters from CIRCL_<sup>2,3</sup>
| Parameter | Statistical Notation | Variable Name in Code | Description | Value |
| ---                             | ---              | ---     | --- | --- |
| Statistical power               | $\pi$            | `power` | Probability of detecting a true effect under $H_A$ | 80% |
| Number of clusters              | $K$              | `K`     | Number of clusters in each treatment arm | 15 |
| Cluster size                    | $m$              | `m`     | Number of individuals in each cluster | 300 |
| Family-wise false positive rate | $\alpha$         | `alpha` | Probability of one or more Type I error(s) | 0.05 |
| Effect for $Y_1$                | $\beta_1^*$      | `beta1` | Estimated intervention effect on the first outcome ($Y_1$) | 10% |
| Effect for $Y_2$                | $\beta_2^*$      | `beta2` | Estimated intervention effect on the second outcome ($Y_2$) | 10% |
| Total variance of $Y_1$         | $\sigma_1^2$     | `varY1` | Total variance of the first outcome, $Y_1$ | 0.23 |
| Total variance of $Y_2$         | $\sigma_2^2$     | `varY2` | Total variance of the second outcome, $Y_2$ | 0.25 |
| Endpoint-specific ICC for $Y_1$ | $\rho_0^{(1)}$   | `rho01` | Correlation for $Y_1$ for two different individuals in the same cluster | 0.025 |
| Endpoint-specific ICC for $Y_2$ | $\rho_0^{(2)}$   | `rho02` | Correlation for $Y_2$ for two different individuals in the same cluster | 0.025 |
| Inter-subject between-endpoint ICC | $\rho_1^{(1,2)}$ | `rho1`  | Correlation between $Y_1$ and $Y_2$ for two different individuals in the same cluster | 0.01 |
| Intra-subject between-endpoint ICC | $\rho_2^{(1,2)}$ | `rho2`  | Correlation between $Y_1$ and $Y_2$ for the same individual | 0.05 |

## Code File Directory

All code necessary to run the calculations is available in the /CodeAppendix folder in the main directory. The R scripts are as follows:
- `CodeAppendix.R`: The main script. This script sources all necessary files and steps through each study design calculation for all five of the methods. Outputs a summary table at the end in order to visualize all the calculations completed.
- `RequiredPackages.R`: This script contains code that automatically checks for the required R packages to run the calculations. If the package is not installed, this script will install them. If they are installed, it will load the package. This is sourced and ran through the main script, `CodeAppendix.R`, so this script does not need to be opened and ran manually.
- `HelperFunctions.R`: This script contains helper functions for various sample size and power calculations. This is sourced and ran through the main script, `CodeAppendix.R`, so this script does not need to be opened and ran manually.
- `powerSampleCal_varCluster_ttest.R`: This script contains code for the conjunctive test<sup>4</sup> written by Siyun Yang, and available (here)[https://github.com/siyunyang/coprimary_CRT]. This is sourced in the main script, and does not need to be opened and ran manually.

## References
1. [Insert publication reference here]
2. Smith JD, Davis P, Kho AN. Community-Driven Health Solutions on Chicago's South Side. Stanf Soc Innov Rev. 2021; 19(3):A27-A9. DOI: 10.48558/85p7-3113
3. Kho A, Daumit GL, Truesdale KP, Brown A, Kilbourne AM, Ladapo J, et al. The National Heart Lung and Blood Institute Disparities Elimination through Coordinated Interventions to Prevent and Control Heart and Lung Disease Alliance. Health Serv Res. 2022; 57 Suppl 1(Suppl 1):20-31. DOI: 10.1111/1475-6773.13983
4. Yang S, Moerbeek M, Taljaard M, Li F. Power analysis for cluster randomized trials with continuous coprimary endpoints. Biometrics. 2022. DOI: 10.1111/biom.13692


