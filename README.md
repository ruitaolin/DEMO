
# **💊 DEMO Design: Dose Exploration, Monitoring, and Optimization using Biological and Clinical Outcomes**

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R->=4.2-blue)](https://cran.r-project.org/)
[![JAGS Required](https://img.shields.io/badge/JAGS-Required-red)](http://mcmc-jags.sourceforge.net/)
<!-- badges: end -->

## Overview

**DEMOdesign** is an R package implementing the **Dose Exploration, Monitoring, and Optimization (DEMO)** design—a Bayesian adaptive Phase I/II trial framework tailored for oncology dose-finding, particularly for immunotherapies and targeted therapies. It selects an **Optimal Therapeutic Dose (OTD)** that maximizes **Restricted Mean Survival Time (RMST)** while ensuring safety, biological activity, and clinical efficacy across three stages:

1. **Sequential Exploration**: Identifies safe, biologically active doses using biomarkers ($Y_B$) and toxicity ($Y_T$).
2. **Randomization and Monitoring**: Refines doses based on response ($Y_R$) and additional safety data.
3. **Randomized Optimization**: Selects the OTD using long-term survival ($Y_S$) via a Weibull model.

This package supports:

> Yang, C.-H., Thall, P. F., & Lin, R. (2024). **DEMO: Dose Exploration, Monitoring, and Optimization Using Biological and Clinical Outcomes**. *Annals of Applied Statistics*. (Under Review).

---

## Features

- **Three-Stage Bayesian Design**:
  - **Stage 1**: BOIN-based exploration with biomarker and toxicity screening.
  - **Stage 2**: Monitors admissible doses for efficacy and safety.
  - **Stage 3**: Optimizes dose selection with RMST.
- **Endpoints**: Integrates 
  - $Y_B$ (biomarker), 
  - $Y_T$ (toxicity), 
  - $Y_R$ (response), 
  - $Y_S$ (survival).

---

## 📦 Installation

Install from GitHub:

``` r
remotes::install_github("cyang728/DEMOdesign")
```

---

## Functions

### `DEMO_design()`

#### Description
Implements a three-stage Bayesian adaptive design for Phase 1-2 oncology trials. It integrates biomarker ($Y_B$), toxicity ($Y_T$), response ($Y_R$), and survival ($Y_S$) outcomes to identify an Optimal Therapeutic Dose (OTD). The OTD is defined as the dose that maximizes the restricted mean survival time (RMST) while meeting criteria for safety, clinical effectiveness, and biological activity.

#### Inputs

- `seed`: Integer. Random seed for reproducibility.
- `doses`: Numeric vector. Dose levels to evaluate.
- `Y_B_sim`: Numeric vector. Mean biomarker values for each dose.
- `sigma2_B_sim`: Numeric. Variance of the biomarker outcome.
- `Y_T_sim`: Numeric vector. Toxicity probabilities for each dose.
- `Y_R_sim`: Numeric vector. Response probabilities for each dose.
- `lambdaT_sim`: Numeric vector. Weibull scale parameters for survival at each dose.
- `shape_sim`: Numeric. Weibull shape parameter.
- `delta1_sim`: Numeric. Effect of toxicity on survival time.
- `delta2_sim`: Numeric. Effect of response on survival time.
- `delta2_sim`: Numeric. Effect of biomarker on survival time.
- `censored_time`: Numeric. Administrative censoring time (e.g., 24 months).
- `RMST_followup`: Numeric. Time horizon for RMST calculation (e.g., 12 months).  
- `cohort_stage1`: Integer. Number of cohorts in Stage 1 (Exploration).  
- `cohort_stage2`: Integer. Number of cohorts in Stage 2 (Monitoring).  
- `cohortsize_stage1`: Integer. Number of patients per cohort in Stage 1.  
- `cohortsize_stage2`: Integer. Number of patients per cohort in Stage 2.  
- `M`: Integer. Maximum number of patients at each dose across the three-stage trial.
- `target_toxicity`: Numeric. Target short-term toxicity rate.  
- `min_acceptable_ORR`: Numeric. Minimum acceptable overall response rate.  
- `min_acceptable_PFS`: Numeric. Minimum acceptable RMST.  
- `c_B`: Numeric. Decision cutoff for early biomarker activity.  
- `c_T`: Numeric. Decision cutoff for short-term toxicity.  
- `c_R`: Numeric. Decision cutoff for tumor response.  
- `c_S`: Numeric. Decision cutoff for survival (RMST).  
- `L1`: Integer. Number of best acceptable doses based on the posterior mean response rate.
- `L2`: Integer. Number of additional acceptable doses to account for plateau scenarios where some doses have similar response rates.
- `kappa`: Numeric. Plateau threshold for including additional acceptable doses for Stage 3.

#### Outputs

A list with the following elements:

- `OTD`: The selected Optimal Therapeutic Dose.
- `trial`: A data frame with patient-level trial data, including assigned dose, biomarker, toxicity, response, and survival outcomes.

### `calibrate_cB_fn()`

#### Description
Calibrates the biomarker cutoff parameter ($c_B$) for the DEMO design by running multiple BOIN-based simulations across various scenarios. It evaluates false active rate (FAR) and false inactive rate (FIR) to select an optimal $c_B$ that minimizes a composite error index, ensuring robust biomarker-based dose selection.

#### Inputs
- `cB_candidate`: Numeric vector. Candidate $c_B$ values to test (e.g., 0.2 to 0.9).
- `ntrial`: Integer. Number of simulation trials (e.g., 1000).
- `doses`: Numeric vector. Dose levels to evaluate.
- `Y_B_sim`: List of numeric vectors. Mean biomarker responses for each simulation scenario.
- `Y_T_sim`: List of numeric vectors. Toxicity probabilities for each simulation scenario.
- `Y_R_sim`: List of numeric vectors. Response probabilities for each simulation scenario.
- `lambdaT_sim`: List of numeric vectors. Weibull scale parameters for each simulation scenario.
- `sigma2_B_sim`: Numeric vector. Biomarker variances for each simulation scenario.
- `delta1_sim`, `delta2_sim`, `delta3_sim`: Numeric vectors. Effects of toxicity, response, and biomarker on survival for each scenario.
- `shape_sim`: Numeric vector. Weibull shape parameters for each simulation scenario.
- `time_C`: Numeric. Administrative censoring time (e.g., 24 months).
- `target`: Numeric. Target toxicity rate (e.g., 0.3).
- `cohortsize`: Integer. Patients per cohort (e.g., 3).
- `ncohort`: Integer. Number of cohorts in BOIN design (e.g., 10).

#### Outputs
A list with the following elements:
- `composite_index`: Matrix of composite error indices (sqrt(FAR^2 + FIR^2)) for each $c_B$ candidate across scenarios.
- `best_cB`: The optimal $c_B$ value minimizing the mean composite error index.


### `calibrate_cT_cR_cS_kappa_fn()`

#### Description
Calibrates the monitoring cutoff parameters ($c_T$, $c_R$, $c_S$) and plateau parameter ($\kappa$) for the DEMO design by running simulations across multiple scenarios. It evaluates the probability of correct selection (PCS) of the OTD and average sample size, selecting optimal values that maximize a composite index (PCS / sqrt(sample size)).


#### Inputs
- `cT_candidate`: Numeric vector. Candidate $c_T$ values (e.g., c(0.5, 0.6, 0.7)).
- `cT_candidate`: Numeric vector. Candidate $c_R$ values (e.g., c(0.6, 0.7, 0.8)).
- `cS_candidate`: Numeric vector. Candidate $c_S$ values (e.g., c(0.7, 0.8, 0.9)).
- `kappa_candidate`: Numeric vector. Candidate $\kappa$ values (e.g., c(0.2, 0.3, 0.4)).
- `monitor_cutoff_B`: Numeric. Fixed biomarker cutoff ($c_B$, e.g., 0.4).
- `ntrial`: Integer. Number of simulation trials (e.g., 1000).
- `true_OTD_sim`: Numeric vector. True OTDs for each scenario (e.g., c(6, 5, 3, 2)).
- `doses`: Numeric vector. Dose levels (e.g., c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85)).
- `Y_B_sim`: List of numeric vectors. Mean biomarker responses per scenario.
- `Y_T_sim`: List of numeric vectors. Toxicity probabilities per scenario.
- `Y_R_sim`: List of numeric vectors. Response probabilities per scenario.
- `lambdaT_sim`: List of numeric vectors. Weibull scale parameters per scenario.
- `sigma2_B_sim`: Numeric vector. Biomarker variances per scenario.
- `delta1_sim`, `delta2_sim`, `delta3_sim`: Numeric vectors. Effects on survival per scenario.
- `shape_sim`: Numeric vector. Weibull shape parameters per scenario.
- `censored_time`: Numeric. Censoring time (e.g., 24 months).
- `RMST_followup`: Numeric. RMST follow-up time (e.g., 12 months).
- `cohort_stage1`, `cohortsize_stage1`: Integers. Stage 1 cohorts and size.
- `cohort_stage2`, `cohortsize_stage2`: Integers. Stage 2 cohorts and size.
- `M`: Integer. Maximum sample size.
- `target_toxicity`: Numeric. Target toxicity rate (e.g., 0.3).
- `min_acceptable_ORR`: Numeric. Minimum response rate (e.g., 0.2).
- `min_acceptable_PFS`: Numeric. Minimum RMST (e.g., 3.0).
- `L1`, `L2`: Integers. Dose selection parameters.


#### Outputs
A list with:
- `best_cT`: Optimal $c_T$ value.
- `best_cR`: Optimal $c_R$ value.
- `best_cS`: Optimal $c_S$ value.
- `best_kappa`: Optimal $\kappa$ value.

---

## 🚀 Example

The following example demonstrates how to use DEMOdesign to simulate a Bayesian adaptive dose-finding trial.

``` r
library(DEMOdesign)

# Calibrate cB with six doses and four scenarios
calib_result <- calibrate_cB_fn(
cB_candidate = 2:9/10,
  ntrial = 1000,
  doses = c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85),
  Y_B_sim = list(
    c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73),
    c(2.01, 2.09, 2.24, 4.46, 5.29, 5.95),
    c(2.24, 4.00, 5.77, 5.99, 6.00, 6.00),
    c(5.04, 5.83, 5.98, 6.00, 6.00, 6.00)
  ),
  Y_T_sim = list(
    c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26),
    c(0.01, 0.03, 0.04, 0.07, 0.14, 0.28),
    c(0.01, 0.02, 0.05, 0.10, 0.27, 0.55),
    c(0.01, 0.06, 0.18, 0.29, 0.51, 0.54)
  ),
  Y_R_sim = list(
    c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47),
    c(0.04, 0.06, 0.09, 0.23, 0.37, 0.44),
    c(0.07, 0.14, 0.32, 0.41, 0.42, 0.44),
    c(0.22, 0.38, 0.41, 0.42, 0.44, 0.45)
  ),
  lambdaT_sim = list(
    c(0.80, 0.60, 0.60, 0.25, 0.20, 0.10),
    c(0.80, 0.40, 0.30, 0.30, 0.20, 0.35),
    c(0.40, 0.10, 0.10, 0.30, 0.30, 0.30),
    c(0.12, 0.10, 0.20, 0.30, 0.30, 0.30)
  ),
  sigma2_B_sim = c(1, 1, 1, 1),
  delta1_sim = c(3, 3, 3, 3),
  delta2_sim = c(-2, -2, -2, -2),
  delta3_sim = c(0, 0, 0, 0),
  shape_sim = c(1.5, 1.5, 1.5, 1.5),
  time_C = 24,
  target = 0.3,
  cohortsize = 3,
  ncohort = 10,
  max_per_dose = 9
)

# View calibrated cB
print(calib_result$best_cB)

# Calibrate monitoring parameters with reduced trials for demonstration
calib_result2 <- calibrate_cT_cR_cS_kappa_fn(
  ntrial = 100,  # Reduced for speed; use 1000+ in practice
  cT_candidate = c(0.5, 0.6),
  cR_candidate = c(0.6, 0.7),
  cS_candidate = c(0.7, 0.8),
  kappa_candidate = c(0.2, 0.3),
  monitor_cutoff_B = 0.4,
  true_OTD_sim = c(6, 5, 3, 2),
  doses = c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85),
  Y_B_sim = list(
    c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73),
    c(2.01, 2.09, 2.24, 4.46, 5.29, 5.95),
    c(2.24, 4.00, 5.77, 5.99, 6.00, 6.00),
    c(5.04, 5.83, 5.98, 6.00, 6.00, 6.00)
  ),
  Y_T_sim = list(
    c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26),
    c(0.01, 0.03, 0.04, 0.07, 0.14, 0.28),
    c(0.01, 0.02, 0.05, 0.10, 0.27, 0.55),
    c(0.01, 0.06, 0.18, 0.29, 0.51, 0.54)
  ),
  Y_R_sim = list(
    c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47),
    c(0.04, 0.06, 0.09, 0.23, 0.37, 0.44),
    c(0.07, 0.14, 0.32, 0.41, 0.42, 0.44),
    c(0.22, 0.38, 0.41, 0.42, 0.44, 0.45)
  ),
  lambdaT_sim = list(
    c(0.80, 0.60, 0.60, 0.25, 0.20, 0.10),
    c(0.80, 0.40, 0.30, 0.30, 0.20, 0.35),
    c(0.40, 0.10, 0.10, 0.30, 0.30, 0.30),
    c(0.12, 0.10, 0.20, 0.30, 0.30, 0.30)
  ),
  sigma2_B_sim = c(1, 1, 1, 1),
  delta1_sim = c(3, 3, 3, 3),
  delta2_sim = c(-2, -2, -2, -2),
  delta3_sim = c(0, 0, 0, 0),
  shape_sim = c(1.5, 1.5, 1.5, 1.5),
  censored_time = 24,
  RMST_followup = 12,
  cohort_stage1 = 10,
  cohortsize_stage1 = 3,
  cohort_stage2 = 3,
  cohortsize_stage2 = 3,
  M = 24,
  target_toxicity = 0.3,
  min_acceptable_ORR = 0.2,
  min_acceptable_PFS = 3.0,
  L1 = 3,
  L2 = 4
)

# View calibrated parameters
print(calib_result2$best_cT)
print(calib_result2$best_cR)
print(calib_result2$best_cS)
print(calib_result2$best_kappa)

# Run DEMO_design with example data
result <- DEMO_design(seed = 1,
                      doses = c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85),
                      Y_B_sim = c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73),  # biomarker means
                      sigma2_B_sim = 1, 
                      Y_T_sim = c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26),  # toxicity rates
                      Y_R_sim = c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47),  # response rates
                      lambdaT_sim = c(0.8, 0.6, 0.6, 0.25, 0.2, 0.1),   # scale for Weibull survival
                      shape_sim = 1.5,      # Weibull shape parameter
                      delta1_sim = 3,       # effect of toxicity on survival
                      delta2_sim = -2,      # effect of response on survival
                      delta3_sim = 0,       # effect of biomarker on survival
                      censored_time = 24    # administrative censoring time 
                      )

# View the optimal therapeutic dose (OTD)
print(result$OTD)

# View summary of patient-level trial data
head(result$trial)
```

### Example Output

#### From `calibrate_cB_fn()` 

##### `$best_cB`

``` 
[1] 0.4
```

#### From `calibrate_cT_cR_cS_kappa_fn()` 

##### `$best_cT`

``` 
[1] 0.6
```

##### `$best_cR`

``` 
[1] 0.8
```

##### `$best_cS`

``` 
[1] 0.9
```

##### `$best_kappa`

``` 
[1] 0.3
```

#### From `DEMO_design()` 

##### `$trial` (first 6 rows)

```txt
   d      Y_B y Y_R        Y_S      event
1  1 1.373546 0   0  0.6038467     1
2  1 2.183643 0   0  1.1267659     1
3  1 1.164371 0   0  0.4746055     1
4  2 2.004233 0   0  1.3692341     1
5  2 4.414653 0   0  0.3783198     1
6  2 2.773593 0   1  5.6059113     1
```

Each row represents a patient with the following outcomes: 
- `d`: Assigned dose 
- `Y_B`: Biomarker outcome 
- `y`: Toxicity outcome
- `Y_R`: Tumor response 
- `Y_S`: Observed survival time 
- `event`: Event indicator (1 = event, 0 = censored) 


##### `$N1`, `$N2`, `$N3` 

``` 
$N1 
[1] 3 3 3 3 6 9 

$N2 
[1] 0 0 0 0 6 9 

$N3 
[1] 0 0 0 0 12 6 
``` 
- `$N1`: Patients in the first stage 
- `$N2`: Patients in the second stage 
- `$N3`: Patients in the third stage


##### `$OTD` 
``` $[1] 6 ```
