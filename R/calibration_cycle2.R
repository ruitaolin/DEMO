#' Calibrate Multiple Monitoring Parameters for Multi-Stage Dose-Finding Design
#'
#' @description
#' \code{calibrate_cT_cR_cS_kappa_fn()} calibrates the monitoring cutoff parameters (cT, cR, cS) and plateau parameter (kappa) for a multi-stage dose-finding trial by:
#' \enumerate{
#'   \item Running multiple simulations using the DEMO design for dose exploration.
#'   \item Evaluating the probability of correct selection (PCS) and sample size across candidate parameter combinations.
#'   \item Returning the optimal cT, cR, cS, and kappa values based on a composite index.
#' }
#'
#' @param cT_candidate A numeric vector of candidate toxicity cutoff (cT) values to test. \strong{Default:} \code{c(0.5, 0.6, 0.7)}
#' @param cR_candidate A numeric vector of candidate response cutoff (cR) values to test. \strong{Default:} \code{c(0.6, 0.7, 0.8)}
#' @param cS_candidate A numeric vector of candidate survival cutoff (cS) values to test. \strong{Default:} \code{c(0.7, 0.8, 0.9)}
#' @param kappa_candidate A numeric vector of candidate plateau parameter (kappa) values to test. \strong{Default:} \code{c(0.2, 0.3, 0.4)}
#' @param monitor_cutoff_B A numeric scalar specifying the biomarker cutoff value. \strong{Default:} \code{0.4}
#' @param ntrial An integer specifying the number of simulation trials. \strong{Default:} \code{1000}
#' @param true_OTD_sim A numeric vector of true optimal therapeutic doses for each simulation scenario. \strong{Default:} \code{c(6, 5, 3, 2)}
#' @param doses A numeric vector of candidate dose levels. \strong{Default:} \code{c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85)}
#' @param Y_B_sim A list of numeric vectors of true mean biomarker responses for each simulation scenario. \strong{Default:} \code{list(c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73), c(2.01, 2.09, 2.24, 4.46, 5.29, 5.95), c(2.24, 4.00, 5.77, 5.99, 6.00, 6.00), c(5.04, 5.83, 5.98, 6.00, 6.00, 6.00))}
#' @param Y_T_sim A list of numeric vectors of true toxicity probabilities for each simulation scenario. \strong{Default:} \code{list(c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26), c(0.01, 0.03, 0.04, 0.07, 0.14, 0.28), c(0.01, 0.02, 0.05, 0.10, 0.27, 0.55), c(0.01, 0.06, 0.18, 0.29, 0.51, 0.54))}
#' @param Y_R_sim A list of numeric vectors of true response probabilities for each simulation scenario. \strong{Default:} \code{list(c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47), c(0.04, 0.06, 0.09, 0.23, 0.37, 0.44), c(0.07, 0.14, 0.32, 0.41, 0.42, 0.44), c(0.22, 0.38, 0.41, 0.42, 0.44, 0.45))}
#' @param lambdaT_sim A list of numeric vectors of baseline hazard scale parameters for each simulation scenario. \strong{Default:} \code{list(c(0.80, 0.60, 0.60, 0.25, 0.20, 0.10), c(0.80, 0.40, 0.30, 0.30, 0.20, 0.35), c(0.40, 0.10, 0.10, 0.30, 0.30, 0.30), c(0.12, 0.10, 0.20, 0.30, 0.30, 0.30))}
#' @param sigma2_B_sim A numeric vector of true biomarker variances for each simulation scenario. \strong{Default:} \code{c(1, 1, 1, 1)}
#' @param delta1_sim,delta2_sim,delta3_sim Numeric vectors of hazard function parameters for time-to-event modeling for each simulation scenario. \strong{Defaults:} \code{c(3, 3, 3, 3)}, \code{c(-2, -2, -2, -2)}, and \code{c(0, 0, 0, 0)}, respectively
#' @param shape_sim A numeric vector of Weibull shape parameters for each simulation scenario. \strong{Default:} \code{c(1.5, 1.5, 1.5, 1.5)}
#' @param censored_time A numeric scalar specifying the administrative censoring time (in months). \strong{Default:} \code{24}
#' @param RMST_followup A numeric scalar specifying the restricted mean survival time follow-up period (in months). \strong{Default:} \code{12}
#' @param cohort_stage1 An integer specifying the number of cohorts in stage 1. \strong{Default:} \code{10}
#' @param cohortsize_stage1 An integer specifying the cohort size in stage 1. \strong{Default:} \code{3}
#' @param cohort_stage2 An integer specifying the number of cohorts in stage 2. \strong{Default:} \code{3}
#' @param cohortsize_stage2 An integer specifying the cohort size in stage 2. \strong{Default:} \code{3}
#' @param M An integer specifying the maximum sample size. \strong{Default:} \code{24}
#' @param target_toxicity A numeric scalar specifying the target toxicity rate. \strong{Default:} \code{0.30}
#' @param min_acceptable_ORR A numeric scalar specifying the minimum acceptable overall response rate. \strong{Default:} \code{0.20}
#' @param min_acceptable_PFS A numeric scalar specifying the minimum acceptable progression-free survival (in months). \strong{Default:} \code{3.0}
#' @param L1,L2 Numeric scalars specifying look-back periods for toxicity and response evaluation (in months). \strong{Defaults:} \code{3} and \code{4}, respectively
#'
#' @details
#' This function evaluates combinations of candidate cT, cR, cS, and kappa values by simulating trials with the DEMO design,
#' calculating the probability of correct selection (PCS) of the optimal therapeutic dose (OTD) and the average sample size,
#' and determining optimal parameter values that maximize a composite index (PCS / sqrt(sample size)). It assumes the existence
#' of a helper function \code{DEMO_design()} for trial simulation.
#'
#' @return A \code{list} with elements:
#' \describe{
#'   \item{\code{best_cT}}{The cT value that maximizes the composite index.}
#'   \item{\code{best_cR}}{The cR value that maximizes the composite index.}
#'   \item{\code{best_cS}}{The cS value that maximizes the composite index.}
#'   \item{\code{best_kappa}}{The kappa value that maximizes the composite index.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Run with reduced parameters for demonstration
#'   set.seed(123)
#'   result <- calibrate_cT_cR_cS_kappa_fn(
#'     ntrial = 1000,
#'     cT_candidate = c(0.5, 0.6, 0.7),
#'     cR_candidate = c(0.6, 0.7, 0.8),
#'     cS_candidate = c(0.7, 0.8, 0.9),
#'     kappa_candidate = c(0.2, 0.3, 0.4)
#'   )
#'
#'   # Check optimal parameter values
#'   print(result$best_cT)
#'   print(result$best_cR)
#'   print(result$best_cS)
#'   print(result$best_kappa)
#' }
#'
#' @export

calibrate_cT_cR_cS_kappa_fn = function(
  cT_candidate = c(0.5, 0.6, 0.7),
  cR_candidate = c(0.6, 0.7, 0.8),
  cS_candidate = c(0.7, 0.8, 0.9),
  kappa_candidate = c(0.2, 0.3, 0.4),
  monitor_cutoff_B =  0.4,
  ntrial = 1000,
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
  target_toxicity = 0.30,
  min_acceptable_ORR = 0.20,
  min_acceptable_PFS = 3.0,
  L1 = 3,
  L2 = 4
){

  diff_monitor_par_mat = expand.grid(cT_candidate, cR_candidate,
                                     cS_candidate, kappa_candidate)

  PCS_scen = Ntol_scen = matrix(0, nrow = NROW(diff_monitor_par_mat), ncol = length(Y_T_sim))

  n_dose = length(doses)

  #num_cB = length(cB_candidate)
  #FP = FN = matrix(NA, ncol = num_cB, nrow = length(Y_B_sim))

  for(i_sim in 1:length(Y_T_sim)){

    true_OTD = true_OTD_sim[i_sim]
    Y_B = Y_B_sim[[i_sim]]
    Y_T = Y_T_sim[[i_sim]]
    Y_R = Y_R_sim[[i_sim]]
    lambdaT = lambdaT_sim[[i_sim]]

    sigma2_B = sigma2_B_sim[i_sim]
    delta1 = delta1_sim[i_sim]
    delta2 = delta2_sim[i_sim]
    delta3 = delta3_sim[i_sim]
    shape = shape_sim[i_sim]

    #results_total = matrix(NA, ncol = n_dose, nrow = length(cB_candidate))

    for(dif_par in 1:dim(diff_monitor_par_mat)[1]){

      PCS  = numeric(ntrial)
      Ntol = numeric(ntrial)

      monitor_cutoff_T = diff_monitor_par_mat[dif_par, 1]
      monitor_cutoff_R = diff_monitor_par_mat[dif_par, 2]
      monitor_cutoff_S = diff_monitor_par_mat[dif_par, 3]
      kappa_plateau = diff_monitor_par_mat[dif_par, 4]

      for(tr in 1:ntrial){

        set.seed(tr)

        # Run DEMO_design with example data
        result = DEMO_design(seed = tr,
                             doses = doses,
                             Y_B_sim = Y_B,                    # biomarker means
                             sigma2_B_sim = sigma2_B,
                             Y_T_sim = Y_T,                    # toxicity rates
                             Y_R_sim = Y_R,                    # response rates
                             lambdaT_sim = lambdaT,            # scale for Weibull survival
                             shape_sim = shape,                # Weibull shape parameter
                             delta1_sim = delta1,              # effect of toxicity on survival
                             delta2_sim = delta2,              # effect of response on survival
                             delta3_sim = delta3,              # effect of biomarker on survival
                             censored_time = censored_time,    # administrative censoring time
                             RMST_followup = RMST_followup,
                             cohort_stage1 = cohort_stage1,
                             cohortsize_stage1 = cohortsize_stage1,
                             cohort_stage2 = cohort_stage2,
                             cohortsize_stage2 = cohortsize_stage2,
                             M = M,
                             target_toxicity = target_toxicity,
                             min_acceptable_ORR = min_acceptable_ORR,
                             min_acceptable_PFS = min_acceptable_PFS,
                             c_B = monitor_cutoff_B,
                             c_T = monitor_cutoff_T,
                             c_R = monitor_cutoff_R,
                             c_S = monitor_cutoff_S,
                             kappa = kappa_plateau,
                             L1 = L1, L2 = L2
        )

        PCS[tr] = result$OTD
        Ntol[tr] = NROW(result$trial)

        #print(tr)
      }

      PCS_scen[dif_par,i_sim] = mean(PCS == true_OTD, na.rm=T)
      Ntol_scen[dif_par,i_sim]= mean(Ntol, na.rm=T)
    }
  }

  PCS_scen = apply(PCS_scen, 1, mean)
  Ntol_scen = apply(Ntol_scen, 1, mean)

  composite_index = PCS_scen / sqrt(Ntol_scen)
  best_cT = diff_monitor_par_mat[which.max(composite_index),1]
  best_cR = diff_monitor_par_mat[which.max(composite_index),2]
  best_cS = diff_monitor_par_mat[which.max(composite_index),3]
  best_kappa = diff_monitor_par_mat[which.max(composite_index),4]

  return(list(
    #composite_index = composite_index,
    best_cT = best_cT,
    best_cR = best_cR,
    best_cS = best_cS,
    best_kappa = best_kappa
  ))
}

#calibrate_cT_cR_cS_kappa_fn(ntrial = 2,
#                            cT_candidate = c(0.5),
#                            cR_candidate = c(0.6),
#                            cS_candidate = c(0.7, 0.9),
#                            kappa_candidate = c(0.2))

