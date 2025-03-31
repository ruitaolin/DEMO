#' A Demonstration of a Multi-Stage Dose-Finding Design
#'
#' @description
#' \code{DEMO_design()} runs a hypothetical multi-stage Phase I/II trial design:
#' \enumerate{
#'   \item Stage 1 uses a BOIN design for early dose exploration and toxicity control.
#'   \item Stage 2 refines the set of \dQuote{admissible} doses based on observed toxicity,
#'         biomarker data, and additional cohorts.
#'   \item Stage 3 performs further expansions on the remaining doses, incorporating a
#'         Weibull survival model, ultimately selecting an \dQuote{optimal} dose based
#'         on efficacy and safety criteria.
#' }
#'
#' @param doses A numeric vector of candidate dose levels. \strong{Default:} \code{c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85)}
#' @param Y_B_sim A numeric vector of true mean biomarker responses at each dose. \strong{Default:} \code{c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73)}
#' @param sigma2_B_sim A numeric scalar for the true variance of the biomarker. \strong{Default:} \code{1}
#' @param Y_T_sim A numeric vector of true short-term toxicity probabilities at each dose. \strong{Default:} \code{c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26)}
#' @param Y_R_sim A numeric vector of true intermediate-term response (ORR) probabilities at each dose. \strong{Default:} \code{c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47)}
#' @param delta1_sim,delta2_sim,delta3_sim Numeric scalars for parameters in the hazard function for time-to-event modeling. \strong{Defaults:} \code{3} and \code{-2}, respectively.
#' @param lambdaT_sim A numeric vector of baseline hazard parameters at each dose for survival modeling. \strong{Default:} \code{c(0.4, 0.1, 0.1, 0.3, 0.3, 0.3)}
#' @param shape_sim A numeric scalar for the Weibull shape parameter. \strong{Default:} \code{1.5}
#' @param censored_time A numeric scalar specifying the censoring time (in months, for example). \strong{Default:} \code{24}
#' @param RMST_followup A numeric scalar specifying the area under the survival curve up to a specific time point. \strong{Default:} \code{12}
#' @param cohort_stage1,cohort_stage2 Integers giving the number of cohorts in Stage 1 (BOIN) and Stage 2 (dose monitoring). \strong{Defaults:} \code{10} and \code{3}
#' @param cohortsize_stage1,cohortsize_stage2 Integers for patients per cohort in Stage 1 and Stage 2. \strong{Defaults:} \code{3} and \code{3}
#' @param M An integer for maximum number of patients at each dose across the three-stage trial. \strong{Default:} \code{24}
#' @param target_toxicity A numeric scalar specifying the target toxicity rate. \strong{Default:} \code{0.30}
#' @param min_acceptable_ORR A numeric scalar specifying the minimum acceptable ORR (overall response). \strong{Default:} \code{0.20}
#' @param min_acceptable_PFS A numeric scalar for minimum acceptable restricted mean survival time (RMST). \strong{Default:} \code{3.0}
#' @param c_B,c_T,c_R,c_S Numeric cutoffs for early biomarker monitoring (\code{c_B}), toxicity (\code{c_T}), ORR (\code{c_R}), and survival (\code{c_S}). \strong{Defaults:} \code{0.30, 0.60, 0.80, 0.90}
#' @param L1,L2 Integers controlling dose selection in Stage 2. \strong{Defaults:} \code{3} and \code{4}
#' @param kappa A numeric scalar for the \dQuote{efficacy plateau} threshold. \strong{Default:} \code{0.30}
#' @param seed An integer for the random seed. \strong{Default:} \code{1}
#'
#' @details
#' This function simulates patient outcomes based on input dose-toxicity, dose-biomarker,
#' and dose-response relationships. It then performs three stages of data-driven decision-making:
#'
#' \itemize{
#'   \item \strong{Stage 1 (BOIN):} Identifies dose ranges unlikely to exceed the \code{target_toxicity} and
#'         refines by biomarker cutoffs.
#'   \item \strong{Stage 2:} Gathers more data on remaining \dQuote{admissible} doses, updating posterior
#'         estimates of toxicity and efficacy, possibly further refining the set.
#'   \item \strong{Stage 3:} Expands only those doses that remain acceptable, fits a Weibull survival model,
#'         and ultimately selects an \emph{optimal therapeutic dose} (OTD) based on multiple endpoints.
#' }
#'
#' @return A \code{list} with elements:
#' \describe{
#'   \item{\code{trial}}{A data frame of patient-level outcomes (dose, biomarker, toxicity, survival, etc.).}
#'   \item{\code{N1, N2, N3}}{Vectors tracking the number of patients enrolled at each dose in Stages 1, 2, and 3.}
#'   \item{\code{OTD}}{The final selected Optimal Therapeutic Dose (or \code{0} if none is selected).}
#' }
#'
#' @examples
#' \dontrun{
#'   # Run with default parameters
#'   set.seed(123)
#'   result <- DEMO_design()
#'
#'   # Check final optimal dose
#'   print(result$OTD)
#'
#'   # Inspect overall trial data
#'   head(result$trial)
#' }
#'
#' @export

DEMO_design = function(doses = c(0.05, 0.10, 0.20, 0.45, 0.65, 0.85),
                       Y_B_sim = c(2.00, 2.01, 2.08, 2.76, 3.75, 4.73),
                       sigma2_B_sim = 1,
                       Y_T_sim = c(0.01, 0.02, 0.03, 0.06, 0.13, 0.26),
                       Y_R_sim = c(0.04, 0.05, 0.08, 0.20, 0.35, 0.47),
                       delta1_sim = 3, delta2_sim = -2, delta3_sim = -0.01,
                       lambdaT_sim = c(0.8, 0.6, 0.6, 0.25, 0.2, 0.1),
                       shape_sim = 1.5,
                       censored_time = 24,
                       RMST_followup = 12,
                       cohort_stage1 = 10,        # Number of cohorts for Stage 1 (dose exploration via BOIN)
                       cohortsize_stage1 = 3,
                       cohort_stage2 = 3,         # Number of additional cohorts for Stage 2 (dose monitoring)
                       cohortsize_stage2 = 3,
                       M = 24,                    # Total sample size for Stage 3 expansions
                       target_toxicity = 0.30,    # Maximum acceptable probability of short-term toxicity
                       min_acceptable_ORR = 0.20, # Minimum acceptable probability of intermediate-term response
                       min_acceptable_PFS = 3.0,  # Minimum acceptable restricted mean survival time
                       c_B = 0.30,    # Cutoff for early biomarker monitoring
                       c_T = 0.60,    # Cutoff for toxicity monitoring
                       c_R = 0.80,    # Cutoff for overall response rate (ORR) monitoring
                       c_S = 0.90,    # Cutoff for survival monitoring
                       L1 = 3,         # The number of best acceptable doses based on the posterior mean ORR
                       L2 = 4,         # The number of additional acceptable doses to account for plateau scenarios
                       kappa = 0.30,  # The tolerance parameter controlling the inclusion of additional acceptable doses
                       seed = 1){

  MTD_identify_file     = system.file("bugs", "MTD_identify.bugs", package = "DEMOdesign")
  stage12_file          = system.file("bugs", "stage12.bugs", package = "DEMOdesign")
  survival_weibull_file = system.file("bugs", "survival_weibull.bugs", package = "DEMOdesign")

  # Check that all vectors have the same length
  if (!(length(doses) == length(Y_B_sim) &&
        length(Y_B_sim) == length(Y_T_sim) &&
        length(Y_T_sim) == length(Y_R_sim) &&
        length(Y_R_sim) == length(lambdaT_sim))) {
    stop("Error: Length mismatch detected. Ensure that doses, Y_B_sim, Y_T_sim, Y_R_sim, and lambdaT_sim have the same length.")
  }

  # Check that sigma2_B_sim is non-negative
  if (sigma2_B_sim < 0) {
    stop("Error: sigma2_B_sim must be non-negative.")
  }

  # Set random seed for reproducibility
  set.seed(seed)

  # -- Basic Setup --
  n.dose = length(doses)

  # Track number of patients at each stage
  N1 = rep(0, n.dose)  # Stage 1
  N2 = rep(0, n.dose)  # Stage 2
  N3 = rep(0, n.dose)  # Stage 3

  # Initialize final selection
  OTD = 0   # Optimal Therapeutic Dose (OTD)

  # --------------------------------------------------------------------
  # STAGE 1: DOSE EXPLORATION USING BOIN
  # --------------------------------------------------------------------
  # Here we run a BOIN-based simulation to guide early dose exploration.
  trial = BOIN_sim(Y_B.true=Y_B_sim,
                   sigma2_B.true = sigma2_B_sim,
                   Y_T.true=Y_T_sim,
                   Y_R.true=Y_R_sim,
                   delta1.true = delta1_sim, delta2.true = delta2_sim, delta3.true = delta3_sim,
                   lambdaT.true = lambdaT_sim,
                   shape.true = shape_sim,
                   time_C = censored_time,
                   target_tox=target_toxicity,
                   p.saf=0.6*target_toxicity,
                   p.tox=1.4*target_toxicity,
                   monitor_cutoff_B = c_B,
                   ncohort=cohort_stage1,
                   cohortsize=cohortsize_stage1,
                   n.earlystop=9,
                   seed=seed)

  # Bayesian Model to Identify MTD (via posterior toxicity probabilities)
  jags_data = list(Npat = length(trial$d_tol),
                   d = doses[trial$d_tol],
                   NTox = trial$y_T_tol)

  jags_obj = R2jags::jags(model.file = MTD_identify_file, data = jags_data,
                          parameters.to.save = c("beta0", "beta1"),
                          n.chains = 3, n.iter = 7000, n.thin = 1,
                          progress.bar = "none", quiet = TRUE)

  post_beta0 = jags_obj$BUGSoutput$sims.matrix[,"beta0"]
  post_beta1 = jags_obj$BUGSoutput$sims.matrix[,"beta1"]

  # Compute posterior probability that pi_T (toxicity rate) > target_toxicity
  p_T_monitor = rep(NA, n.dose)
  for(j in 1:n.dose){
    tmp = pi_T_mean(beta0 = post_beta0,
                    beta1 = post_beta1,
                    beta2 = 0,            # not used at this stage
                    YB    = 0,            # no biomarker effect here
                    dosage = doses[j])
    p_T_monitor[j] = mean(tmp > target_toxicity)
  }

  # Remove the overly toxic doses
  MTD = max(which(p_T_monitor <= c_T))
  if(is.infinite(MTD)) MTD = 0

  # Number of DLT events & patients from Stage 1
  DLT_rate = sum(trial$y_T_tol)
  N1 = N1 + trial$summary[,3]

  # Stage 1: Identify change point for biomarker Y_B
  dat_all = data.frame(d = trial$d_tol,
                       Y_B = trial$y_B_tol,
                       y = trial$y_T_tol,
                       Y_R = trial$y_R_tol,
                       Y_S = trial$y_S_tol,
                       event = trial$event_tol)

  # Only evaluate biomarker data within [tau, MTD] range
  if(MTD > 0){
    dat_tmp = dat_all[dat_all$d <= MTD,]
    tau_hat = tau_ms(dat=dat_tmp, monitor_cutoff_B=c_B)
  }else{
    tau_hat = 0
  }

  # Define the "admissible set" at the end of Stage 1
  # Admissible doses are from tau_hat through MTD (if both are > 0)
  if(tau_hat <= MTD){
    admissible_set = tau_hat:MTD
    admissible_set = admissible_set[admissible_set!=0]
  }else{
    admissible_set = 0
  }

  # --------------------------------------------------------------------
  # STAGE 2: DOSE MONITORING
  # --------------------------------------------------------------------
  # We recruit additional cohorts to "admissible" doses to gather more data
  # about toxicity and biomarker. We keep updating the posterior estimates
  # after each small group is enrolled and re-check the MTD and biomarker cutoffs.
  N2_each_dose = cohort_stage2 * cohortsize_stage2
  interim = rep(N2_each_dose, n.dose)
  n2_stage = 0

  while(n2_stage < N2_each_dose & sum(admissible_set==0)==0 & length(admissible_set)!=0){

    n2_stage = n2_stage + 3
    index_dose = which(interim[admissible_set] > 0)
    if(length(index_dose)==0) break

    # Enroll 3 patients at each dose in the index_dose set
    for(ind_admst in seq_along(index_dose)){

      ind_tmp = index_dose[ind_admst]
      Y_out = outcome(dose.ind = admissible_set[ind_tmp], doses=doses,
                      coh.size = 3,
                      Y_B_mean.true=Y_B_sim,
                      sigma2_B.true=sigma2_B_sim,
                      Y_T_mean.true=Y_T_sim,
                      Y_R_mean.true=Y_R_sim,
                      delta1.true = delta1_sim, delta2.true = delta2_sim, delta3.true = delta3_sim,
                      lambdaT.true = lambdaT_sim,
                      shape.true = shape_sim,
                      time_C = censored_time)

      dat_tmp = cbind(d = rep(admissible_set[ind_tmp], 3),
                      data.frame(Y_out))
      colnames(dat_tmp) = colnames(dat_all)
      dat_all = rbind(dat_all, dat_tmp)
    }

    N2[admissible_set[index_dose]] = N2[admissible_set[index_dose]] + 3
    interim[admissible_set[index_dose]] = interim[admissible_set[index_dose]] - 3

    # Re-fit Bayesian model for updated toxicity estimates
    jags_data = list(Npat = length(dat_all$d),
                     d = doses[dat_all$d],
                     NTox = dat_all$y)

    jags_obj = R2jags::jags(model.file = MTD_identify_file, data = jags_data,
                            parameters.to.save = c("beta0", "beta1"),
                            n.chains = 3, n.iter = 7000, n.thin = 1,
                            progress.bar = "none", quiet = TRUE)

    post_beta0 = jags_obj$BUGSoutput$sims.matrix[,"beta0"]
    post_beta1 = jags_obj$BUGSoutput$sims.matrix[,"beta1"]

    # Compute updated posterior toxicity probabilities
    for(j in 1:n.dose){
      tmp = pi_T_mean(beta0 = post_beta0, beta1 = post_beta1, beta2 = 0,
                      YB = 0, dosage = doses[j])
      p_T_monitor[j] = mean(tmp > target_toxicity)
    }

    # Update MTD and admissible set
    MTD = min(max(admissible_set), max(which(p_T_monitor <= c_T)))
    if(is.infinite(MTD)) MTD = 0

    if(MTD > 0){
      dat_tmp = dat_all[dat_all$d <= MTD,]
      tau_hat = tau_ms(dat=dat_tmp, monitor_cutoff_B = c_B)
    }else{
      tau_hat = 0
    }

    # Update the new "admissible set"
    if(tau_hat <= MTD){
      admissible_set = tau_hat:MTD
      admissible_set = admissible_set[admissible_set!=0]
    }else{
      admissible_set = 0
    }
  } # end of Stage 2 while-loop

  # --------------------------------------------------------------------
  # Consolidate Stage 1 & 2 data, move to a combined Bayesian model
  # that includes toxicity, biomarker, and ORR for final expansions.
  # --------------------------------------------------------------------
  jags_data = list(Npat = dim(dat_all)[1],
                   YB = dat_all$Y_B,
                   d = doses[dat_all$d],
                   NEff = dat_all$Y_R,
                   NTox = dat_all$y)

  jags_obj = R2jags::jags(model.file = stage12_file, data = jags_data,
                          parameters.to.save = c("alpha0", "alpha1", "alpha2", "alpha3",
                                                 "beta0", "beta1", "beta2",
                                                 "gamma0", "gamma1", "gamma2", "gamma3"),
                          n.chains = 3, n.iter = 7000, n.thin = 1,
                          progress.bar = "none", quiet = TRUE)

  # Extract posterior draws for the emax-model parameters (biomarker)
  post_alpha0 = jags_obj$BUGSoutput$sims.matrix[,"alpha0"]
  post_alpha1 = jags_obj$BUGSoutput$sims.matrix[,"alpha1"]
  post_alpha2 = jags_obj$BUGSoutput$sims.matrix[,"alpha2"]
  post_alpha3 = jags_obj$BUGSoutput$sims.matrix[,"alpha3"]
  post_new = data.frame(post_alpha0=post_alpha0, post_alpha1=post_alpha1,
                        post_alpha2=post_alpha2, post_alpha3=post_alpha3)
  post_YB_mean = sapply(doses, function(x) post_YB_function(dat_alpha=post_new, d = x))

  # Extract posterior draws for toxicity
  post_beta0 = jags_obj$BUGSoutput$sims.matrix[,"beta0"]
  post_beta1 = jags_obj$BUGSoutput$sims.matrix[,"beta1"]
  post_beta2 = jags_obj$BUGSoutput$sims.matrix[,"beta2"]

  # Compute posterior p_T (prob. of toxicity) at each dose
  p_T_monitor = rep(NA, n.dose)
  post_T = rep(NA, n.dose)
  for(j in 1:n.dose){
    tmp = pi_T_mean(beta0 = post_beta0, beta1 = post_beta1, beta2 = post_beta2,
                    YB = post_YB_mean[j], dosage = doses[j])
    p_T_monitor[j] = mean(tmp > target_toxicity)
    post_T[j] = mean(tmp)
  }

  # Extract posterior draws for efficacy (ORR)
  post_gamma0 = jags_obj$BUGSoutput$sims.matrix[,"gamma0"]
  post_gamma1 = jags_obj$BUGSoutput$sims.matrix[,"gamma1"]
  post_gamma2 = jags_obj$BUGSoutput$sims.matrix[,"gamma2"]
  post_gamma3 = jags_obj$BUGSoutput$sims.matrix[,"gamma3"]

  p_R_monitor = rep(NA, n.dose)
  post_R = rep(NA, n.dose)
  for(j in 1:n.dose){
    tmp = pi_R_mean(gamma0 = post_gamma0, gamma1 = post_gamma1,
                    gamma2 = post_gamma2, gamma3 = post_gamma3,
                    YB = post_YB_mean[j], dosage = doses[j])
    p_R_monitor[j] = mean(tmp < min_acceptable_ORR)
    post_R[j] = mean(tmp)
  }

  # Eliminate doses if toxicity is too high or ORR is too poor
  elmi_monitor_T = which(p_T_monitor > c_T)
  elmi_monitor_R = which(p_R_monitor > c_R)

  # If any dose has unacceptable toxicity, eliminate doses from that dose upward
  if(length(elmi_monitor_T)!=0){
    elmi_monitor_T = min(elmi_monitor_T)
    admissible_set = admissible_set[!admissible_set %in% (elmi_monitor_T:n.dose)]
  }

  # If any dose has poor ORR, remove it from the set
  if(length(elmi_monitor_R)!=0){
    admissible_set = admissible_set[!admissible_set %in% elmi_monitor_R]
  }

  # If still more than L1 doses remain, keep only top L by ORR
  if(length(admissible_set) > L1){
    len_adm = length(admissible_set)
    original_set = admissible_set[which(len_adm - rank(post_R[admissible_set]) < L1)]
  }else{
    original_set = admissible_set
  }

  # Next: incorporate the “efficacy plateau” concept.
  # We find all doses whose ORR is at least (1-kappa)*max{ P(pi_R(d) > min_acceptable_ORR) }.
  max_pR = max(1-p_R_monitor)
  plateau_set = which((1-p_R_monitor) > max_pR*(1-kappa))
  plateau_set = plateau_set[plateau_set %in% admissible_set]
  if(length(plateau_set) > L2){
    len_adm = length(plateau_set)
    plateau_set = plateau_set[which(len_adm - rank(post_R[plateau_set]) < L2)]
  }

  # Combine both sets W = C1 \cup C2
  admissible_set = sort(unique(c(original_set, plateau_set)) )


  # --------------------------------------------------------------------
  # STAGE 3: DOSE OPTIMIZATION
  # --------------------------------------------------------------------
  # In Stage 3, we expand the sample size for the remaining "admissible" doses,
  # then refine final selection based on toxicity, ORR, and PFS (survival).
  dat_stageIII = dat_all[dat_all$d %in% admissible_set,]

  if(sum(admissible_set==0)==0 & length(admissible_set)!=0 & dim(dat_stageIII)[1]!=0){

    # ---------------------------
    # INTERIM 1 of Stage 3
    # ---------------------------
    # Expand half of the total leftover sample (M) in each admissible dose.
    dat_tmp = dat_all %>% dplyr::group_by(d) %>% dplyr::summarise(tol = dplyr::n())
    num_dose = dat_tmp$tol[admissible_set]
    num_expand = (M - num_dose) / 2; num_expand = ceiling(num_expand)
    index_dose = which(num_expand > 0)
    N3[admissible_set[index_dose]] = N3[admissible_set[index_dose]] +
      num_expand[index_dose]

    # Enroll these expansions
    for(ind_admst in seq_along(index_dose)){

      ind_tmp = index_dose[ind_admst]
      Y_out = outcome(dose.ind = admissible_set[ind_tmp], doses=doses,
                      coh.size = num_expand[ind_tmp],
                      Y_B_mean.true=Y_B_sim,
                      sigma2_B.true=sigma2_B_sim,
                      Y_T_mean.true=Y_T_sim,
                      Y_R_mean.true=Y_R_sim,
                      delta1.true = delta1_sim, delta2.true = delta2_sim, delta3.true = delta3_sim,
                      lambdaT.true = lambdaT_sim,
                      shape.true = shape_sim,
                      time_C = censored_time)

      dat_tmp = cbind(d = rep(admissible_set[ind_tmp], num_expand[ind_tmp]),
                      data.frame(Y_out))

      colnames(dat_tmp) = colnames(dat_all)
      dat_all = rbind(dat_all, dat_tmp)
    }

    # Re-estimate posterior for toxicity & ORR
    dat_all = dat_all[order(dat_all$d),]
    jags_data = list(Npat = dim(dat_all)[1],
                     YB = dat_all$Y_B,
                     d = doses[dat_all$d],
                     NEff = dat_all$Y_R,
                     NTox = dat_all$y)

    jags_obj = R2jags::jags(model.file = stage12_file, data = jags_data,
                            parameters.to.save = c("alpha0", "alpha1", "alpha2", "alpha3", "tau",
                                                   "beta0", "beta1", "beta2",
                                                   "gamma0", "gamma1", "gamma2", "gamma3"),
                            n.chains = 3, n.iter = 7000, n.thin = 1,
                            progress.bar = "none", quiet = T)

    post_alpha0 = jags_obj$BUGSoutput$sims.matrix[,"alpha0"]
    post_alpha1 = jags_obj$BUGSoutput$sims.matrix[,"alpha1"]
    post_alpha2 = jags_obj$BUGSoutput$sims.matrix[,"alpha2"]
    post_alpha3 = jags_obj$BUGSoutput$sims.matrix[,"alpha3"]
    post_tau = jags_obj$BUGSoutput$sims.matrix[,"tau"]

    # posterior for coefficients of toxicity
    post_beta0 = jags_obj$BUGSoutput$sims.matrix[,"beta0"]
    post_beta1 = jags_obj$BUGSoutput$sims.matrix[,"beta1"]
    post_beta2 = jags_obj$BUGSoutput$sims.matrix[,"beta2"]

    # posterior for coefficients of overall response rate (ORR)
    post_gamma0 = jags_obj$BUGSoutput$sims.matrix[,"gamma0"]
    post_gamma1 = jags_obj$BUGSoutput$sims.matrix[,"gamma1"]
    post_gamma2 = jags_obj$BUGSoutput$sims.matrix[,"gamma2"]
    post_gamma3 = jags_obj$BUGSoutput$sims.matrix[,"gamma3"]

    post_new = data.frame(post_alpha0=post_alpha0, post_alpha1=post_alpha1,
                          post_alpha2=post_alpha2, post_alpha3=post_alpha3)
    post_YB_mean = sapply(doses, function(x) post_YB_function(dat_alpha=post_new, d = x))

    # Updated toxicity posterior
    p_T_monitor = rep(NA, n.dose)
    post_T = rep(NA, n.dose)
    for(j in 1:n.dose){
      tmp = pi_T_mean(beta0 = post_beta0, beta1 = post_beta1, beta2 = post_beta2,
                      YB = post_YB_mean[j], dosage = doses[j])
      p_T_monitor[j] = mean(tmp > target_toxicity)
      post_T[j] = mean(tmp)
    }

    # Updated ORR posterior
    p_R_monitor = rep(NA, n.dose)
    post_R = rep(NA, n.dose)
    for(j in 1:n.dose){
      tmp = pi_R_mean(gamma0 = post_gamma0, gamma1 = post_gamma1,
                      gamma2 = post_gamma2, gamma3 = post_gamma3,
                      YB = post_YB_mean[j], dosage = doses[j])
      p_R_monitor[j] = mean(tmp < min_acceptable_ORR)
      post_R[j] = mean(tmp)
    }

    # Now incorporate survival (PFS) via Weibull model
    YS_jags = dat_all$Y_S
    is.na(YS_jags) = (dat_all$event==0)
    is.censored = 1 - dat_all$event
    YS_cen = dat_all$Y_S + (1-dat_all$event)

    dat_tmp = dat_all %>% group_by(d) %>% summarise(tol = n())

    jags_data = list(Npat = dim(dat_all)[1],
                     d = rep(dat_tmp$d, dat_tmp$tol),
                     J = length(unique(dat_all$d)),
                     YB = dat_all$Y_B,
                     YT = dat_all$y,
                     YE = dat_all$Y_R,
                     YS = YS_jags,
                     YS_cen = YS_cen,
                     is.censored = is.censored)

    jags_obj = R2jags::jags(model.file = survival_weibull_file, data = jags_data,
                            parameters.to.save = c("delta1", "delta2", "delta3", "lambdaT", "alpha"),
                            n.chains = 3, n.iter = 7000, n.thin = 1,
                            progress.bar = "none", quiet = T)

    # Posterior for Weibull parameters
    post_shape  = jags_obj$BUGSoutput$sims.matrix[,"alpha"]
    post_delta1 = jags_obj$BUGSoutput$sims.matrix[,"delta1"]
    post_delta2 = jags_obj$BUGSoutput$sims.matrix[,"delta2"]
    post_delta3 = jags_obj$BUGSoutput$sims.matrix[,"delta3"]

    p_S_monitor = rep(0, n.dose)
    n_tmp = NCOL(jags_obj$BUGSoutput$sims.matrix) - 5
    for(j in 1:n_tmp){
      post_lambdaT = jags_obj$BUGSoutput$sims.matrix[,5+j]
      RMST_ind = RMST_BUGS(n_gen=10, dose=doses[j], t_S = RMST_followup,
                           alpha0=post_alpha0, alpha1=post_alpha1,
                           alpha2=post_alpha2, alpha3=post_alpha3, tau=post_tau,
                           beta0=post_beta0, beta1=post_beta1, beta2=post_beta2,
                           gamma0=post_gamma0, gamma1=post_gamma1, gamma2=post_gamma2,
                           gamma3=post_gamma3,
                           lambda=exp(mean(post_lambdaT)),
                           shape=mean(post_shape),
                           delta1=mean(post_delta1),
                           delta2=mean(post_delta2),
                           delta3=mean(post_delta3))

      p_S_monitor[j] = mean(RMST_ind < min_acceptable_PFS)
    }

    # Eliminate doses based on updated cutoff criteria
    elmi_monitor_T = which(p_T_monitor > c_T)
    elmi_monitor_R = which(p_R_monitor > c_R)
    elmi_monitor_S = which(p_S_monitor > c_S)

    if(length(elmi_monitor_T)!=0){
      elmi_monitor_T = max(elmi_monitor_T)
      admissible_set = admissible_set[!admissible_set %in% (elmi_monitor_T:n.dose)]
    }
    if(length(elmi_monitor_R)!=0){
      admissible_set = admissible_set[!admissible_set %in% elmi_monitor_R]
    }
    if(length(elmi_monitor_S)!=0){
      admissible_set = admissible_set[!admissible_set %in% elmi_monitor_S]
    }
    # Interim End

    # ---------------------------
    # INTERIM 2 of Stage 3
    # ---------------------------
    if(length(admissible_set) != 0){

      dat_tmp = dat_all %>% group_by(d) %>% summarise(tol = n())
      num_dose = dat_tmp$tol[admissible_set]
      num_expand = M - num_dose
      index_dose = which(num_expand > 0)
      N3[admissible_set[index_dose]] = N3[admissible_set[index_dose]] +
        num_expand[index_dose]

      # Recruit these expansions
      for(ind_admst in seq_along(index_dose)){

        ind_tmp = index_dose[ind_admst]
        Y_out = outcome(dose.ind = admissible_set[ind_tmp], doses=doses,
                        coh.size = num_expand[ind_tmp],
                        Y_B_mean.true=Y_B_sim,
                        sigma2_B.true= sigma2_B_sim,
                        Y_T_mean.true=Y_T_sim,
                        Y_R_mean.true=Y_R_sim,
                        delta1.true = delta1_sim, delta2.true = delta2_sim, delta3.true = delta3_sim,
                        lambdaT.true = lambdaT_sim,
                        shape.true = shape_sim,
                        time_C = censored_time)

        dat_tmp = cbind(d = rep(admissible_set[ind_tmp], num_expand[ind_tmp]),
                        data.frame(Y_out))

        colnames(dat_tmp) = colnames(dat_all)
        dat_all = rbind(dat_all, dat_tmp)
      }

      # Re-estimate ptox and peff
      dat_all = dat_all[order(dat_all$d),]
      jags_data = list(Npat = dim(dat_all)[1],
                       YB = dat_all$Y_B,
                       d = doses[dat_all$d],
                       NEff = dat_all$Y_R,
                       NTox = dat_all$y)

      jags_obj = R2jags::jags(model.file = stage12_file, data = jags_data,
                              parameters.to.save = c("alpha0", "alpha1", "alpha2", "alpha3", "tau",
                                                     "beta0", "beta1", "beta2",
                                                     "gamma0", "gamma1", "gamma2", "gamma3"),
                              n.chains = 3, n.iter = 7000, n.thin = 1,
                              progress.bar = "none", quiet = T)

      post_alpha0 = jags_obj$BUGSoutput$sims.matrix[,"alpha0"]
      post_alpha1 = jags_obj$BUGSoutput$sims.matrix[,"alpha1"]
      post_alpha2 = jags_obj$BUGSoutput$sims.matrix[,"alpha2"]
      post_alpha3 = jags_obj$BUGSoutput$sims.matrix[,"alpha3"]
      post_tau = jags_obj$BUGSoutput$sims.matrix[,"tau"]

      # posterior for coefficients of toxicity
      post_beta0 = jags_obj$BUGSoutput$sims.matrix[,"beta0"]
      post_beta1 = jags_obj$BUGSoutput$sims.matrix[,"beta1"]
      post_beta2 = jags_obj$BUGSoutput$sims.matrix[,"beta2"]

      # posterior for coefficients of efficacy
      post_gamma0 = jags_obj$BUGSoutput$sims.matrix[,"gamma0"]
      post_gamma1 = jags_obj$BUGSoutput$sims.matrix[,"gamma1"]
      post_gamma2 = jags_obj$BUGSoutput$sims.matrix[,"gamma2"]
      post_gamma3 = jags_obj$BUGSoutput$sims.matrix[,"gamma3"]

      post_new = data.frame(post_alpha0=post_alpha0, post_alpha1=post_alpha1,
                            post_alpha2=post_alpha2, post_alpha3=post_alpha3)
      post_YB_mean = sapply(doses, function(x) post_YB_function(dat_alpha=post_new, d = x))

      p_T_monitor = rep(NA, n.dose)
      post_T = rep(NA, n.dose)
      for(j in 1:n.dose){
        tmp = pi_T_mean(beta0 = post_beta0, beta1 = post_beta1, beta2 = post_beta2,
                        YB = post_YB_mean[j], dosage = doses[j])
        p_T_monitor[j] = mean(tmp > target_toxicity)
        post_T[j] = mean(tmp)
      }

      p_R_monitor = rep(NA, n.dose)
      post_R = rep(NA, n.dose)
      for(j in 1:n.dose){
        tmp = pi_R_mean(gamma0 = post_gamma0, gamma1 = post_gamma1,
                        gamma2 = post_gamma2, gamma3 = post_gamma3,
                        YB = post_YB_mean[j], dosage = doses[j])
        p_R_monitor[j] = mean(tmp < min_acceptable_ORR)
        post_R[j] = mean(tmp)
      }

      # Re-estimate Survival
      YS_jags = dat_all$Y_S
      is.na(YS_jags) = (dat_all$event==0)
      is.censored = 1 - dat_all$event
      YS_cen = dat_all$Y_S + (1-dat_all$event)

      dat_tmp = dat_all %>% group_by(d) %>% summarise(tol = n())

      jags_data = list(Npat = dim(dat_all)[1],
                       d = rep(dat_tmp$d, dat_tmp$tol),
                       J = length(unique(dat_all$d)),
                       YB = dat_all$Y_B,
                       YT = dat_all$y,
                       YE = dat_all$Y_R,
                       YS = YS_jags,
                       YS_cen = YS_cen,
                       is.censored = is.censored)

      jags_obj = R2jags::jags(model.file = survival_weibull_file, data = jags_data,
                              parameters.to.save = c("delta1", "delta2", "delta3", "lambdaT", "alpha"),
                              n.chains = 3, n.iter = 7000, n.thin = 1,
                              progress.bar = "none", quiet = T)

      # posterior for coefficients of weibull distribution
      post_shape  = jags_obj$BUGSoutput$sims.matrix[,"alpha"]
      post_delta1 = jags_obj$BUGSoutput$sims.matrix[,"delta1"]
      post_delta2 = jags_obj$BUGSoutput$sims.matrix[,"delta2"]
      post_delta3 = jags_obj$BUGSoutput$sims.matrix[,"delta3"]

      RMST_tol = rep(-1, n.dose)
      p_S_monitor = rep(0, n.dose)
      n_tmp = NCOL(jags_obj$BUGSoutput$sims.matrix) - 5
      for(j in 1:n_tmp){
        post_lambdaT = jags_obj$BUGSoutput$sims.matrix[,5+j]
        RMST_ind = RMST_BUGS(n_gen=10, dose=doses[j], t_S = RMST_followup,
                             alpha0=post_alpha0, alpha1=post_alpha1,
                             alpha2=post_alpha2, alpha3=post_alpha3, tau=post_tau,
                             beta0=post_beta0, beta1=post_beta1, beta2=post_beta2,
                             gamma0=post_gamma0, gamma1=post_gamma1, gamma2=post_gamma2,
                             gamma3=post_gamma3,
                             lambda=exp(mean(post_lambdaT)),
                             shape=mean(post_shape),
                             delta1=mean(post_delta1),
                             delta2=mean(post_delta2),
                             delta3=mean(post_delta3))

        p_S_monitor[j] = mean(RMST_ind < min_acceptable_PFS)
        RMST_tol[j] = mean(RMST_ind)
      }

      elmi_monitor_T = which(p_T_monitor > c_T)
      elmi_monitor_R = which(p_R_monitor > c_R)
      elmi_monitor_S = which(p_S_monitor > c_S)

      if(length(elmi_monitor_T)!=0){
        elmi_monitor_T = max(elmi_monitor_T)
        admissible_set = admissible_set[!admissible_set %in% (elmi_monitor_T:n.dose)]
      }
      if(length(elmi_monitor_R)!=0){
        admissible_set = admissible_set[!admissible_set %in% elmi_monitor_R]
      }
      if(length(elmi_monitor_S)!=0){
        admissible_set = admissible_set[!admissible_set %in% elmi_monitor_S]
      }

      # Determine Optimal Therapeutic Dose (OTD)
      if(length(admissible_set)!=0){
        OTD = admissible_set[which.max(RMST_tol[admissible_set])]
      }
    }
  }

  dat_all = dat_all[order(dat_all$d),]
  rownames(dat_all) = NULL

  result_list = list(
    trial = dat_all,
    N1    = N1,
    N2    = N2,
    N3    = N3,
    OTD   = OTD
  )

  return(result_list)
}
