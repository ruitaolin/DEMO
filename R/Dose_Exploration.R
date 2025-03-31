outcome = function(doses, dose.ind, coh.size, Y_B_mean.true, sigma2_B.true,
                   Y_T_mean.true, Y_R_mean.true,
                   delta1.true, delta2.true, delta3.true,
                   lambdaT.true, shape.true,
                   time_C = 24){

  Y_B = rnorm(coh.size, Y_B_mean.true[dose.ind], sqrt(sigma2_B.true))
  Y_T = as.numeric(runif(coh.size)<Y_T_mean.true[dose.ind])
  Y_R = as.numeric(runif(coh.size)<Y_R_mean.true[dose.ind])

  time_T = rweibull(coh.size, shape=shape.true, scale=(lambdaT.true[dose.ind]*exp(delta1.true*Y_T+delta2.true*Y_R+delta3.true*Y_B))^(-1/shape.true) )
  time_T = time_T
  Y_S = pmin(time_T, time_C)
  event = as.numeric(time_T <= time_C)

  list(Y_B=Y_B, Y_T=Y_T, Y_R=Y_R, Y_S=Y_S, event=event)
}

# BOIN
BOIN_sim = function(Y_B.true=Y_B_mean.true, sigma2_B.true=sigma2_B,
                    Y_T.true=Y_T,
                    Y_R.true=Y_R,
                    delta1.true=delta1,
                    delta2.true=delta2,
                    delta3.true=delta3,
                    lambdaT.true=lambdaT,
                    shape.true=shape,
                    time_C = time_C,
                    target_tox=target,
                    ncohort=10, cohortsize=3, startdose=1, n.earlystop=cohortsize*ncohort,
                    p.saf=0.6*target, p.tox=1.4*target, cutoff.eli=0.95, extrasafe=FALSE,
                    titration=F, offset=0.05, seed=1,
                    max_per_dose = 9,
                    monitor_cutoff_B = 0.30){
  set.seed(seed)
  ndose = length(Y_B.true)
  npts = ncohort*cohortsize

  if (cohortsize > 1) {
    temp = BOIN::get.boundary(target_tox, ncohort, cohortsize, n.earlystop=cohortsize*ncohort,
                              p.saf, p.tox, cutoff.eli, extrasafe)$full_boundary_tab
  }else{
    temp = BOIN::get.boundary(target_tox, ncohort, cohortsize, n.earlystop=cohortsize*ncohort,
                              p.saf, p.tox, cutoff.eli, extrasafe)$boundary_tab
  }

  b.e = temp[2, ]
  b.d = temp[3, ]
  b.elim = temp[4, ]

  d_tol = y_T_tol = y_R_tol = y_B_tol = NULL
  y_S_tol = event_tol = NULL
  y_T = rep(0, ndose)
  y_R = rep(0, ndose)
  n = rep(0, ndose)
  earlystop = 0
  d = startdose
  elimi = rep(0, ndose)
  ft = TRUE #flag used to determine whether or not to add cohortsize-1 patients to a dose for the first time when titration is triggered.
  if (titration){
    z = (runif(ndose) < p.true)
    if (sum(z) == 0) {
      d = ndose
      n[1:ndose] = 1
    }else{
      d = which(z == 1)[1]
      n[1:d] = 1
      y_T[d] = 1
    }
  }

  for (i in 1:ncohort){

    d_tol = c(d_tol, rep(d, cohortsize))

    if (titration && n[d] < cohortsize && ft){
      ft=FALSE
      y_T[d] = y_T[d] + sum(runif(cohortsize - 1) < Y_T.true[d])
      n[d] = n[d] + cohortsize - 1
    }else{
      Y_out = outcome(dose.ind = d, doses=doses, coh.size = cohortsize,
                      Y_B_mean.true=Y_B.true,
                      sigma2_B.true=sigma2_B.true,
                      Y_T_mean.true=Y_T.true,
                      Y_R_mean.true=Y_R.true,
                      delta1.true=delta1.true,
                      delta2.true=delta2.true,
                      delta3.true=delta3.true,
                      lambdaT.true=lambdaT.true,
                      time_C = time_C,
                      shape.true=shape.true)

      newcohort = Y_out$Y_T;
      y_T_tol = c(y_T_tol, newcohort)
      y_B_tol = c(y_B_tol, Y_out$Y_B)
      y_R_tol = c(y_R_tol, Y_out$Y_R)
      y_S_tol = c(y_S_tol, Y_out$Y_S)
      event_tol = c(event_tol, Y_out$event)

      if((sum(n)+cohortsize) >= npts){
        nremain = npts - sum(n);
        y_T[d] = y_T[d] + sum(newcohort[1:nremain]);
        n[d] = n[d] + nremain;

        y_R[d] = y_R[d] + sum(Y_out$Y_R[1:nremain]);
        break;
      }else{
        y_T[d] = y_T[d] + sum(newcohort);
        n[d] = n[d] + cohortsize;

        y_R[d] = y_R[d] + sum(Y_out$Y_R);
      }
    }

    if(!is.na(b.elim[n[d]])){
      if(y_T[d] >= b.elim[n[d]]){
        elimi[d:ndose] = 1
        if (d == 1){
          earlystop = 1
          break
        }
      }
      if(extrasafe){
        if (d == 1 && n[1] >= 3){
          if(1 - pbeta(target_tox, y_T[1] + 1, n[1] - y_T[1] + 1) > cutoff.eli - offset){
            earlystop = 1
            break
          }
        }
      }
    }

    if(i == ceiling(ncohort/2) ){
      dat_all = data.frame(d = d_tol,
                           Y_B = y_B_tol)
      tau_hat = tau_ms(dat=dat_all, monitor_cutoff_B=monitor_cutoff_B)
      no_eff = tau_hat - 1
      if(no_eff > 0) elimi[1:no_eff] = 1
    }

    if(n[d] >= max_per_dose) elimi[d:ndose] = 1

    if(n[d]>=n.earlystop && ((y_T[d]>b.e[n[d]] && y_T[d]<b.d[n[d]])||
                             (d==1 && y_T[d]>=b.d[n[d]]) ||
                             ((d==ndose||elimi[d+1]==1) && y_T[d]<=b.e[n[d]]))
    ) break;

    if (y_T[d] <= b.e[n[d]] && d != ndose){
      if (elimi[d + 1] == 0)
        d = d + 1
    }else if (y_T[d] >= b.d[n[d]] && d != 1){
      if (elimi[d - 1] == 0)
        d = d - 1
    }else{
      d = d
    }
  }
  return(list(summary=data.frame(y_T=y_T, Y_R=y_R, n=n),
              d_tol=d_tol, y_B_tol=y_B_tol,
              y_T_tol=y_T_tol, y_R_tol=y_R_tol,
              y_S_tol=y_S_tol, event_tol=event_tol))
}

pi_R_mean = function(gamma0, gamma1, gamma2, gamma3, dosage, YB){
  eta = gamma0 + gamma1*dosage + gamma2*dosage^2 + gamma3*YB
  pi_ORR = exp(eta)/(1+exp(eta))
  return(pi_ORR)
}

pi_T_mean = function(beta0, beta1, beta2, dosage, YB){
  eta = beta0 + exp(beta1)*dosage + exp(beta2)*YB
  pi_TOX = exp(eta)/(1+exp(eta))
  return(pi_TOX)
}

tau_ms = function(dat = dat_phaseI, monitor_cutoff_B = 0.3){

  max_dose = max(dat$d)
  post_prob = rep(0, max_dose)

  if(max_dose == 1) return(1)

  for(ind in 2:max_dose){

    y1 = dat$Y_B[dat$d < ind]; y2 = dat$Y_B[dat$d >= ind]
    y1_bar = mean(y1); y2_bar = mean(y2)
    n1 = length(y1); n2 = length(y2); n0 = 0.1
    alpha0 = beta0 = 0.01
    alpha1 = alpha0 + (n1+n2)/2
    beta1  = beta0 + 1/2*sum((y1-y1_bar)^2) + 1/2*sum((y2-y2_bar)^2) +
      n1*n0/(2*(n1+n0))*(y1_bar-0.0)^2 + n2*n0/(2*(n2+n0))*(y2_bar-0.5)^2
    post_prob[ind] = sqrt(1/(n1+n0)) * sqrt(1/(n2*n0)) * gamma(alpha1) / (beta1^alpha1)
  }

  post_prob = post_prob/sum(post_prob)
  tau_hat = which.max(post_prob)
  if(max(post_prob) <= monitor_cutoff_B) tau_hat = 1
  return(tau_hat)
}

post_YB_function = function(dat_alpha=post_alpha, d=1){
  post_YB = dat_alpha$post_alpha0 + dat_alpha$post_alpha1*d^(dat_alpha$post_alpha3) /
    (dat_alpha$post_alpha2^(dat_alpha$post_alpha3) + d^(dat_alpha$post_alpha3))
  mean(post_YB)
}

RMST_TRUE = function(lambda, alpha, delta1, delta2, delta3, YB, YT, YR, t_S=12){

  surv_fn = function(t){
    lambdaT = lambda*exp(delta1*YT + delta2*YR + delta3*YB)
    cum_hazard = lambdaT*t^(alpha)
    return(exp(-cum_hazard))
  }

  tmp = integrate(surv_fn, 0, t_S)
  return(tmp$value)
}

RMST_BUGS = function(n_gen=10, t_S=12, dose,
                     alpha0, alpha1, alpha2, alpha3, tau,
                     beta0, beta1, beta2,
                     gamma0, gamma1, gamma2, gamma3,
                     lambda=exp(mean(post_lambdaT)),
                     shape=mean(post_shape),
                     delta1=mean(post_delta1),
                     delta2=mean(post_delta2),
                     delta3=mean(post_delta3)){

  YB = YT = YE = rep(-1, n_gen * length(alpha0))

  # Generate data
  for(iter in 1:length(alpha0)){

    YB_mean = alpha0[iter] + alpha1[iter]*dose^(alpha3[iter]) / (alpha2[iter]^alpha3[iter] + dose^alpha3[iter])
    YB[ ((iter-1)*n_gen+1):(iter*n_gen) ] = rnorm(n_gen, YB_mean, 1/sqrt(tau[iter]))

    eta = beta0[iter] + exp(beta1[iter])*dose + exp(beta2[iter])*YB[ ((iter-1)*n_gen+1):(iter*n_gen) ]
    eta[eta>  100] = 100
    eta[eta< -100] = -100
    pi_T = exp(eta)/(1+exp(eta))
    YT[ ((iter-1)*n_gen+1):(iter*n_gen) ] = as.numeric(runif(n_gen)<pi_T)

    eta = gamma0[iter] + exp(gamma1[iter])*dose + gamma2[iter]*dose^2 + exp(gamma3[iter])*YB[ ((iter-1)*n_gen+1):(iter*n_gen) ]
    eta[eta>  100] = 100
    eta[eta< -100] = -100
    pi_E = exp(eta)/(1+exp(eta))
    YE[ ((iter-1)*n_gen+1):(iter*n_gen) ] = as.numeric(runif(n_gen)<pi_E)
  }

  RMST_fn = function(YT_i, YE_i, YB_i){

    surv_fn = function(t){
      lambdaT = lambda*exp(delta1*YT_i + delta2*YE_i + delta3*YB_i)
      cum_hazard = lambdaT*t^(shape)
      return(exp(-cum_hazard))
    }

    RMST_i = integrate(surv_fn, 0, t_S)$value
    return(RMST_i)
  }

  tmp = mapply(RMST_fn, YT_i=YT, YE_i=YE, YB_i=YB)

  return(tmp)
}
