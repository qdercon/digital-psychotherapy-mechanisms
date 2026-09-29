
data {
  int nTimes;                                                         // number of time points (sessions)
  int nPpts;                                                          // number of participants
  array[nPpts] int condition;                                         // intervention condition for each participant
  array[nPpts] int completed;                                         // whether participant completed time 2
  int nTrials_max;                                                    // maximum number of trials per participant per session
  array[nPpts, nTimes] int nT_ppts;                                   // actual number of trials per participant per session
  array[nPpts, nTimes, nTrials_max] real rew1;                        // array of option1 rewards, across participants, sessions, trials
  array[nPpts, nTimes, nTrials_max] real eff1;                        // array of option1 efforts, across participants, sessions, trials
  array[nPpts, nTimes, nTrials_max] real rew2;                        // array of option2 rewards, across participants, sessions, trials
  array[nPpts, nTimes, nTrials_max] real eff2;                        // array of option2 efforts, across participants, sessions, trials
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> choice01;  // chosen option (0=route 1, 1=route 2, -1 missing) per ppts, session, trial
}

parameters {
  // group-level correlation matrix (cholesky factor for faster computation)
  cholesky_factor_corr[nTimes] R_chol_rewSens;
  cholesky_factor_corr[nTimes] R_chol_effSens;

  // group-level parameters
  // means
  vector[nTimes] mu_rewSens;
  vector[nTimes] mu_effSens;
  // sds
  vector<lower=0>[nTimes] sigma_rewSens;           // sigma must be positive for normal dist
  vector<lower=0>[nTimes] sigma_effSens;           // sigma must be positive for normal dist

  // individual-level parameters  (raw/untransformed values)
  matrix[nTimes, nPpts] rewSens_raw;
  matrix[nTimes, nPpts] effSens_raw;

  // group-level intervention effects
  real rewSens_int_c;
  real rewSens_int_nc;
  real effSens_int_c;
  real effSens_int_nc;

  // group-level effects of non-completion
  real rewSens_ncpl;
  real effSens_ncpl;
}

transformed parameters {
  // individual-level parameter off-sets (for non-centered parameterization)
  matrix[nTimes, nPpts] rewSens_tilde;
  matrix[nTimes, nPpts] effSens_tilde;

  // individual-level parameters
  matrix[nPpts, nTimes] rewSens;
  matrix[nPpts, nTimes] effSens;

  // construct individual offsets (for non-centered parameterization)
  // (to get the covariance Cholesky factor from the correlation Cholesky factor,
  // we need to multiply it by a diagonal matrix constructed from the variances
  // of the individual variates).
  rewSens_tilde = diag_pre_multiply(sigma_rewSens, R_chol_rewSens) * rewSens_raw;
  effSens_tilde = diag_pre_multiply(sigma_effSens, R_chol_effSens) * effSens_raw;

  // compute individual-level parameters from non-centered parameterizationfor (p in 1:nPpts) {
  for (p in 1:nPpts) {
    // baseline non-completion offset, both sessions (as in additive)
    if (completed[p] == 0) {
      rewSens_tilde[:, p] += rewSens_ncpl;
      effSens_tilde[:, p] += effSens_ncpl;
    }
    // time 1
    rewSens[p, 1] = -1 + Phi_approx(mu_rewSens[1] + rewSens_tilde[1, p]) * 4;
    effSens[p, 1] = -1 + Phi_approx(mu_effSens[1] + effSens_tilde[1, p]) * 10;

    // time 2 — intervention effect stratified by completion, applied to active arm
    real rew_int_p = 0;
    real eff_int_p = 0;
    if (condition[p] == 1) {          // active arm only (as in additive's condition gate)
      if (completed[p] == 1) {
        rew_int_p = rewSens_int_c;    // completer active-arm effect
        eff_int_p = effSens_int_c;
      } else {
        rew_int_p = rewSens_int_nc;   // non-completer active-arm effect (the relaxation)
        eff_int_p = effSens_int_nc;
      }
    }
    rewSens[p, 2] = -1 + Phi_approx(mu_rewSens[2] + rewSens_tilde[2, p] + rew_int_p) * 4;
    effSens[p, 2] = -1 + Phi_approx(mu_effSens[2] + effSens_tilde[2, p] + eff_int_p) * 10;
  }
}

model {
  // priors on cholesky factor of correlation (i.e. test-retest) matrix
  R_chol_rewSens ~ lkj_corr_cholesky(1);
  R_chol_effSens ~ lkj_corr_cholesky(1);

  // define priors on distribution of group-level parameters
  // means
  mu_rewSens ~ normal(0, 1);
  mu_effSens ~ normal(0, 1);
  // sds
  sigma_rewSens ~ cauchy(0,1);
  sigma_effSens ~ cauchy(0,1);

  // define priors on individual participant deviations from group parameter values
  to_vector(rewSens_raw) ~ normal(0, 1);
  to_vector(effSens_raw) ~ normal(0, 1);

  // intervention effects
  rewSens_int_c ~ normal(0, 1);
  effSens_int_c ~ normal(0, 1);
  rewSens_int_nc ~ normal(0, 1);
  effSens_int_nc ~ normal(0, 1);

  // non-completion baseline diff
  rewSens_ncpl ~ normal(0, 1);
  effSens_ncpl ~ normal(0, 1);

  // loop over observations
  for (p in 1:nPpts) {
    if (nT_ppts[p, 1] > 0) {
      vector[nT_ppts[p, 1]] v1_t1;
      vector[nT_ppts[p, 1]] v2_t1;
      
      v1_t1 = rewSens[p, 1] * to_vector(rew1[p, 1, 1:nT_ppts[p, 1]]) - effSens[p, 1] * to_vector(eff1[p, 1, 1:nT_ppts[p, 1]]);
      v2_t1 = rewSens[p, 1] * to_vector(rew2[p, 1, 1:nT_ppts[p, 1]]) - effSens[p, 1] * to_vector(eff2[p, 1, 1:nT_ppts[p, 1]]);
      choice01[p, 1, 1:nT_ppts[p, 1]] ~ bernoulli_logit(v2_t1 - v1_t1);
    }

    if (completed[p] == 1 && nT_ppts[p, 2] > 0) {
      vector[nT_ppts[p, 2]] v1_t2;
      vector[nT_ppts[p, 2]] v2_t2;

      v1_t2 = rewSens[p, 2] * to_vector(rew1[p, 2, 1:nT_ppts[p, 2]]) - effSens[p, 2] * to_vector(eff1[p, 2, 1:nT_ppts[p, 2]]);
      v2_t2 = rewSens[p, 2] * to_vector(rew2[p, 2, 1:nT_ppts[p, 2]]) - effSens[p, 2] * to_vector(eff2[p, 2, 1:nT_ppts[p, 2]]);
      choice01[p, 2, 1:nT_ppts[p, 2]] ~ bernoulli_logit(v2_t2 - v1_t2);
    }
  }
}

generated quantities {
  // test-retest correlations
  corr_matrix[2] R_rewSens;
  corr_matrix[2] R_effSens;

  // log-likelihoods
  array[nPpts, nTimes] real sum_log_lik_ppt;
  vector[nPpts] log_lik;

	// reconstruct correlation matrix from cholesky factors
  R_effSens = multiply_lower_tri_self_transpose(R_chol_effSens);
  R_rewSens = multiply_lower_tri_self_transpose(R_chol_rewSens);

  // changes in parameters
  real delta_rewSens = mu_rewSens[2] - mu_rewSens[1];
  real delta_effSens = mu_effSens[2] - mu_effSens[1];

  // difference between interventions, in completers and non-completers
  real rewSens_int = rewSens_int_c - rewSens_int_nc;
  real effSens_int = effSens_int_c - effSens_int_nc;

  // changes in parameters (individual-level)
  vector[nPpts] delta_rewSens_p;
  vector[nPpts] delta_effSens_p;

  // posterior predictions
  array[nPpts, nTimes, nTrials_max] int y_pred;

  // generate posterior predictions and get likelihood
  for (p in 1:nPpts) {
    delta_rewSens_p[p] = rewSens[p, 2] - rewSens[p, 1];
    delta_effSens_p[p] = effSens[p, 2] - effSens[p, 1];
    sum_log_lik_ppt[p, 1] = 0;
    sum_log_lik_ppt[p, 2] = 0;

    if (nT_ppts[p, 1] > 0) {
      vector[nT_ppts[p, 1]] v1_t1;
      vector[nT_ppts[p, 1]] v2_t1;
      
      v1_t1 = rewSens[p, 1] * to_vector(rew1[p, 1, 1:nT_ppts[p, 1]]) - effSens[p, 1] * to_vector(eff1[p, 1, 1:nT_ppts[p, 1]]);
      v2_t1 = rewSens[p, 1] * to_vector(rew2[p, 1, 1:nT_ppts[p, 1]]) - effSens[p, 1] * to_vector(eff2[p, 1, 1:nT_ppts[p, 1]]);
      
      for (t in 1:nT_ppts[p, 1]) {
        sum_log_lik_ppt[p, 1] += bernoulli_logit_lpmf(choice01[p,1,t] | (v2_t1[t] - v1_t1[t]));
      }
      
      y_pred[p, 1, 1:nT_ppts[p, 1]] = bernoulli_rng(inv_logit(v2_t1 - v1_t1));
    }

    if (completed[p] == 1 && nT_ppts[p, 2] > 0) {
      vector[nT_ppts[p, 2]] v1_t2;
      vector[nT_ppts[p, 2]] v2_t2;

      v1_t2 = rewSens[p, 2] * to_vector(rew1[p, 2, 1:nT_ppts[p, 2]]) - effSens[p, 2] * to_vector(eff1[p, 2, 1:nT_ppts[p, 2]]);
      v2_t2 = rewSens[p, 2] * to_vector(rew2[p, 2, 1:nT_ppts[p, 2]]) - effSens[p, 2] * to_vector(eff2[p, 2, 1:nT_ppts[p, 2]]);

      y_pred[p, 2, 1:nT_ppts[p, 2]] = bernoulli_rng(inv_logit(v2_t2 - v1_t2));
      
      for (t in 1:nT_ppts[p, 2]) {
        sum_log_lik_ppt[p, 2] += bernoulli_logit_lpmf(choice01[p,2,t] | (v2_t2[t] - v1_t2[t])); 
      }
    }
    log_lik[p] = sum(sum_log_lik_ppt[p, :]);
  }
}
