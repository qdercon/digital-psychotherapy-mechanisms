
data {
  int nTimes;                                                              // number of time points (sessions)
  int nPpts;                                                               // number of participants
  array[nPpts] int condition;                                              // intervention condition for each participant
  array[nPpts] int completed;                                              // whether participant completed time 2
  int nTrials_max;                                                         // maximum number of trials per participant per session
  array[nPpts, nTimes] int nT_ppts;                                        // actual number of trials per participant per session
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> internal_neg;   // responses for each participant per session
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> internal_pos;   // responses for each participant per session
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> global_neg;     // responses for each participant per session
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> global_pos;     // responses for each participant per session
}

parameters {
  // group-level correlation matrix (cholesky factor for faster computation)
  // across timepoints and two parameters
  cholesky_factor_corr[4] R_chol_theta_neg;
  cholesky_factor_corr[4] R_chol_theta_pos;

  // group-level parameters
  // means for each parameter and timepoint
  vector[nTimes] mu_internal_theta_neg;
  vector[nTimes] mu_internal_theta_pos;
  vector[nTimes] mu_global_theta_neg;
  vector[nTimes] mu_global_theta_pos;

  // sds for each parameter and timepoint
  vector<lower=0>[4] pars_sigma_neg;
  vector<lower=0>[4] pars_sigma_pos;

  // individual-level parameters (raw/untransformed values)
  matrix[4, nPpts] pars_pr_neg;
  matrix[4, nPpts] pars_pr_pos;

  // group-level effects of active intervention at t2 (group level)
  real theta_int_internal_neg;
  real theta_int_internal_pos;
  real theta_int_global_neg;
  real theta_int_global_pos;

  // group-level effects of non-completion (t1 and t2)
  real theta_ncompl_internal_neg;
  real theta_ncompl_internal_pos;
  real theta_ncompl_global_neg;
  real theta_ncompl_global_pos;
}

transformed parameters {
  // individual-level parameter off-sets (for non-centered parameterization)
  matrix[4, nPpts] pars_tilde_neg;
  matrix[4, nPpts] pars_tilde_pos;

  // individual-level parameters (transformed values)
  matrix[nPpts, nTimes] theta_internal_neg;
  matrix[nPpts, nTimes] theta_global_neg;
  matrix[nPpts, nTimes] theta_internal_pos;
  matrix[nPpts, nTimes] theta_global_pos;

  // construct individual offsets (for non-centered parameterization)
  // with potential correlation for parameters across timepoints
  // and potential correlation between internal-global parameter estimates
  pars_tilde_neg = diag_pre_multiply(pars_sigma_neg, R_chol_theta_neg) * pars_pr_neg;
  pars_tilde_pos = diag_pre_multiply(pars_sigma_pos, R_chol_theta_pos) * pars_pr_pos;

  // compute individual-level parameters from non-centered parameterization
  for (p in 1:nPpts) {
    // negative events
    // time 1
    theta_internal_neg[p, 1] = mu_internal_theta_neg[1] + pars_tilde_neg[1, p];
    theta_global_neg[p, 1]   = mu_global_theta_neg[1]   + pars_tilde_neg[2, p];

    // time 2
    if (condition[p] == 1) {
      // for active intervention participants
      theta_internal_neg[p, 2] = mu_internal_theta_neg[2] + pars_tilde_neg[3, p] + theta_int_internal_neg;
      theta_global_neg[p, 2]   = mu_global_theta_neg[2]   + pars_tilde_neg[4, p] + theta_int_global_neg;
    } else {
      // for control intervention participants
      theta_internal_neg[p, 2] = mu_internal_theta_neg[2] + pars_tilde_neg[3, p];
      theta_global_neg[p, 2]   = mu_global_theta_neg[2]   + pars_tilde_neg[4, p];
    }

    // positive events
    // time 1
    theta_internal_pos[p, 1] = mu_internal_theta_pos[1] + pars_tilde_pos[1, p];
    theta_global_pos[p, 1]   = mu_global_theta_pos[1]   + pars_tilde_pos[2, p];

    // time 2
    if (condition[p] == 1) {
      // for active intervention participants
      theta_internal_pos[p, 2] = mu_internal_theta_pos[2] + pars_tilde_pos[3, p] + theta_int_internal_pos;
      theta_global_pos[p, 2]   = mu_global_theta_pos[2]   + pars_tilde_pos[4, p] + theta_int_global_pos;
    } else {
      // for control intervention participants
      theta_internal_pos[p, 2] = mu_internal_theta_pos[2] + pars_tilde_pos[3, p];
      theta_global_pos[p, 2]   = mu_global_theta_pos[2]   + pars_tilde_pos[4, p];
    }

    // add non-completion effects (equivalent to adding to the offsets)
    if (completed[p] == 0) {
      theta_internal_neg[p, :] += theta_ncompl_internal_neg;
      theta_internal_pos[p, :] += theta_ncompl_internal_pos;
      theta_global_neg[p, :]   += theta_ncompl_global_neg;
      theta_global_pos[p, :]   += theta_ncompl_global_pos;
    }
  }
}

model {
  // uniform [0,1] priors on cholesky factor of correlation matrix
  R_chol_theta_neg ~ lkj_corr_cholesky(1);
  R_chol_theta_pos ~ lkj_corr_cholesky(1);

  // define priors on distribution of group-level parameters
  // means
  mu_internal_theta_neg ~ normal(0,1);
  mu_internal_theta_pos ~ normal(0,1);
  mu_global_theta_neg   ~ normal(0,1);
  mu_global_theta_pos   ~ normal(0,1);

  // sds of the individual deviations. these are declared <lower=0> but had no prior at all,
  // i.e. an implicit improper uniform on (0, inf); half-cauchy matches every rew-eff model
  pars_sigma_neg ~ cauchy(0,1);
  pars_sigma_pos ~ cauchy(0,1);

  // define priors on individual participant deviations from group parameter values
  to_vector(pars_pr_neg) ~ normal(0,1);
  to_vector(pars_pr_pos) ~ normal(0,1);

  // priors on group-level effects of active intervention on theta values at t2
  theta_int_internal_neg ~ normal(0,1);
  theta_int_internal_pos ~ normal(0,1);
  theta_int_global_neg   ~ normal(0,1);
  theta_int_global_pos   ~ normal(0,1);

  // priors on group-level effects of non-completion on theta values at t1 and t2
  theta_ncompl_internal_neg ~ normal(0,1);
  theta_ncompl_internal_pos ~ normal(0,1);
  theta_ncompl_global_neg   ~ normal(0,1);
  theta_ncompl_global_pos   ~ normal(0,1);

  // loop over observations
  for (p in 1:nPpts) {
    if (nT_ppts[p, 1] > 0) {
      target += bernoulli_logit_lpmf(internal_neg[p, 1, 1:nT_ppts[p, 1]] | theta_internal_neg[p, 1]);
      target += bernoulli_logit_lpmf(internal_pos[p, 1, 1:nT_ppts[p, 1]] | theta_internal_pos[p, 1]);
      target += bernoulli_logit_lpmf(global_neg[p, 1, 1:nT_ppts[p, 1]] | theta_global_neg[p, 1]);
      target += bernoulli_logit_lpmf(global_pos[p, 1, 1:nT_ppts[p, 1]] | theta_global_pos[p, 1]);
    }
    if (nT_ppts[p, 2] > 0) {
      target += bernoulli_logit_lpmf(internal_neg[p, 2, 1:nT_ppts[p, 2]] | theta_internal_neg[p, 2]);
      target += bernoulli_logit_lpmf(internal_pos[p, 2, 1:nT_ppts[p, 2]] | theta_internal_pos[p, 2]);
      target += bernoulli_logit_lpmf(global_neg[p, 2, 1:nT_ppts[p, 2]] | theta_global_neg[p, 2]);
      target += bernoulli_logit_lpmf(global_pos[p, 2, 1:nT_ppts[p, 2]] | theta_global_pos[p, 2]);
    }
  }
}

generated quantities {
  // test-retest correlations
  corr_matrix[4] R_theta_neg = multiply_lower_tri_self_transpose(R_chol_theta_neg);
  corr_matrix[4] R_theta_pos = multiply_lower_tri_self_transpose(R_chol_theta_pos);

  // success probability estimates for each individual
  matrix[nPpts, nTimes] p_internal_neg = inv_logit(theta_internal_neg);
  matrix[nPpts, nTimes] p_internal_pos = inv_logit(theta_internal_pos);
  matrix[nPpts, nTimes] p_global_neg = inv_logit(theta_global_neg);
  matrix[nPpts, nTimes] p_global_pos = inv_logit(theta_global_pos);

  // differences in group-level parameters
  real delta_internal_neg = mu_internal_theta_neg[2] - mu_internal_theta_neg[1];
  real delta_internal_pos = mu_internal_theta_pos[2] - mu_internal_theta_pos[1];
  real delta_global_neg   = mu_global_theta_neg[2] - mu_global_theta_neg[1];
  real delta_global_pos   = mu_global_theta_pos[2] - mu_global_theta_pos[1];

  // differences in individual-level parameters
  vector[nPpts] delta_internal_neg_p = theta_internal_neg[:, 2] - theta_internal_neg[:, 1];
  vector[nPpts] delta_internal_pos_p = theta_internal_pos[:, 2] - theta_internal_pos[:, 1];
  vector[nPpts] delta_global_neg_p   = theta_global_neg[:, 2] - theta_global_neg[:, 1];
  vector[nPpts] delta_global_pos_p   = theta_global_pos[:, 2] - theta_global_pos[:, 1];

  // log likelihoods
  array[nPpts, nTimes] real sum_log_lik_ppt;
  vector[nPpts] log_lik;

  // generate posterior predictions
  for (p in 1:nPpts) {
    sum_log_lik_ppt[p, 1] = 0;
    sum_log_lik_ppt[p, 2] = 0;
    if (nT_ppts[p, 1] > 0) {
      sum_log_lik_ppt[p, 1] += bernoulli_logit_lpmf(internal_neg[p, 1, 1:nT_ppts[p, 1]] | theta_internal_neg[p, 1]);
      sum_log_lik_ppt[p, 1] += bernoulli_logit_lpmf(internal_pos[p, 1, 1:nT_ppts[p, 1]] | theta_internal_pos[p, 1]);
      sum_log_lik_ppt[p, 1] += bernoulli_logit_lpmf(global_neg[p, 1, 1:nT_ppts[p, 1]] | theta_global_neg[p, 1]);
      sum_log_lik_ppt[p, 1] += bernoulli_logit_lpmf(global_pos[p, 1, 1:nT_ppts[p, 1]] | theta_global_pos[p, 1]);
    }
    if (nT_ppts[p, 2] > 0) {
      sum_log_lik_ppt[p, 2] += bernoulli_logit_lpmf(internal_neg[p, 2, 1:nT_ppts[p, 2]] | theta_internal_neg[p, 2]);
      sum_log_lik_ppt[p, 2] += bernoulli_logit_lpmf(internal_pos[p, 2, 1:nT_ppts[p, 2]] | theta_internal_pos[p, 2]);
      sum_log_lik_ppt[p, 2] += bernoulli_logit_lpmf(global_neg[p, 2, 1:nT_ppts[p, 2]] | theta_global_neg[p, 2]);
      sum_log_lik_ppt[p, 2] += bernoulli_logit_lpmf(global_pos[p, 2, 1:nT_ppts[p, 2]] | theta_global_pos[p, 2]);
    }
    log_lik[p] = sum(sum_log_lik_ppt[p, :]);
  }
}
