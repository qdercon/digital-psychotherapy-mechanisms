// joint moderated mediation, causal-attribution tendencies -> PHQ-9, with an ANCOVA symptom submodel.
//
// the ATTRIBUTION submodel is byte-identical to causAttr-multisess-ITT-symptoms-joint-decomp.stan
// (including the newly-added pars_sigma_* priors). only the symptom half differs, and the reason is a
// misspecification rather than a refinement:
//
// the old submodel put a participant random intercept (ranef_p ~ normal(0, sigma_p)) on BOTH PHQ-9
// rows with one i.i.d. residual sigma_phq. that is compound symmetry -- it forces Var(t1) == Var(t2).
// measured on this data (96 complete cases) PHQ-9 goes sd 3.29 -> 3.86, a variance ratio of 1.37,
// because enrolment screened on baseline PHQ-9 but not on follow-up. test-retest r is only 0.41.
//
// this version instead conditions on the observed baseline: the outcome is the t2 LEVEL and t1 enters
// as a freely-estimated covariate, where the random intercept forced that slope to whatever the
// compound-symmetry variance ratio implied. one row per participant, so there is no within-person
// covariance left to model -- that is the point of the parameterisation.
//
// consequences for what can be reported, all of them structural:
//   * there is ONE direct effect (direct_arm = beta_phq_group), not a per-arm pair. an ANCOVA has no
//     within-arm change term, so direct_ba / direct_cr are not defined here.
//   * likewise one total_arm, not total_ba / total_cr, and no prop_med_* (each mediator's share of a
//     per-arm total that no longer exists, and unstable whenever that total was near 0).
//   * a participant missing EITHER timepoint contributes nothing to the symptom submodel. they still
//     contribute fully to the attribution submodel, exactly as the primary questionnaire ANCOVA does.
//
// PHQ9 is expected t1-centred and divided by sd(t2) -- prep_data(qnr_df = , z_sd_ref = "t2") -- so the
// coefficients are on the same scale as the questionnaire ANCOVA forest and total_arm can be checked
// against it. condition == 1 = CR here (OPPOSITE to the rewEff models), so beta_phq_group and every
// _group term is the CR - BA contrast, and total_arm must be NEGATED before comparing with a
// questionnaire forest drawn as BA - CR.

data {
  int nTimes;                                                              // number of time points (sessions)
  int nPpts;                                                               // number of participants
  array[nPpts] int condition;                                              // intervention condition for each participant (0 = BA, 1 = CR)
  array[nPpts] int completed;                                              // whether participant completed time 2
  int nTrials_max;                                                         // maximum number of trials per participant per session
  array[nPpts, nTimes] int nT_ppts;                                        // actual number of trials per participant per session
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> internal_neg;   // responses for each participant per session
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> internal_pos;   // responses for each participant per session
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> global_neg;     // responses for each participant per session
  array[nPpts, nTimes, nTrials_max] int<lower=-1, upper=1> global_pos;     // responses for each participant per session
  array[nPpts, nTimes] real PHQ9;                                          // PHQ-9 z-score per participant per session (use -999 for missing)
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

  // PHQ9 ANCOVA parameters -- shared
  real alpha_phq;             // intercept
  real beta_t1_phq;           // baseline (t1) PHQ-9 slope, freely estimated
  real beta_phq_group;        // c' path: baseline-adjusted arm contrast (CR - BA), net of all 4 mediators
  real<lower=0> sigma_phq;    // residual sd about the baseline-adjusted mean (no random intercept)

  // PHQ9 regression parameters -- internal-negative mediator
  real beta_phq_int_negB;         // between-person (baseline) association
  real beta_phq_int_negB_group;   // group difference in between-person association
  real beta_phq_int_negW;         // within-person (change) association
  real beta_phq_int_negW_group;   // group difference in within-person association

  // PHQ9 regression parameters -- internal-positive mediator
  real beta_phq_int_posB;
  real beta_phq_int_posB_group;
  real beta_phq_int_posW;
  real beta_phq_int_posW_group;

  // PHQ9 regression parameters -- global-negative mediator
  real beta_phq_glob_negB;
  real beta_phq_glob_negB_group;
  real beta_phq_glob_negW;
  real beta_phq_glob_negW_group;

  // PHQ9 regression parameters -- global-positive mediator
  real beta_phq_glob_posB;
  real beta_phq_glob_posB_group;
  real beta_phq_glob_posW;
  real beta_phq_glob_posW_group;
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

  // within-person deviations from baseline (0 at t=1 by construction) -- all four mediators
  matrix[nPpts, nTimes] theta_internal_neg_within;
  matrix[nPpts, nTimes] theta_internal_pos_within;
  matrix[nPpts, nTimes] theta_global_neg_within;
  matrix[nPpts, nTimes] theta_global_pos_within;
  for (p in 1:nPpts) {
    for (t in 1:nTimes) {
      theta_internal_neg_within[p, t] = theta_internal_neg[p, t] - theta_internal_neg[p, 1];
      theta_internal_pos_within[p, t] = theta_internal_pos[p, t] - theta_internal_pos[p, 1];
      theta_global_neg_within[p, t]   = theta_global_neg[p, t]   - theta_global_neg[p, 1];
      theta_global_pos_within[p, t]   = theta_global_pos[p, t]   - theta_global_pos[p, 1];
    }
  }

  // grand mean of each mediator's BASELINE level, subtracted in the symptom likelihood below.
  // a pure reparameterisation of the likelihood, but it makes beta_phq_group a comparison at the
  // pooled-mean baseline theta rather than an extrapolation to theta = 0 on the logit scale (which
  // is a 50% attribution rate, well outside where internal-positive actually sits).
  real theta_internal_neg_t1_mean = mean(theta_internal_neg[:, 1]);
  real theta_internal_pos_t1_mean = mean(theta_internal_pos[:, 1]);
  real theta_global_neg_t1_mean   = mean(theta_global_neg[:, 1]);
  real theta_global_pos_t1_mean   = mean(theta_global_pos[:, 1]);
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

  // weakly-informative priors on the STANDARDISED PHQ-9 ANCOVA: with PHQ-9 divided by sd(t2),
  // coefficients are on ~unit scale so default N(0,1) is genuinely weakly informative rather than
  // shrinking associations toward 0. beta_t1_phq is a slope between two administrations of the SAME
  // instrument -- dimensionless and expected around 0.4-0.8 -- so it is centred there rather than at
  // 0, matching how the primary questionnaire ANCOVA treats its t1_c term.
  alpha_phq ~ normal(0, 1);
  beta_t1_phq ~ normal(0.6, 0.5);
  beta_phq_group ~ normal(0, 1);
  beta_phq_int_negB ~ normal(0, 1);
  beta_phq_int_negB_group ~ normal(0, 1);
  beta_phq_int_negW ~ normal(0, 1);
  beta_phq_int_negW_group ~ normal(0, 1);
  beta_phq_int_posB ~ normal(0, 1);
  beta_phq_int_posB_group ~ normal(0, 1);
  beta_phq_int_posW ~ normal(0, 1);
  beta_phq_int_posW_group ~ normal(0, 1);
  beta_phq_glob_negB ~ normal(0, 1);
  beta_phq_glob_negB_group ~ normal(0, 1);
  beta_phq_glob_negW ~ normal(0, 1);
  beta_phq_glob_negW_group ~ normal(0, 1);
  beta_phq_glob_posB ~ normal(0, 1);
  beta_phq_glob_posB_group ~ normal(0, 1);
  beta_phq_glob_posW ~ normal(0, 1);
  beta_phq_glob_posW_group ~ normal(0, 1);
  // sigma_phq is now a residual about a baseline-adjusted mean, so it is SMALLER than the old
  // model's marginal sigma_phq -- roughly sd(t2) * sqrt(1 - r^2) ~ 0.91 on this scale. kept as
  // cauchy(0,1) for consistency with every other sd prior in this repo (the reference
  // implementation uses normal(0, 0.5); at n ~ 96 the likelihood dominates either way, and the
  // heavier tail is the less informative of the two).
  sigma_phq ~ cauchy(0, 1);

  // loop over observations
  for (p in 1:nPpts) {
    // choice model (identical to causAttr-multisess-ITT-intervention-additive.stan)
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

    // PHQ9 ANCOVA -- t2 level conditional on t1, all four mediators simultaneously (each adjusted
    // for the other three). requires BOTH timepoints observed, so a participant missing either
    // contributes nothing here while still contributing to the attribution model above.
    if (PHQ9[p, 1] != -999 && PHQ9[p, 2] != -999) {
      PHQ9[p, 2] ~ normal(
        alpha_phq +
        beta_t1_phq * PHQ9[p, 1] +
        beta_phq_group * condition[p] +
        (beta_phq_int_negB + beta_phq_int_negB_group * condition[p])
          * (theta_internal_neg[p, 1] - theta_internal_neg_t1_mean) +
        (beta_phq_int_negW + beta_phq_int_negW_group * condition[p]) * theta_internal_neg_within[p, 2] +
        (beta_phq_int_posB + beta_phq_int_posB_group * condition[p])
          * (theta_internal_pos[p, 1] - theta_internal_pos_t1_mean) +
        (beta_phq_int_posW + beta_phq_int_posW_group * condition[p]) * theta_internal_pos_within[p, 2] +
        (beta_phq_glob_negB + beta_phq_glob_negB_group * condition[p])
          * (theta_global_neg[p, 1] - theta_global_neg_t1_mean) +
        (beta_phq_glob_negW + beta_phq_glob_negW_group * condition[p]) * theta_global_neg_within[p, 2] +
        (beta_phq_glob_posB + beta_phq_glob_posB_group * condition[p])
          * (theta_global_pos[p, 1] - theta_global_pos_t1_mean) +
        (beta_phq_glob_posW + beta_phq_glob_posW_group * condition[p]) * theta_global_pos_within[p, 2],
        sigma_phq
      );
    }
  }
}

generated quantities {
  // test-retest correlations
  corr_matrix[4] R_theta_neg = multiply_lower_tri_self_transpose(R_chol_theta_neg);
  corr_matrix[4] R_theta_pos = multiply_lower_tri_self_transpose(R_chol_theta_pos);

  // differences in group-level parameters
  real delta_internal_neg = mu_internal_theta_neg[2] - mu_internal_theta_neg[1];
  real delta_internal_pos = mu_internal_theta_pos[2] - mu_internal_theta_pos[1];
  real delta_global_neg   = mu_global_theta_neg[2] - mu_global_theta_neg[1];
  real delta_global_pos   = mu_global_theta_pos[2] - mu_global_theta_pos[1];

  // per-participant within-person change in each parameter (for raw-data plots)
  vector[nPpts] delta_internal_neg_p = theta_internal_neg[:, 2] - theta_internal_neg[:, 1];
  vector[nPpts] delta_internal_pos_p = theta_internal_pos[:, 2] - theta_internal_pos[:, 1];
  vector[nPpts] delta_global_neg_p   = theta_global_neg[:, 2] - theta_global_neg[:, 1];
  vector[nPpts] delta_global_pos_p   = theta_global_pos[:, 2] - theta_global_pos[:, 1];

  // --- joint moderated mediation: symptom change via all four attribution-style
  //     parameters in parallel (X = time, M1..M4 = within-person delta theta_*,
  //     Y = standardised PHQ-9). Each b-path is adjusted for the other three
  //     mediators (all four enter the same PHQ9 regression), unlike separate
  //     single-mediator models where each mediator's b-path would omit the
  //     other three. a-paths are ITT (all participants, using the model-imputed
  //     t2 theta for non-completers).
  //     NB: for causal-attribution data, condition == 1 = CR, condition == 0 = BA
  //     -- the OPPOSITE convention to the reward-effort models (see prep_data()
  //     and modelling_causal-attr.R) -- so here the base beta_phq_*W term is BA
  //     and the +_group term adds CR, and a/indirect "cr"/"ba" labels below
  //     follow that same condition==1-is-CR mapping.
  real a_int_neg_cr = 0;
  real a_int_neg_ba = 0;
  real a_int_pos_cr = 0;
  real a_int_pos_ba = 0;
  real a_glob_neg_cr = 0;
  real a_glob_neg_ba = 0;
  real a_glob_pos_cr = 0;
  real a_glob_pos_ba = 0;
  int n_cr = 0;
  int n_ba = 0;
  for (p in 1:nPpts) {
    real d_in = theta_internal_neg[p, 2] - theta_internal_neg[p, 1];
    real d_ip = theta_internal_pos[p, 2] - theta_internal_pos[p, 1];
    real d_gn = theta_global_neg[p, 2] - theta_global_neg[p, 1];
    real d_gp = theta_global_pos[p, 2] - theta_global_pos[p, 1];
    if (condition[p] == 1) {
      a_int_neg_cr += d_in;
      a_int_pos_cr += d_ip;
      a_glob_neg_cr += d_gn;
      a_glob_pos_cr += d_gp;
      n_cr += 1;
    } else {
      a_int_neg_ba += d_in;
      a_int_pos_ba += d_ip;
      a_glob_neg_ba += d_gn;
      a_glob_pos_ba += d_gp;
      n_ba += 1;
    }
  }
  a_int_neg_cr /= n_cr;
  a_int_neg_ba /= n_ba;
  a_int_pos_cr /= n_cr;
  a_int_pos_ba /= n_ba;
  a_glob_neg_cr /= n_cr;
  a_glob_neg_ba /= n_ba;
  a_glob_pos_cr /= n_cr;
  a_glob_pos_ba /= n_ba;

  // arm means of the BASELINE theta level, for the total_arm reconstruction below
  real theta_int_neg_t1_cr = 0;
  real theta_int_neg_t1_ba = 0;
  real theta_int_pos_t1_cr = 0;
  real theta_int_pos_t1_ba = 0;
  real theta_glob_neg_t1_cr = 0;
  real theta_glob_neg_t1_ba = 0;
  real theta_glob_pos_t1_cr = 0;
  real theta_glob_pos_t1_ba = 0;
  for (p in 1:nPpts) {
    if (condition[p] == 1) {
      theta_int_neg_t1_cr  += theta_internal_neg[p, 1];
      theta_int_pos_t1_cr  += theta_internal_pos[p, 1];
      theta_glob_neg_t1_cr += theta_global_neg[p, 1];
      theta_glob_pos_t1_cr += theta_global_pos[p, 1];
    } else {
      theta_int_neg_t1_ba  += theta_internal_neg[p, 1];
      theta_int_pos_t1_ba  += theta_internal_pos[p, 1];
      theta_glob_neg_t1_ba += theta_global_neg[p, 1];
      theta_glob_pos_t1_ba += theta_global_pos[p, 1];
    }
  }
  theta_int_neg_t1_cr  /= n_cr;
  theta_int_pos_t1_cr  /= n_cr;
  theta_glob_neg_t1_cr /= n_cr;
  theta_glob_pos_t1_cr /= n_cr;
  theta_int_neg_t1_ba  /= n_ba;
  theta_int_pos_t1_ba  /= n_ba;
  theta_glob_neg_t1_ba /= n_ba;
  theta_glob_pos_t1_ba /= n_ba;

  // sensitivity companion to the a-paths above: the same arm means restricted to participants who
  // actually have t2 attribution trials, i.e. whose delta theta is data-driven rather than drawn
  // from the population model. a_*_cr/ba are the ITT quantities (all randomised, missing t2 imputed
  // by the hierarchical structure + theta_ncompl_*); these say how much of that leans on the
  // imputation. only 96 of 188 return here, so the gap is worth reporting rather than assuming
  // away. report both, do not silently pick one.
  real a_int_neg_cr_obs = 0;
  real a_int_neg_ba_obs = 0;
  real a_int_pos_cr_obs = 0;
  real a_int_pos_ba_obs = 0;
  real a_glob_neg_cr_obs = 0;
  real a_glob_neg_ba_obs = 0;
  real a_glob_pos_cr_obs = 0;
  real a_glob_pos_ba_obs = 0;
  {
    int n_cr_obs = 0;
    int n_ba_obs = 0;
    for (p in 1:nPpts) {
      if (nT_ppts[p, 2] > 0) {
        real d_in = theta_internal_neg[p, 2] - theta_internal_neg[p, 1];
        real d_ip = theta_internal_pos[p, 2] - theta_internal_pos[p, 1];
        real d_gn = theta_global_neg[p, 2] - theta_global_neg[p, 1];
        real d_gp = theta_global_pos[p, 2] - theta_global_pos[p, 1];
        if (condition[p] == 1) {
          a_int_neg_cr_obs += d_in;
          a_int_pos_cr_obs += d_ip;
          a_glob_neg_cr_obs += d_gn;
          a_glob_pos_cr_obs += d_gp;
          n_cr_obs += 1;
        } else {
          a_int_neg_ba_obs += d_in;
          a_int_pos_ba_obs += d_ip;
          a_glob_neg_ba_obs += d_gn;
          a_glob_pos_ba_obs += d_gp;
          n_ba_obs += 1;
        }
      }
    }
    a_int_neg_cr_obs /= n_cr_obs;
    a_int_pos_cr_obs /= n_cr_obs;
    a_glob_neg_cr_obs /= n_cr_obs;
    a_glob_pos_cr_obs /= n_cr_obs;
    a_int_neg_ba_obs /= n_ba_obs;
    a_int_pos_ba_obs /= n_ba_obs;
    a_glob_neg_ba_obs /= n_ba_obs;
    a_glob_pos_ba_obs /= n_ba_obs;
  }

  // b-paths: within-person change -> PHQ-9 slope, each adjusted for the other three
  // mediators (base term = BA [condition==0], +_group term adds CR [condition==1])
  real b_int_neg_ba = beta_phq_int_negW;
  real b_int_neg_cr = beta_phq_int_negW + beta_phq_int_negW_group;
  real b_int_pos_ba = beta_phq_int_posW;
  real b_int_pos_cr = beta_phq_int_posW + beta_phq_int_posW_group;
  real b_glob_neg_ba = beta_phq_glob_negW;
  real b_glob_neg_cr = beta_phq_glob_negW + beta_phq_glob_negW_group;
  real b_glob_pos_ba = beta_phq_glob_posW;
  real b_glob_pos_cr = beta_phq_glob_posW + beta_phq_glob_posW_group;

  // direct (c') path: ONE quantity, not the old model's per-arm pair -- an ANCOVA's outcome is the
  // t2 LEVEL, so there is no within-arm change term to report per arm. this is the baseline-adjusted
  // arm contrast (CR - BA) at the pooled-mean baseline theta, net of all four mediators' change.
  // NOT in general equal to total_arm - (indirect diffs) when baseline theta is arm-imbalanced --
  // see total_arm below, which is the exact reconstruction.
  real direct_arm = beta_phq_group;

  // indirect effects (a * b) per mediator, and each mediator's moderated-mediation index (CR - BA)
  real indirect_int_neg_cr = a_int_neg_cr * b_int_neg_cr;
  real indirect_int_neg_ba = a_int_neg_ba * b_int_neg_ba;
  real indirect_int_pos_cr = a_int_pos_cr * b_int_pos_cr;
  real indirect_int_pos_ba = a_int_pos_ba * b_int_pos_ba;
  real indirect_glob_neg_cr = a_glob_neg_cr * b_glob_neg_cr;
  real indirect_glob_neg_ba = a_glob_neg_ba * b_glob_neg_ba;
  real indirect_glob_pos_cr = a_glob_pos_cr * b_glob_pos_cr;
  real indirect_glob_pos_ba = a_glob_pos_ba * b_glob_pos_ba;

  real index_mod_med_int_neg  = indirect_int_neg_cr  - indirect_int_neg_ba;
  real index_mod_med_int_pos  = indirect_int_pos_cr  - indirect_int_pos_ba;
  real index_mod_med_glob_neg = indirect_glob_neg_cr - indirect_glob_neg_ba;
  real index_mod_med_glob_pos = indirect_glob_pos_cr - indirect_glob_pos_ba;

  // total (marginal) arm effect on the baseline-adjusted t2 outcome, reconstructed EXACTLY.
  // direct_arm + the indirect diffs is not it on its own: the baseline theta terms are centred at
  // the POOLED mean, so averaging the model's own prediction over each arm's actual (possibly
  // imbalanced) baseline theta leaves the four between-person terms below, which vanish only if
  // baseline theta happens to be balanced across arms. it is not, for this design. the within-person
  // part of the arm difference is algebraically (a_cr * b_cr - a_ba * b_ba) = index_mod_med_*,
  // hence those terms.
  real total_arm = direct_arm
    + beta_phq_int_negB * (theta_int_neg_t1_cr - theta_int_neg_t1_ba)
    + beta_phq_int_negB_group * (theta_int_neg_t1_cr - theta_internal_neg_t1_mean)
    + beta_phq_int_posB * (theta_int_pos_t1_cr - theta_int_pos_t1_ba)
    + beta_phq_int_posB_group * (theta_int_pos_t1_cr - theta_internal_pos_t1_mean)
    + beta_phq_glob_negB * (theta_glob_neg_t1_cr - theta_glob_neg_t1_ba)
    + beta_phq_glob_negB_group * (theta_glob_neg_t1_cr - theta_global_neg_t1_mean)
    + beta_phq_glob_posB * (theta_glob_pos_t1_cr - theta_glob_pos_t1_ba)
    + beta_phq_glob_posB_group * (theta_glob_pos_t1_cr - theta_global_pos_t1_mean)
    + index_mod_med_int_neg
    + index_mod_med_int_pos
    + index_mod_med_glob_neg
    + index_mod_med_glob_pos;

  // deliberately NO prop_med_*: those were each mediator's share of a per-arm total effect, and an
  // ANCOVA has no per-arm total. they were also unstable whenever that total was near 0.

  // log-likelihoods -- ATTRIBUTION MODEL ONLY, as in the original model, so LOO comparisons
  // are on the task submodel and are not silently changed by the symptom reformulation.
  array[nPpts, nTimes] real sum_log_lik_ppt;
  vector[nPpts] log_lik;

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
