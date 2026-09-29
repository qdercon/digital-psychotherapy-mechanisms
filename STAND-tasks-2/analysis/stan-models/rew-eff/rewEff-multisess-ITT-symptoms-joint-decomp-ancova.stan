// joint moderated mediation, reward/effort sensitivity -> PHQ-9, with an ANCOVA symptom submodel.
//
// the CHOICE submodel is byte-identical to rewEff-multisess-ITT-symptoms-joint-decomp.stan. only the
// symptom half differs, and the reason is a misspecification rather than a refinement:
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
//     contribute fully to the choice submodel, exactly as the primary questionnaire ANCOVA treats them.
//
// PHQ9 is expected t1-centred and divided by sd(t2) -- prep_data(qnr_df = , z_sd_ref = "t2") -- so the
// coefficients are on the same scale as the questionnaire ANCOVA forest and total_arm can be checked
// against it. condition == 1 = BA throughout (OPPOSITE to the causAttr models), so beta_phq_group and
// every _group term is the BA - CR contrast.

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
  array[nPpts, nTimes] real PHQ9;                                      // PHQ-9 score per participant per session (use -999 for missing)
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
  real rewSens_int;
  real effSens_int;

  // group-level effects of non-completion
  real rewSens_ncpl;
  real effSens_ncpl;

  // PHQ9 ANCOVA parameters -- shared
  real alpha_phq;             // intercept
  real beta_t1_phq;           // baseline (t1) PHQ-9 slope, freely estimated
  real beta_phq_group;        // c' path: baseline-adjusted arm contrast (BA - CR), net of both mediators
  real<lower=0> sigma_phq;    // residual sd about the baseline-adjusted mean (no random intercept)

  // PHQ9 regression parameters -- reward-sensitivity mediator
  real beta_phq_rewB;         // between-person (baseline) association
  real beta_phq_rewB_group;   // group difference in between-person association
  real beta_phq_rewW;         // within-person (change) association
  real beta_phq_rewW_group;   // group difference in within-person association

  // PHQ9 regression parameters -- effort-sensitivity mediator
  real beta_phq_effB;         // between-person (baseline) association
  real beta_phq_effB_group;   // group difference in between-person association
  real beta_phq_effW;         // within-person (change) association
  real beta_phq_effW_group;   // group difference in within-person association
}

transformed parameters {
  // individual-level parameter off-sets (for non-centered parameterization)
  matrix[nTimes, nPpts] rewSens_tilde;
  matrix[nTimes, nPpts] effSens_tilde;

  // individual-level parameters
  matrix[nPpts, nTimes] rewSens;
  matrix[nPpts, nTimes] effSens;

  // construct individual offsets (for non-centered parameterization)
  rewSens_tilde = diag_pre_multiply(sigma_rewSens, R_chol_rewSens) * rewSens_raw;
  effSens_tilde = diag_pre_multiply(sigma_effSens, R_chol_effSens) * effSens_raw;

  for (p in 1:nPpts) {
    if (completed[p] == 0) {
      rewSens_tilde[:, p] += rewSens_ncpl;
      effSens_tilde[:, p] += effSens_ncpl;
    }
    // time 1
    rewSens[p, 1] = -1 + Phi_approx(mu_rewSens[1] + rewSens_tilde[1, p]) * 4;
    effSens[p, 1] = -1 + Phi_approx(mu_effSens[1] + effSens_tilde[1, p]) * 10;

    // time 2
    if (condition[p] == 1) {
      rewSens[p, 2] = -1 + Phi_approx(mu_rewSens[2] + rewSens_tilde[2, p] + rewSens_int) * 4;
      effSens[p, 2] = -1 + Phi_approx(mu_effSens[2] + effSens_tilde[2, p] + effSens_int) * 10;
    } else {
      rewSens[p, 2] = -1 + Phi_approx(mu_rewSens[2] + rewSens_tilde[2, p]) * 4;
      effSens[p, 2] = -1 + Phi_approx(mu_effSens[2] + effSens_tilde[2, p]) * 10;
    }
  }

  // within-person deviations from baseline (0 at t=1 by construction) -- both mediators
  matrix[nPpts, nTimes] rewSens_within;
  matrix[nPpts, nTimes] effSens_within;
  for (p in 1:nPpts) {
    for (t in 1:nTimes) {
      rewSens_within[p, t] = rewSens[p, t] - rewSens[p, 1];
      effSens_within[p, t] = effSens[p, t] - effSens[p, 1];
    }
  }

  // grand mean of each mediator's BASELINE level, subtracted in the symptom likelihood below.
  // a pure reparameterisation of the likelihood, but it makes beta_phq_group a comparison at the
  // pooled-mean baseline sensitivity rather than an extrapolation to a sensitivity of 0 (which is
  // outside the range rewSens/effSens actually occupy).
  real rewSens_t1_mean = mean(rewSens[:, 1]);
  real effSens_t1_mean = mean(effSens[:, 1]);
}

model {
  // priors
  R_chol_rewSens ~ lkj_corr_cholesky(1);
  R_chol_effSens ~ lkj_corr_cholesky(1);

  mu_rewSens ~ normal(0, 1);
  mu_effSens ~ normal(0, 1);
  sigma_rewSens ~ cauchy(0,1);
  sigma_effSens ~ cauchy(0,1);

  to_vector(rewSens_raw) ~ normal(0, 1);
  to_vector(effSens_raw) ~ normal(0, 1);

  rewSens_int ~ normal(0,1);
  effSens_int ~ normal(0,1);
  rewSens_ncpl ~ normal(0,1);
  effSens_ncpl ~ normal(0,1);

  // weakly-informative priors on the STANDARDISED PHQ-9 ANCOVA: with PHQ-9 divided by sd(t2),
  // coefficients are on ~unit scale so default N(0,1) is genuinely weakly informative rather than
  // shrinking associations toward 0. beta_t1_phq is a slope between two administrations of the SAME
  // instrument -- dimensionless and expected around 0.4-0.8 -- so it is centred there rather than at
  // 0, matching how the primary questionnaire ANCOVA treats its t1_c term.
  alpha_phq ~ normal(0, 1);
  beta_t1_phq ~ normal(0.6, 0.5);
  beta_phq_group ~ normal(0, 1);
  beta_phq_rewB ~ normal(0, 1);
  beta_phq_rewB_group ~ normal(0, 1);
  beta_phq_rewW ~ normal(0, 1);
  beta_phq_rewW_group ~ normal(0, 1);
  beta_phq_effB ~ normal(0, 1);
  beta_phq_effB_group ~ normal(0, 1);
  beta_phq_effW ~ normal(0, 1);
  beta_phq_effW_group ~ normal(0, 1);
  // sigma_phq is now a residual about a baseline-adjusted mean, so it is SMALLER than the old
  // model's marginal sigma_phq -- roughly sd(t2) * sqrt(1 - r^2) ~ 0.91 on this scale. kept as
  // cauchy(0,1) for consistency with every other sd prior in this repo (the reference
  // implementation uses normal(0, 0.5); at n ~ 96 the likelihood dominates either way, and the
  // heavier tail is the less informative of the two).
  sigma_phq ~ cauchy(0, 1);

  // loop over observations
  for (p in 1:nPpts) {
    // choice model
    if (nT_ppts[p, 1] > 0) {
      vector[nT_ppts[p, 1]] v1_t1 = rewSens[p, 1] * to_vector(rew1[p, 1, 1:nT_ppts[p, 1]]) - effSens[p, 1] * to_vector(eff1[p, 1, 1:nT_ppts[p, 1]]);
      vector[nT_ppts[p, 1]] v2_t1 = rewSens[p, 1] * to_vector(rew2[p, 1, 1:nT_ppts[p, 1]]) - effSens[p, 1] * to_vector(eff2[p, 1, 1:nT_ppts[p, 1]]);
      target += bernoulli_logit_lpmf(choice01[p, 1, 1:nT_ppts[p, 1]] | v2_t1 - v1_t1);
    }

    if (completed[p] == 1 && nT_ppts[p, 2] > 0) {
      vector[nT_ppts[p, 2]] v1_t2 = rewSens[p, 2] * to_vector(rew1[p, 2, 1:nT_ppts[p, 2]]) - effSens[p, 2] * to_vector(eff1[p, 2, 1:nT_ppts[p, 2]]);
      vector[nT_ppts[p, 2]] v2_t2 = rewSens[p, 2] * to_vector(rew2[p, 2, 1:nT_ppts[p, 2]]) - effSens[p, 2] * to_vector(eff2[p, 2, 1:nT_ppts[p, 2]]);
      target += bernoulli_logit_lpmf(choice01[p, 2, 1:nT_ppts[p, 2]] | v2_t2 - v1_t2);
    }

    // PHQ9 ANCOVA -- t2 level conditional on t1, both mediators simultaneously (each adjusted for
    // the other). requires BOTH timepoints observed, so a participant missing either contributes
    // nothing here while still contributing to the choice model above.
    if (PHQ9[p, 1] != -999 && PHQ9[p, 2] != -999) {
      PHQ9[p, 2] ~ normal(
        alpha_phq +
        beta_t1_phq * PHQ9[p, 1] +
        beta_phq_group * condition[p] +
        (beta_phq_rewB + beta_phq_rewB_group * condition[p]) * (rewSens[p, 1] - rewSens_t1_mean) +
        (beta_phq_rewW + beta_phq_rewW_group * condition[p]) * rewSens_within[p, 2] +
        (beta_phq_effB + beta_phq_effB_group * condition[p]) * (effSens[p, 1] - effSens_t1_mean) +
        (beta_phq_effW + beta_phq_effW_group * condition[p]) * effSens_within[p, 2],
        sigma_phq
      );
    }
  }
}

generated quantities {
  corr_matrix[2] R_rewSens = multiply_lower_tri_self_transpose(R_chol_rewSens);
  corr_matrix[2] R_effSens = multiply_lower_tri_self_transpose(R_chol_effSens);
  real delta_rewSens = mu_rewSens[2] - mu_rewSens[1];
  real delta_effSens = mu_effSens[2] - mu_effSens[1];

  // per-participant within-person change in each parameter (for raw-data plots)
  vector[nPpts] delta_rewSens_p;
  vector[nPpts] delta_effSens_p;
  for (p in 1:nPpts) {
    delta_rewSens_p[p] = rewSens[p, 2] - rewSens[p, 1];
    delta_effSens_p[p] = effSens[p, 2] - effSens[p, 1];
  }

  // --- joint moderated mediation: symptom change via BOTH reward- and effort-
  //     sensitivity change in parallel (M1 = within-person delta rewSens, M2 =
  //     within-person delta effSens, Y = standardised PHQ-9 at t2 adjusted for t1).
  //     Each b-path is adjusted for the other mediator (both enter the same PHQ9
  //     regression), unlike the single-mediator rew-decomp/eff-decomp models where
  //     each mediator's b-path omits the other. a-paths are ITT (all participants,
  //     using the model-imputed t2 values for non-completers). condition == 1 = BA.
  real a_rew_cr = 0;
  real a_rew_ba = 0;
  real a_eff_cr = 0;
  real a_eff_ba = 0;
  int n_cr = 0;
  int n_ba = 0;
  // arm means of the BASELINE mediator level, for the total_arm reconstruction below
  real rew_t1_cr = 0;
  real rew_t1_ba = 0;
  real eff_t1_cr = 0;
  real eff_t1_ba = 0;
  for (p in 1:nPpts) {
    real drew = rewSens[p, 2] - rewSens[p, 1];
    real deff = effSens[p, 2] - effSens[p, 1];
    if (condition[p] == 1) {
      a_rew_ba += drew;
      a_eff_ba += deff;
      rew_t1_ba += rewSens[p, 1];
      eff_t1_ba += effSens[p, 1];
      n_ba += 1;
    } else {
      a_rew_cr += drew;
      a_eff_cr += deff;
      rew_t1_cr += rewSens[p, 1];
      eff_t1_cr += effSens[p, 1];
      n_cr += 1;
    }
  }
  a_rew_cr /= n_cr;
  a_rew_ba /= n_ba;
  a_eff_cr /= n_cr;
  a_eff_ba /= n_ba;
  rew_t1_cr /= n_cr;
  rew_t1_ba /= n_ba;
  eff_t1_cr /= n_cr;
  eff_t1_ba /= n_ba;

  // sensitivity companion to the a-paths above: the same arm means restricted to participants who
  // actually have t2 task data, i.e. whose delta is data-driven rather than drawn from the
  // population model. a_*_ba/cr are the ITT quantities (all randomised, missing t2 imputed by the
  // hierarchical structure + *_ncpl); these say how much of that leans on the imputation. only 96
  // of 188 return here, so the gap is worth reporting rather than assuming away. report both, do
  // not silently pick one.
  real a_rew_cr_obs = 0;
  real a_rew_ba_obs = 0;
  real a_eff_cr_obs = 0;
  real a_eff_ba_obs = 0;
  {
    int n_cr_obs = 0;
    int n_ba_obs = 0;
    for (p in 1:nPpts) {
      if (completed[p] == 1 && nT_ppts[p, 2] > 0) {
        real drew = rewSens[p, 2] - rewSens[p, 1];
        real deff = effSens[p, 2] - effSens[p, 1];
        if (condition[p] == 1) {
          a_rew_ba_obs += drew;
          a_eff_ba_obs += deff;
          n_ba_obs += 1;
        } else {
          a_rew_cr_obs += drew;
          a_eff_cr_obs += deff;
          n_cr_obs += 1;
        }
      }
    }
    a_rew_cr_obs /= n_cr_obs;
    a_rew_ba_obs /= n_ba_obs;
    a_eff_cr_obs /= n_cr_obs;
    a_eff_ba_obs /= n_ba_obs;
  }

  // b-paths: within-person change -> PHQ-9 slope, each adjusted for the other mediator
  real b_rew_cr = beta_phq_rewW;
  real b_rew_ba = beta_phq_rewW + beta_phq_rewW_group;
  real b_eff_cr = beta_phq_effW;
  real b_eff_ba = beta_phq_effW + beta_phq_effW_group;

  // direct (c') path: ONE quantity, not the old model's per-arm pair -- an ANCOVA's outcome is the
  // t2 LEVEL, so there is no within-arm change term to report per arm. this is the baseline-adjusted
  // arm contrast (BA - CR) at the pooled-mean baseline sensitivity, net of both mediators' change.
  // NOT in general equal to total_arm - (indirect diffs) when baseline sensitivity is arm-imbalanced
  // -- see total_arm below, which is the exact reconstruction.
  real direct_arm = beta_phq_group;

  // indirect effects (a * b) per mediator, and each mediator's moderated-mediation index (BA - CR)
  real indirect_rew_cr = a_rew_cr * b_rew_cr;
  real indirect_rew_ba = a_rew_ba * b_rew_ba;
  real indirect_eff_cr = a_eff_cr * b_eff_cr;
  real indirect_eff_ba = a_eff_ba * b_eff_ba;
  real index_mod_med_rew = indirect_rew_ba - indirect_rew_cr;
  real index_mod_med_eff = indirect_eff_ba - indirect_eff_cr;

  // total (marginal) arm effect on the baseline-adjusted t2 outcome, reconstructed EXACTLY.
  // direct_arm + the indirect diffs is not it on its own: the baseline mediator terms are centred at
  // the POOLED mean, so averaging the model's own prediction over each arm's actual (possibly
  // imbalanced) baseline sensitivity leaves the two between-person terms below, which vanish only if
  // baseline sensitivity happens to be balanced across arms. the within-person part of the arm
  // difference is algebraically (a_ba * b_ba - a_cr * b_cr) = index_mod_med_*, hence those terms.
  real total_arm = direct_arm
    + beta_phq_rewB * (rew_t1_ba - rew_t1_cr)
    + beta_phq_rewB_group * (rew_t1_ba - rewSens_t1_mean)
    + beta_phq_effB * (eff_t1_ba - eff_t1_cr)
    + beta_phq_effB_group * (eff_t1_ba - effSens_t1_mean)
    + index_mod_med_rew
    + index_mod_med_eff;

  // deliberately NO prop_med_*: those were each mediator's share of a per-arm total effect, and an
  // ANCOVA has no per-arm total. they were also unstable whenever that total was near 0.

  // log-likelihoods -- CHOICE MODEL ONLY, as in the original model, so LOO comparisons are on the
  // task submodel and are not silently changed by the symptom reformulation.
  array[nPpts, nTimes] real sum_log_lik_ppt;
  vector[nPpts] log_lik;

  for (p in 1:nPpts) {
    sum_log_lik_ppt[p, 1] = 0;
    sum_log_lik_ppt[p, 2] = 0;

    if (nT_ppts[p, 1] > 0) {
      vector[nT_ppts[p, 1]] v1_t1 = rewSens[p, 1] * to_vector(rew1[p, 1, 1:nT_ppts[p, 1]]) - effSens[p, 1] * to_vector(eff1[p, 1, 1:nT_ppts[p, 1]]);
      vector[nT_ppts[p, 1]] v2_t1 = rewSens[p, 1] * to_vector(rew2[p, 1, 1:nT_ppts[p, 1]]) - effSens[p, 1] * to_vector(eff2[p, 1, 1:nT_ppts[p, 1]]);
      sum_log_lik_ppt[p, 1] = bernoulli_logit_lpmf(choice01[p, 1, 1:nT_ppts[p, 1]] | v2_t1 - v1_t1);
    }

    if (completed[p] == 1 && nT_ppts[p, 2] > 0) {
      vector[nT_ppts[p, 2]] v1_t2 = rewSens[p, 2] * to_vector(rew1[p, 2, 1:nT_ppts[p, 2]]) - effSens[p, 2] * to_vector(eff1[p, 2, 1:nT_ppts[p, 2]]);
      vector[nT_ppts[p, 2]] v2_t2 = rewSens[p, 2] * to_vector(rew2[p, 2, 1:nT_ppts[p, 2]]) - effSens[p, 2] * to_vector(eff2[p, 2, 1:nT_ppts[p, 2]]);
      sum_log_lik_ppt[p, 2] = bernoulli_logit_lpmf(choice01[p, 2, 1:nT_ppts[p, 2]] | v2_t2 - v1_t2);
    }
    log_lik[p] = sum(sum_log_lik_ppt[p, :]);
  }
}
