################################################################################
###### Test-retest covariance, Δ-reliability and mediator collinearity ##########
################################################################################
#
# Three things, none of which belongs in either task's own modelling script:
#
# 1. TEST-RETEST. Every model here builds its individual-level parameters as
#    diag_pre_multiply(sigma, R_chol) * raw and emits only the CORRELATION matrix, so the covariance
#    is reconstructed in R as diag(sigma) R diag(sigma) -- per draw, so the summary carries the joint
#    uncertainty of the correlation and both sds. Written out so this study's retest reliability can
#    be compared directly with the sibling paper's (causality_training_mh/analyses/outputs/
#    attr_theta_retest_corr.csv, which uses the identical 4x4 parameterisation).
#
# 2. Δ-RELIABILITY. lambda = tau^2/(tau^2 + sigma^2) on each per-participant change score: the share
#    of between-participant variance in a Δ that is signal rather than shrinkage/measurement noise.
#    A mediator with low lambda cannot carry much of an indirect effect however large its a-path is.
#
# 3. MEDIATOR COLLINEARITY. The comparison that justifies this repo's JOINT (mutually-adjusted)
#    mediation models against the sibling repo's one-fit-per-mediator choice. That repo abandoned its
#    joint 4-mediator attribution model precisely because its four Δθ domains ran r ~ 0.8 on posterior
#    means (0.49-0.70 per draw), which flipped b-path signs. Computed PER DRAW, so shrinkage and
#    measurement-error uncertainty propagate -- correlating posterior means understates it.
#    Reported within each task (among that model's own mediators) and ACROSS tasks (attribution Δθ
#    against Δ reward/effort sensitivity), the latter aligned on subID, never on id: prep_data()
#    numbers participants from each task's own long file, so the two id vectors are different people.
#
# A NOTE TO CARRY INTO THE WRITE-UP, not silently fixed here: the attribution model's 4x4 correlation
# carries internal<->global covariance, but the reward-effort model's two 2x2 Choleskys are
# INDEPENDENT, so any rewSens<->effSens covariance is structurally forced to zero at the group level.
# The per-draw Δrew<->Δeff correlation below is therefore an association between two quantities that
# model assumed uncorrelated. Making it a joint 4x4 [rew_t1, eff_t1, rew_t2, eff_t2] would fix that
# but changes the model's structure and its results -- out of scope here.

source("model_fns.R")

pth_t1 <- "STAND-tasks-1/analysis/"
pth_t2 <- "STAND-tasks-2/analysis/"
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 1. data, ids and the participant frames ------------------------------------------
quest_t1 <- read_stand_data(pth_t1, "STAND-tasks-1-self-report-data")
quest_t2 <- read_stand_data(pth_t2, "STAND-tasks-2-self-report-data")
data_long_reff <- read_stand_data(pth_t2, "stand-tasks-both-rew-eff-task-data-long")
data_long_ca <- read_stand_data(pth_t2, "stand-tasks-both-causal-attr-task-data-long")

reff_args <- list(
  df = data_long_reff, type = "effort", excl_eff_task = "standard",
  quest_t1 = quest_t1, quest_t2 = quest_t2, excl_quest = "standard"
)
ca_args <- list(
  df = data_long_ca, type = "attribution",
  quest_t1 = quest_t1, quest_t2 = quest_t2, excl_quest = "standard"
)
stan_ls_reff <- do.call(prep_data, reff_args)
stan_ls_ca <- do.call(prep_data, ca_args)
ids_reff <- do.call(prep_data, c(reff_args, list(ret_ids = TRUE)))
ids_ca <- do.call(prep_data, c(ca_args, list(ret_ids = TRUE)))

# condition coding is INVERTED between the two tasks (prep_data(): effort BA = 1, attribution CR = 1).
# getting this backwards would label every arm-specific correlation with the wrong arm, so the arm
# index vectors are built from each task's own coding and then cross-checked against subID.
arm_ids_reff <- list(
  BA = which(stan_ls_reff$condition == 1L), CR = which(stan_ls_reff$condition == 0L)
)
arm_ids_ca <- list(
  BA = which(stan_ls_ca$condition == 0L), CR = which(stan_ls_ca$condition == 1L)
)
for (a in c("BA", "CR")) {
  subs_reff <- ids_reff$subID[arm_ids_reff[[a]]]
  subs_ca <- ids_ca$subID[arm_ids_ca[[a]]]
  shared <- intersect(subs_reff, subs_ca)
  stopifnot(length(shared) > 0)
  # a participant in arm `a` by the reward-effort coding must be in arm `a` by the attribution coding
  stopifnot(all(shared %in% subs_reff), all(shared %in% subs_ca))
}
cat(sprintf(
  "participants: reward-effort n = %d, attribution n = %d, shared n = %d\n",
  stan_ls_reff$nPpts, stan_ls_ca$nPpts,
  length(intersect(ids_reff$subID, ids_ca$subID))
))

## 2. the correlation-matrix / sd draws --------------------------------------------
# the existing *_draws.rds caches predate these variables, so they are extracted once and cached
# separately. `source_label` and the source file's mtime are recorded in the outputs, because the
# attribution fit is being re-sampled (for the pars_sigma_* priors) and a stale extraction here would
# otherwise be invisible.
reff_itt <- "rewEff-multisess-ITT-intervention-additive"
ca_itt <- "causAttr-multisess-ITT-intervention-additive"

extract_or_cache <- function(cache_path, vars, extract_fn, source_path) {
  if (file.exists(cache_path)) {
    cat(sprintf("  loaded cached draws: %s\n", basename(cache_path)))
    return(readRDS(cache_path))
  }
  if (!all(file.exists(source_path))) {
    warning(
      "cannot extract ", paste(vars, collapse = "/"), ": source not found (",
      paste(basename(source_path), collapse = ", "), ") -- skipping", call. = FALSE
    )
    return(NULL)
  }
  cat(sprintf("  extracting %s from %s ...\n", paste(vars, collapse = ", "),
              paste(basename(source_path), collapse = ", ")))
  d <- extract_fn()
  saveRDS(d, cache_path)
  d
}

reff_src <- paste0(pth_t2, "stan-fits/", reff_itt, "-16000.rds")
reff_cov_draws <- extract_or_cache(
  paste0(pth_t2, "stan-fits/", reff_itt, "_cov_draws.rds"),
  c("R_rewSens", "R_effSens", "sigma_rewSens", "sigma_effSens"),
  function() {
    readRDS(reff_src)$draws(
      format = "df", variables = c("R_rewSens", "R_effSens", "sigma_rewSens", "sigma_effSens")
    )
  },
  reff_src
)

# the attribution ITT fit is stored as raw cmdstan CSVs (~2.1 GB each), not a fit object, so the four
# variables are read selectively rather than by loading the whole thing
ca_src <- paste0(pth_t2, "stan-fits/", ca_itt, "-", 1:4, ".csv")
ca_cov_draws <- extract_or_cache(
  paste0(pth_t2, "stan-fits/", ca_itt, "_cov_draws.rds"),
  c("R_theta_neg", "R_theta_pos", "pars_sigma_neg", "pars_sigma_pos"),
  function() {
    posterior::as_draws_df(cmdstanr::read_cmdstan_csv(
      ca_src, variables = c("R_theta_neg", "R_theta_pos", "pars_sigma_neg", "pars_sigma_pos")
    )$post_warmup_draws)
  },
  ca_src
)

## 3. test-retest correlations and covariances -------------------------------------
# the 4 dimensions of the attribution matrices are ordered [internal_t1, global_t1, internal_t2,
# global_t2] PER VALENCE (see the .stan's transformed parameters), so [1,3] and [2,4] are the
# test-retest entries and [1,2]/[3,4] the within-timepoint internal<->global ones
ca_theta_labs <- c("internal_t1", "global_t1", "internal_t2", "global_t2")
reff_labs <- c("t1", "t2")

retest_corr <- dplyr::bind_rows(
  if (!is.null(reff_cov_draws)) dplyr::bind_rows(
    get_corr_matrix_tbl(reff_cov_draws, "R_rewSens", reff_labs) |>
      dplyr::mutate(task = "reward-effort", parameter = "reward sensitivity"),
    get_corr_matrix_tbl(reff_cov_draws, "R_effSens", reff_labs) |>
      dplyr::mutate(task = "reward-effort", parameter = "effort sensitivity")
  ),
  if (!is.null(ca_cov_draws)) dplyr::bind_rows(
    get_corr_matrix_tbl(ca_cov_draws, "R_theta_neg", ca_theta_labs) |>
      dplyr::mutate(task = "attribution", parameter = "theta (negative events)"),
    get_corr_matrix_tbl(ca_cov_draws, "R_theta_pos", ca_theta_labs) |>
      dplyr::mutate(task = "attribution", parameter = "theta (positive events)")
  )
) |>
  dplyr::mutate(
    pair_type = dplyr::case_when(
      paste(label_i, label_j) %in% c("t1 t2", "internal_t1 internal_t2", "global_t1 global_t2") ~
        "test-retest",
      sub("_t[12]$", "", label_i) != sub("_t[12]$", "", label_j) &
        sub(".*_", "", label_i) == sub(".*_", "", label_j) ~ "within-timepoint",
      TRUE ~ "cross (different parameter and timepoint)"
    ),
    .after = label_j
  ) |>
  dplyr::relocate(task, parameter)

retest_cov <- dplyr::bind_rows(
  if (!is.null(reff_cov_draws)) dplyr::bind_rows(
    get_cov_matrix_tbl(reff_cov_draws, "R_rewSens", "sigma_rewSens", reff_labs) |>
      dplyr::mutate(task = "reward-effort", parameter = "reward sensitivity"),
    get_cov_matrix_tbl(reff_cov_draws, "R_effSens", "sigma_effSens", reff_labs) |>
      dplyr::mutate(task = "reward-effort", parameter = "effort sensitivity")
  ),
  if (!is.null(ca_cov_draws)) dplyr::bind_rows(
    get_cov_matrix_tbl(ca_cov_draws, "R_theta_neg", "pars_sigma_neg", ca_theta_labs) |>
      dplyr::mutate(task = "attribution", parameter = "theta (negative events)"),
    get_cov_matrix_tbl(ca_cov_draws, "R_theta_pos", "pars_sigma_pos", ca_theta_labs) |>
      dplyr::mutate(task = "attribution", parameter = "theta (positive events)")
  )
) |>
  dplyr::relocate(task, parameter)

if (nrow(retest_corr)) readr::write_csv(retest_corr, file.path(out_dir, "task_retest_corr.csv"))
if (nrow(retest_cov)) readr::write_csv(retest_cov, file.path(out_dir, "task_retest_cov.csv"))

## 4. per-participant Δ draws, reliability, and within-task collinearity ------------
rf_i_draws <- readRDS(paste0(pth_t2, "stan-fits/", reff_itt, "_indiv_draws.rds"))
ca_i_draws <- readRDS(paste0(pth_t2, "stan-fits/", ca_itt, "_indiv_draws.rds"))

reff_delta_vars <- c(
  `Δ reward sensitivity` = "delta_rewSens_p", `Δ effort sensitivity` = "delta_effSens_p"
)
ca_delta_vars <- c(
  `Δθ internal-negative` = "delta_internal_neg_p", `Δθ internal-positive` = "delta_internal_pos_p",
  `Δθ global-negative` = "delta_global_neg_p", `Δθ global-positive` = "delta_global_pos_p"
)

# lambda on each Δ: computed from the per-participant posterior mean and sd of the SAME draws the
# mediation models use as mediators
reliability_rows <- function(draws, vars, task) {
  dplyr::bind_rows(lapply(names(vars), function(nm) {
    cols <- grep(paste0("^", vars[[nm]], "\\[\\d+\\]$"), names(draws), value = TRUE)
    stopifnot(length(cols) > 1)
    m <- as.matrix(draws[, cols])
    tibble::tibble(
      task = task, parameter = nm, n_ppts = ncol(m),
      lambda = reliability(colMeans(m), apply(m, 2, stats::sd)),
      sd_between = stats::sd(colMeans(m)), mean_within_sd = mean(apply(m, 2, stats::sd))
    )
  }))
}
delta_reliability <- dplyr::bind_rows(
  reliability_rows(rf_i_draws, reff_delta_vars, "reward-effort"),
  reliability_rows(ca_i_draws, ca_delta_vars, "attribution")
)
readr::write_csv(delta_reliability, file.path(out_dir, "delta_reliability.csv"))

# within-task collinearity among each joint model's own mediators
collin_within <- dplyr::bind_rows(
  mediator_collinearity(rf_i_draws, reff_delta_vars, group_ids = arm_ids_reff,
                        label = "reward-effort (within task)"),
  mediator_collinearity(ca_i_draws, ca_delta_vars, group_ids = arm_ids_ca,
                        label = "attribution (within task)")
)

## 5. cross-task collinearity -------------------------------------------------------
shared_subs <- intersect(ids_reff$subID, ids_ca$subID)
collin_cross <- dplyr::bind_rows(lapply(c("pooled", "BA", "CR"), function(g) {
  subs <- if (g == "pooled") shared_subs else intersect(shared_subs, ids_reff$subID[arm_ids_reff[[g]]])
  cross_task_collinearity(
    ca_i_draws, ca_delta_vars, ids_ca, rf_i_draws, reff_delta_vars, ids_reff,
    sub_ids = subs, group = g
  )
})) |>
  dplyr::mutate(comparison = "attribution vs reward-effort (cross task)", .before = 1)

collinearity <- dplyr::bind_rows(collin_within, collin_cross)
readr::write_csv(collinearity, file.path(out_dir, "mediator_collinearity.csv"))

## 6. reporting ---------------------------------------------------------------------
cat("\n========== TEST-RETEST (individual-level parameters) ==========\n\n")
if (nrow(retest_corr)) {
  rt <- retest_corr |> dplyr::filter(pair_type == "test-retest")
  for (j in seq_len(nrow(rt))) {
    r <- rt[j, ]
    cat(sprintf("  %-14s %-24s %-12s ~ %-12s r = %.3f [%.3f, %.3f]\n",
                r$task, r$parameter, r$label_i, r$label_j, r$mean, r$hdi_95_lo, r$hdi_95_hi))
  }
  cat("\n  within-timepoint (cross-parameter) correlations:\n")
  wt <- retest_corr |> dplyr::filter(pair_type == "within-timepoint")
  for (j in seq_len(nrow(wt))) {
    r <- wt[j, ]
    cat(sprintf("  %-14s %-24s %-12s ~ %-12s r = %.3f [%.3f, %.3f]\n",
                r$task, r$parameter, r$label_i, r$label_j, r$mean, r$hdi_95_lo, r$hdi_95_hi))
  }
} else {
  cat("  (no correlation-matrix draws available -- see the warnings above)\n")
}

cat("\n========== RELIABILITY OF THE CHANGE SCORES (lambda) ==========\n")
cat("lambda = tau^2/(tau^2 + sigma^2): the share of between-participant variance in a Delta that is\n")
cat("signal. a mediator with low lambda cannot carry much of an indirect effect.\n\n")
print(delta_reliability |> dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 3))),
      n = Inf)

cat("\n\n========== MEDIATOR COLLINEARITY ==========\n")
cat("per-draw correlations across participants, so shrinkage/measurement-error uncertainty\n")
cat("propagates. this is the evidence for or against the joint (mutually-adjusted) mediation\n")
cat("models used here: the sibling study moved to one fit per mediator when its four attribution\n")
cat("domains reached r ~ 0.5-0.7 per draw and b-path signs began to flip.\n\n")
print(
  collinearity |>
    dplyr::select(tidyselect::any_of(c("comparison", "group", "par_x", "par_y", "n_ppts", "mean",
                                       "hdi_95_lo", "hdi_95_hi", "pd"))) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 3))),
  n = Inf
)

cat(sprintf(
  "\nmax |r| among the attribution model's four mediators (pooled): %.3f\n",
  max(abs(collin_within$mean[collin_within$comparison == "attribution (within task)" &
                               collin_within$group == "pooled"]))
))
cat(sprintf(
  "max |r| between the reward-effort model's two mediators (pooled): %.3f\n",
  max(abs(collin_within$mean[collin_within$comparison == "reward-effort (within task)" &
                               collin_within$group == "pooled"]))
))

cat("\noutputs written to ", out_dir, "/: task_retest_corr.csv, task_retest_cov.csv, ",
    "delta_reliability.csv, mediator_collinearity.csv\n", sep = "")
