## MODELLING -------------------------------------------------------------------

library(patchwork)
source("model_fns.R")

pth_t1 <- "STAND-tasks-1/analysis/"
pth_t2 <- "STAND-tasks-2/analysis/"
dir.create("outputs", showWarnings = FALSE)
dir.create(paste0(pth_t2, "stan-fits"), showWarnings = FALSE)

## Load questionnaire data for exclusions
quest_t1 <- read_stand_data(pth_t1, "STAND-tasks-1-self-report-data")
quest_t2 <- read_stand_data(pth_t2, "STAND-tasks-2-self-report-data")
qqnrs_t2 <- read_stand_data(pth_t2, "qqnrs_t2")

data_long_reff <- read_stand_data(pth_t2, "stand-tasks-both-rew-eff-task-data-long")

stan_ls_reff <- prep_data(
  data_long_reff, type = "effort", excl_eff_task = "standard",
  quest_t1 = quest_t1, quest_t2 = quest_t2, excl_quest = "standard"
)

# MAIN TEXT =====================================================================
## model-free analyses ---------------------------------------------------------

mf_df <- prep_data(
  data_long_reff, type = "effort", ret_df = TRUE, excl_eff_task = "standard",
  quest_t1 = quest_t1, quest_t2 = quest_t2, excl_quest = "standard"
) |>
  dplyr::mutate(
    effPropDiff = abs(trialEffortPropMax1 - trialEffortPropMax2),
    higherEffRew = higherEffChosen * higherRewChosen
  )

rew_plot <- reff_plot(
  mf_df,
  plot_type = "reward",
  mode = "randomised",
  label_colour = "#b75347",
  ln_width = 1.2,
  err_width = 0.3,
  dodge_width = 0.45,
  point_size = 3.5,
  legend_pos = "bottom",
  legend_fnt_sc = 13,
  fnt_sz = 1.6,
  brew_col = "Archambault",
  col_nums = c(1, 4, 3)
)

eff_plot <- reff_plot(
  mf_df,
  plot_type = "effort",
  mode = "randomised",
  label_colour = "#2f70a1",
  ln_width = 1.2,
  err_width = 4,
  dodge_width = 6,
  point_size = 3.5,
  legend_pos = "bottom",
  legend_fnt_sc = 13,
  fnt_sz = 1.6,
  brew_col = "Archambault",
  col_nums = c(1, 4, 3)
)

mf_plots <- wrap_elements(
  rew_plot + eff_plot +
    plot_layout(nrow = 1, guides = "collect") &
    plot_annotation(tag_levels = list(c("A", "B"))) &
    ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
    )
)

# model-free analyses of choice of higher effort by group, timepoint, and symptoms

mf_symp_df <- mf_df |>
  dplyr::left_join(
    qqnrs_t2 |> dplyr::select(subID, PHQ9_total_t1, PHQ9_total_t2), by = c("subID")
  ) |>
  dplyr::mutate(
    intCond = factor(intCond, levels = c("CR", "BA")),
    PHQ9_total = dplyr::if_else(timepoint == 0, PHQ9_total_t1, PHQ9_total_t2)
  ) |>
  dplyr::select(subID, timepoint, intCond, trialNo, higherEffRew, non_completed, PHQ9_total)

higherEffRew_mod <- brms::brm(
  higherEffRew ~ intCond * timepoint + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = mf_symp_df,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(higherEffRew_mod, rel = "BA")

## baseline in non-completers?
bsl_data <- mf_symp_df |> dplyr::filter(timepoint == 0)

higherEffRew_ncpl <- brms::brm(
  higherEffRew ~ non_completed + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = bsl_data,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(higherEffRew_ncpl, rel = NULL, var = "non_completed")

## model-based analyses ---------------------------------------------------------

rf_itt <- "rewEff-multisess-ITT-intervention-additive"
mod_rf_itt <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_itt, ".stan"))

## fit models with cmdstan
fit_rf_itt <- mod_rf_itt$sample(
  data = stan_ls_reff,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
fit_rf_itt$save_object(file = paste0(pth_t2, "stan-fits/", rf_itt, "-16000.rds"))
# fit_rf_itt <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt, "-16000.rds"))

# check we have ~10000 effective samples for the parameter of interest
fit_rf_itt$summary(variables = "effSens_int")
fit_rf_loo <- fit_rf_itt$loo(variables = "log_lik", cores = 4)
# saveRDS(fit_rf_loo, file = paste0(pth_t2, "stan-fits/", rf_itt, "-loo.rds"))

## Extract draws
rf_vars <- c(
  "mu_effSens[1]", "mu_effSens[2]", "mu_rewSens[1]", "mu_rewSens[2]",
  "effSens_int", "rewSens_int", "effSens_ncpl", "rewSens_ncpl",
  "delta_effSens", "delta_rewSens"
)
draws_rf_itt <- fit_rf_itt$draws(format = "df", variables = rf_vars)
saveRDS(draws_rf_itt, paste0(pth_t2, "stan-fits/", rf_itt, "_draws.rds"))
# draws_rf_itt <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt, "_draws.rds"))

long_draws_rf <- draws_rf_itt |> tidybayes::gather_draws(!!!rlang::syms(rf_vars))
tidy_rf_itt <- extract_draws(draws_rf_itt, long_draws_rf, hdi = 0.95)
print(tidy_rf_itt$hdi)
write.csv(tidy_rf_itt$hdi, file = paste0("outputs/", rf_itt, "_hdi.csv"))

# probability of direction
pd_es <- tidy_rf_itt$df |> dplyr::filter(.variable == "effSens_int") |> dplyr::pull(.value)
pd_rs_cr <- tidy_rf_itt$df |> dplyr::filter(.variable == "delta_rewSens") |> dplyr::pull(.value)
pd_rs_ba <- tidy_rf_itt$df |> dplyr::filter(.variable == "delta_rewSens_int") |> dplyr::pull(.value)
pd_ncpl_es <- tidy_rf_itt$df |> dplyr::filter(.variable == "effSens_ncpl") |> dplyr::pull(.value)
pd_ncpl_rs <- tidy_rf_itt$df |> dplyr::filter(.variable == "rewSens_ncpl") |> dplyr::pull(.value)

bayestestR::p_direction(pd_es, method = "direct")
bayestestR::p_direction(pd_rs_cr, method = "direct")
bayestestR::p_direction(pd_rs_ba, method = "direct")

bayestestR::p_direction(pd_ncpl_es, method = "direct")
bayestestR::p_direction(pd_ncpl_rs, method = "direct")

## Individuals' changes in parameters, by group

rf_i_pars <- c("rewSens", "effSens", "delta_effSens_p", "delta_rewSens_p")
rf_i_draws <- fit_rf_itt$draws(format = "df", variables = rf_i_pars)
saveRDS(rf_i_draws, paste0(pth_t2, "stan-fits/", rf_itt, "_indiv_draws.rds"))
# rf_i_draws <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt, "_indiv_draws.rds"))
rf_hdi <- extract_indiv_pars(rf_i_draws, stan_ls_reff)

## combined group-level results: pre->post change per arm + their difference
rew_sens_plts <- plot_reff_change(
  draws_rf_itt,
  diff_lab = "behav. activ.\nvs. cogn. restr.",
  flip_order = TRUE,
  param = "rewSens",
  param_nm = "reward sensitivity",
  col_nums = c(1, 4),
  fnt_sz = 1.8
)
eff_sens_plts <- plot_reff_change(
  draws_rf_itt,
  diff_lab = "behav. activ.\nvs. cogn. restr.",
  flip_order = TRUE,
  suppress_labels = TRUE,
  param = "effSens",
  param_nm = "effort sensitivity",
  col_nums = c(1, 4),
  fnt_sz = 1.8
)

rew_sens_change <- wrap_elements(
  rew_sens_plts$change_panel &
    plot_annotation(tag_levels = list("C")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
    )
)
eff_sens_change <- wrap_elements(
  eff_sens_plts$change_panel &
    plot_annotation(tag_levels = list("D")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = "Open Sans", size = 28, hjust = 1.5
      )
    )
)

figure_3_abcd <-
  wrap_elements(
    mf_plots / (
      rew_sens_change + eff_sens_change +
        plot_layout(ncol = 2, widths = c(0.56, 0.44)) &
        ggplot2::theme(
          plot.margin = ggplot2::margin(t = 0, b = 0, l = 0, r = 0)
        )
    ) + plot_layout(nrow = 2)
  )

## regression of change in parameters by group
# summary(
#   lm(
#     mean ~ group, weights = hdi_width,
#     data = rf_hdi[rf_hdi$variable == "delta_effSens" & rf_hdi$timepoint == "t2-t1", ]
#   )
# )

# reward sensitivity
# mod_rewSens_indiv <- lme4::lmer(
#   mean ~ group * t + (1 | id),
#   weights = 1 / hdi_width,
#   data = rf_hdi[rf_hdi$variable == "rewSens" & rf_hdi$t != "change", ]
# )
# broom.mixed::tidy(
#   as(mod_rewSens_indiv, "merModLmerTest"),
#   effects = "fixed",
#   conf.method = "boot",
#   conf.int = TRUE
# )

# # effort sensitivity
# mod_effSens_indiv <- lme4::lmer(
#   mean ~ factor(group, levels = c("CR", "BA")) * t + (1 | id),
#   weights = 1 / hdi_width,
#   data = rf_hdi[rf_hdi$variable == "effSens" & rf_hdi$t != "change", ]
# )
# broom.mixed::tidy(
#   as(mod_effSens_indiv, "merModLmerTest"),
#   effects = "fixed",
#   conf.method = "boot",
#   conf.int = TRUE
# )

# plot individual means and sds for each parameter, by group
# for t1 vs t2
# rew <- plot_indiv_pars(rf_hdi, "rewSens", "reward sensitivity", compl = FALSE, col_nums = c(1, 4))
# eff <- plot_indiv_pars(rf_hdi, "effSens", "effort sensitivity", compl = FALSE, col_nums = c(1, 4))

# rew_eff_indiv <- wrap_elements(
#   rew + eff +
#     plot_layout(nrow = 2, axis_title = "collect", guides = "collect") &
#     ggplot2::theme(
#       legend.position = "bottom",
#       legend.text = ggplot2::element_text(size = 18 * 1.2)
#     )
# ) &
#   plot_annotation(tag_levels = list("D")) &
#   ggplot2::theme(
#     plot.tag = ggplot2::element_text(family = "Open Sans", size = 28, face = "bold")
#   )

# non-completer vs completer baseline boxplots
# rew_eff_compl <-
#   plot_indiv_pars(
#     rf_hdi,
#     var = c("rewSens", "effSens"),
#     var_nm = c("**reward sensitivity**", "**effort sensitivity**"),
#     type = "box", tpt = "t1", sz_scale = c(2, 6),
#     clr = "completed", fnt_sz = 1.2, clr_levs = c(0, 1),
#     clr_labs = c("non-completers", "completers"),
#     compl = FALSE, col_nums = c(3, 6)
#   ) +
#   ggplot2::scale_y_continuous(
#     trans = scales::pseudo_log_trans(sigma = 0.7),
#     breaks = c(-1, 0, 1, 2, 4, 8)
#   )

# rew_eff_boxp_compl <- rew_eff_compl +
#   plot_layout(nrow = 1, guides = "collect") &
#   ggplot2::theme(legend.position = "bottom") &
#   plot_annotation(tag_levels = list("C")) &
#   ggplot2::theme(
#     plot.tag = ggplot2::element_text(family = "Open Sans", size = 28),
#     plot.tag.position = c(0, 1.02)
#   )

# ca_qn <- ca_hdi |>
#   dplyr::left_join(ids, by = c("id", "completed")) |>
#   dplyr::filter(variable == "delta_global_pos" & completed) |>
#   dplyr::left_join(qqnrs_t2, by = c("group", "completed", "subID")) |>
#   dplyr::select(id, group, mean, PHQ9_total_t2)

# ca_qn_ba <- ca_qn |> dplyr::filter(group == "BA")
# ca_qn_cr <- ca_qn |> dplyr::filter(group == "CR")

# cor(ca_qn$mean, ca_qn$PHQ9_total_t2, method = "pearson")
# cor(ca_qn_ba$mean, ca_qn_ba$PHQ9_total_t2, method = "pearson")
# cor(ca_qn_cr$mean, ca_qn_cr$PHQ9_total_t2, method = "pearson")

## heatmap of correlations between change in parameters and change in symptoms
symptom_corrs <- get_symptom_corrs(
  data_long = data_long_reff,
  model_type = "effort",
  draws = rf_i_draws,
  stan_ls = stan_ls_reff,
  qqnrs = qqnrs_t2,
  delta_par_pattern = "delta_(.*)_p\\[(.*)\\]",
  par_labels = c(
    "rewSens" = "reward sensitivity", "effSens" = "effort sensitivity"
  ),
  qns_to_invert = c("ERQCR_delta", "BADS_delta")
)
# print(symptom_corrs$panels$difference, n = 36)

ba_heatmap <- plot_correlation_heatmap(
  cor_summary_df = symptom_corrs$summary_df,
  int_group = "BA",
  par_labels = c(
    "rewSens" = "reward<br>sensitivity", "effSens" = "effort<br>sensitivity"
  ),
  qn_labels = c(
    "PHQ9" = "ΔPHQ-9", "DAS" = "ΔDAS", "miniSPIN" = "ΔminiSPIN",
    "BADS" = "ΔBADS", "ERQCR" = "ΔERQ-CR", "AMI" = "ΔAMI-BA"
  ),
  cor_lims = c(-0.35, 0.35),
  htmp_title = "behavioural activation",
  title_clr = "#88a0dc",
  label_type = "r_only",
  txt_sz = 7.5,
  fnt_sz = 1.2
)

cr_heatmap <- plot_correlation_heatmap(
  cor_summary_df = symptom_corrs$summary_df,
  int_group = "CR",
  par_labels = c(
    "rewSens" = "reward<br>sensitivity", "effSens" = "effort<br>sensitivity"
  ),
  qn_labels = c(
    "PHQ9" = "ΔPHQ-9", "DAS" = "ΔDAS", "miniSPIN" = "ΔminiSPIN",
    "BADS" = "ΔBADS", "ERQCR" = "ΔERQ-CR", "AMI" = "ΔAMI-BA"
  ),
  cor_lims = c(-0.35, 0.35),
  htmp_title = "cognitive restructuring",
  title_clr = "#ed968c",
  label_type = "r_only",
  txt_sz = 7.5,
  fnt_sz = 1.2
)

# pre-/post-intervention level correlations (parameter value vs. symptom
# score within the same timepoint) are available via symptom_corrs$summary_df
# too -- pass timepoint = "t1" / "t2" (default remains "change"), and drop the
# Δ from qn_labels since these aren't change scores, e.g.:
# ba_heatmap_t1 <- plot_correlation_heatmap(
#   cor_summary_df = symptom_corrs$summary_df,
#   int_group = "BA", timepoint = "t1",
#   par_labels = c("rewSens" = "reward<br>sensitivity", "effSens" = "effort<br>sensitivity"),
#   qn_labels = c(
#     "PHQ9" = "PHQ-9", "DAS" = "DAS", "miniSPIN" = "miniSPIN",
#     "BADS" = "BADS", "ERQCR" = "ERQ-CR", "AMI" = "AMI-BA"
#   ),
#   cor_lims = c(-0.35, 0.35),
#   htmp_title = "behavioural activation (timepoint 1)",
#   title_clr = "#88a0dc", label_type = "r_only", txt_sz = 6.5, fnt_sz = 1.2
# )

figure_3e <- wrap_elements(
  ba_heatmap + cr_heatmap +
    plot_layout(nrow = 2, guides = "collect") &
    ggplot2::theme(legend.position = "right") &
    plot_annotation(tag_levels = list("E")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
    )
)
figure_3abcde <- figure_3_abcd + figure_3e + plot_layout(nrow = 1, widths = c(0.69, 0.31))

## mediation analysis: arm -> parameter change -> PHQ9 --------------------
## z_sd_ref = "t2" puts PHQ-9 in units of the COMPLETE-CASE sd(t2), i.e. the same b_std the
## questionnaire ANCOVA forest is drawn on (qnr_analysis_ancova.R), so this model's direct_arm and
## total_arm can be checked against it below. NB not a pure relabel -- the .stan priors are fixed, so
## the divisor changes how informative they are, and fits on the old "pooled" scaling are not
## comparable. See prep_data()'s qnr_df block.
stan_ls_symptoms <- prep_data(
  data_long_reff, type = "effort", qnr_df = qqnrs_t2, excl_eff_task = "standard",
  quest_t1 = quest_t1, quest_t2 = quest_t2, excl_quest = "standard", z_sd_ref = "t2"
)

## joint mediation model: both mediators adjusted for each other, with an ANCOVA symptom submodel.
## the previous model (...-joint-decomp, no suffix) put a participant random intercept on BOTH PHQ-9
## rows, which imposes compound symmetry -- it forces Var(t1) == Var(t2), and PHQ-9 violates that
## here (sd 3.29 -> 3.86, variance ratio 1.37, because enrolment screened on baseline PHQ-9 but not
## on follow-up). The ANCOVA version conditions on the observed baseline instead. The CHOICE submodel
## is byte-identical between the two; see the new file's header for what that changes in the
## reportable quantities (one direct_arm instead of a per-arm pair, one total_arm, no prop_med_*).
joint_rew_med_mod_name <- "rewEff-multisess-ITT-symptoms-joint-decomp-ancova"
mod_rew_joint_med <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", joint_rew_med_mod_name, ".stan"))

fit_rew_joint_med <- mod_rew_joint_med$sample(
  data = stan_ls_symptoms,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
fit_rew_joint_med$save_object(file = paste0(pth_t2, "stan-fits/", joint_rew_med_mod_name, "-mediation-16000.rds"))
# rew_joint_med_rds <- paste0(pth_t2, "stan-fits/", joint_rew_med_mod_name, "-mediation-16000.rds")
# if (!file.exists(rew_joint_med_rds)) {
#   stop(
#     "no cached fit at ", rew_joint_med_rds, ".\n",
#     "The symptom submodel was reformulated as an ANCOVA, so the old ...-joint-decomp fit is stale ",
#     "and is deliberately NOT reused (its draws have direct_ba/direct_cr, which the ANCOVA does not ",
#     "define). Uncomment the $sample() block above and re-run.",
#     call. = FALSE
#   )
# }
# fit_rew_joint_med <- readRDS(rew_joint_med_rds)

joint_med_vars <- c(
  "alpha_phq", "beta_t1_phq", "beta_phq_group", "sigma_phq",
  "beta_phq_rewB", "beta_phq_rewB_group", "beta_phq_rewW", "beta_phq_rewW_group",
  "beta_phq_effB", "beta_phq_effB_group", "beta_phq_effW", "beta_phq_effW_group",
  "delta_rewSens_p", "delta_effSens_p",
  "a_rew_cr", "a_rew_ba", "a_eff_cr", "a_eff_ba",
  # observed-t2-only companions to the ITT a-paths: only 96 of 188 return, so the ITT a-path leans
  # heavily on the model-imputed t2 values, and the gap between the two is worth reporting
  "a_rew_cr_obs", "a_rew_ba_obs", "a_eff_cr_obs", "a_eff_ba_obs",
  "b_rew_cr", "b_rew_ba", "b_eff_cr", "b_eff_ba",
  # ANCOVA: ONE direct effect and ONE total, not a per-arm pair -- the outcome is the t2 LEVEL, so
  # there is no within-arm change term to split per arm, and prop_med_* no longer exists
  "direct_arm", "total_arm",
  "indirect_rew_cr", "indirect_rew_ba", "indirect_eff_cr", "indirect_eff_ba",
  "index_mod_med_rew", "index_mod_med_eff"
)
# confirm the fit actually carries every quantity read below, BEFORE any of it is used -- a stale
# cache or an edited .stan otherwise surfaces as a confusing error deep in a plotting call
missing_vars <- setdiff(sub("\\[.*", "", joint_med_vars), fit_rew_joint_med$metadata()$stan_variables)
if (length(missing_vars)) {
  stop("fit is missing generated quantities: ", paste(missing_vars, collapse = ", "), call. = FALSE)
}
fit_rew_joint_med$summary(variables = joint_med_vars)

draws_joint_med <- fit_rew_joint_med$draws(format = "df", variables = joint_med_vars) |>
  dplyr::mutate(
    phq_rewSens_t1_group_cr = beta_phq_rewB,
    phq_rewSens_t1_group_ba = beta_phq_rewB + beta_phq_rewB_group,
    phq_rewSens_t1_diff = beta_phq_rewB_group,
    phq_rewSens_t2_change_cr = beta_phq_rewW,
    phq_rewSens_t2_change_ba = beta_phq_rewW + beta_phq_rewW_group,
    phq_rewSens_t2_change_diff = beta_phq_rewW_group,
    phq_effSens_t1_group_cr = beta_phq_effB,
    phq_effSens_t1_group_ba = beta_phq_effB + beta_phq_effB_group,
    phq_effSens_t1_diff = beta_phq_effB_group,
    phq_effSens_t2_change_cr = beta_phq_effW,
    phq_effSens_t2_change_ba = beta_phq_effW + beta_phq_effW_group,
    phq_effSens_t2_change_diff = beta_phq_effW_group
  )

# difference in the within-person PHQ-9 <-> reward-sensitivity association between BA and CR
get_hdi_pd(draws_joint_med, "beta_phq_rewW_group")
get_hdi_pd(draws_joint_med, "beta_phq_effW_group")

## does reward-sensitivity change mediate symptom change differently in BA vs CR?
get_hdi_pd(draws_joint_med, "index_mod_med_rew")

## ...and does effort-sensitivity change still not mediate, adjusting for reward-sensitivity change?
get_hdi_pd(draws_joint_med, "index_mod_med_eff")

## direct / indirect / total decomposition, with the model's own total_arm checked against the
## mediator-free questionnaire ANCOVA. condition == 1 = BA here, so direct_arm and total_arm are
## already BA - CR, the same direction as the questionnaire forest -- no sign correction.
reff_med_decomp <- mediation_direct_total(
  draws_joint_med,
  suffixes = c("rew", "eff"),
  mediator_nms = c(rew = "reward sensitivity", eff = "effort sensitivity"),
  arm_dir = "BA - CR",
  out_path = "outputs/mediation_direct_total_reff.csv"
)
print(reff_med_decomp, n = Inf)

# plot the symptom decomposition panels for reward and effort sensitivity, side by side
rew_symp_decomp_plts <- plot_sympt_decomp(
  draws_joint_med,
  stan_ls = stan_ls_symptoms,
  prefix = "phq_rewSens",
  param_nm = "reward sensitivity",
  diff_lab = "behav. activ.\nvs. cogn. restr.",
  flip_order = TRUE,
  col_nums = c(1, 4), fnt_sz = 1.6
)
eff_symp_decomp_plts <- plot_sympt_decomp(
  draws_joint_med,
  stan_ls = stan_ls_symptoms,
  prefix = "phq_effSens",
  param_nm = "effort sensitivity",
  diff_lab = "behav. activ.\nvs. cogn. restr.",
  flip_order = TRUE,
  suppress_labels = FALSE,
  col_nums = c(1, 4), fnt_sz = 1.6
)

rew_symp_change <- wrap_elements(
  wrap_elements(rew_symp_decomp_plts$scatter_panel) +
    wrap_elements(
      rew_symp_decomp_plts$within_panel +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 0, b = 0, l = 0, r = 0))
    ) +
    plot_layout(widths = c(0.45, 0.55)) &
    plot_annotation(tag_levels = list("F")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", face = "bold", size = 28)
    )
)
eff_symp_change <- wrap_elements(
  wrap_elements(eff_symp_decomp_plts$scatter_panel) +
    wrap_elements(
      eff_symp_decomp_plts$within_panel +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 0, b = 0, l = 0, r = 0))
    ) +
    plot_layout(widths = c(0.45, 0.55)) &
    plot_annotation(tag_levels = list("G")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", face = "bold", size = 28)
    )
)
figure_3fg <- wrap_elements(
  rew_symp_change + eff_symp_change + plot_layout(nrow = 1)
)

reff_mediator_cols <- stats::setNames(
  MetBrewer::met.brewer("Hokusai2")[c(1, 5)], c("rew", "eff")
)

## plot the two-mediator path diagram
joint_med_plts <- plot_joint_mediation(
  draws_joint_med,
  mediator1_suffix = "rew", mediator2_suffix = "eff",
  mediator1_nm = "reward sensitivity",
  mediator2_nm = "effort sensitivity",
  mediator_cols = reff_mediator_cols,
  diff_lab = "behav. activ.\nvs. cogn. restr.",
  flip_order = TRUE,
  col_nums = c(1, 4), fnt_sz = 1.5
)
## inset the a1/b1 (upper, reward) and a2/b2 (lower, effort) path plots by their
## respective arrows, the residual (c') path plot in the slot freed up above its
## now-nudged-down line, then stack the combined indirect-effect panel underneath
joint_med_diagram <- wrap_elements(
  joint_med_plts$diagram +
    patchwork::inset_element(
      joint_med_plts$a1_path, left = 0.025, bottom = 0.67, right = 0.3, top = 0.98, align_to = "panel"
    ) +
    patchwork::inset_element(
      joint_med_plts$b1_path, left = 0.7, bottom = 0.67, right = 0.975, top = 0.98, align_to = "panel"
    ) +
    patchwork::inset_element(
      joint_med_plts$a2_path, left = 0.025, bottom = 0.02, right = 0.3, top = 0.33, align_to = "panel"
    ) +
    patchwork::inset_element(
      joint_med_plts$b2_path, left = 0.7, bottom = 0.02, right = 0.975, top = 0.33, align_to = "panel"
    ) +
    patchwork::inset_element(
      joint_med_plts$direct_path, left = 0.34, bottom = 0.43, right = 0.6, top = 0.72, align_to = "panel"
    )
)

figure_3h <- wrap_elements(
  joint_med_diagram + wrap_elements(joint_med_plts$indirect) +
    plot_layout(widths = c(0.64, 0.36)) +
    plot_annotation(tag_levels = list("H")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", face = "bold", size = 28)
    )
)

# Figure 3 all together
fig_3 <-  figure_3abcde / figure_3fg / figure_3h + plot_layout(nrow = 3, heights = c(0.48, 0.22, 0.3))
fig_3

# non-completer vs completers at baseline --------------------------------------

## baseline symptom/questionnaire scores, completers vs. non-completers
ncpl_symp_meta <- tibble::tribble(
  ~measure,   ~col,
  "PHQ-9",    "PHQ9_total",
  "DAS",      "DAS_total",
  "miniSPIN", "miniSPIN_total",
  "BADS",     "BADS_total",
  "AMI-BA",   "AMI_behavActiv",
  "ERQ-CR",   "ERQCR_total"
)
ncpl_symp_df <- quest_t1 |>
  dplyr::mutate(non_completer = !subID %in% quest_t2$subID) |>
  dplyr::select(subID, non_completer, tidyselect::all_of(ncpl_symp_meta$col)) |>
  tidyr::pivot_longer(
    cols = tidyselect::all_of(ncpl_symp_meta$col), names_to = "col", values_to = "score"
  ) |>
  dplyr::left_join(ncpl_symp_meta, by = "col") |>
  dplyr::mutate(measure = factor(measure, levels = ncpl_symp_meta$measure)) |>
  dplyr::filter(!is.na(score)) |>
  dplyr::group_by(measure) |>
  dplyr::mutate(score_z = as.numeric(scale(score))) |>
  dplyr::ungroup()

figure_5a_symptoms <- wrap_elements(
  plot_ncpl_symptom_bars(
    ncpl_symp_df,
    measure_col = "measure", value_col = "score_z", completer_col = "non_completer",
    error_type = "se", fnt_sz = 1.6, col_nums = c(2, 5), legend_pos = "none"
  ) +
    plot_annotation(tag_levels = list("A")) &
    ggplot2::theme(plot.tag = ggplot2::element_text(family = "Open Sans", size = 28))
)
## NB: reward-effort's figure_5ab and this script's figure_5cd already claim tags
## A-D below -- re-letter those (or drop this panel's own "A" tag) once you've decided
## where figure_5a_symptoms slots into the assembled fig_5.

ncpl_rew_bsl <- reff_plot(
  mf_df,
  plot_type = "reward",
  mode = "completer_status",
  label_colour = "#b75347",
  ln_width = 1.2,
  err_width = 0.3,
  dodge_width = 0.45,
  point_size = 3.5,
  shape_types = c(25, 23),
  legend_pos = "bottom",
  legend_fnt_sc = 14,
  fnt_sz = 1.5,
  brew_col = "Archambault",
  col_nums = c(2, 5)
)
ncpl_eff_bsl <- reff_plot(
  mf_df,
  plot_type = "effort",
  mode = "completer_status",
  label_colour = "#2f70a1",
  ln_width = 1.2,
  err_width = 4,
  dodge_width = 6,
  point_size = 3.5,
  shape_types = c(25, 23),
  legend_pos = "bottom",
  legend_fnt_sc = 14,
  fnt_sz = 1.5,
  brew_col = "Archambault",
  col_nums = c(2, 5)
)

mf_ncpl_reff <- wrap_elements(
  ncpl_rew_bsl + ncpl_eff_bsl + plot_layout(ncol = 2, guides = "collect") &
    plot_annotation(tag_levels = list(c("B", ""))) &
    ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
    )
)

rew_ncpl_plt <- plot_ncpl_diff(
  draws_rf_itt,
  param = "rewSens",
  param_nm = "reward sensitivity",
  diff_lab = "completers vs.\nnon-completers",
  col_nums = c(2, 5, 3),
  fnt_sz = 1.6
) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())

eff_ncpl_plt <- plot_ncpl_diff(
  draws_rf_itt,
  param = "effSens",
  param_nm = "effort sensitivity",
  diff_lab = "non-completers\nvs. completers",
  flip_order = TRUE,
  suppress_labels = TRUE,
  col_nums = c(2, 5, 3),
  fnt_sz = 1.6
) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())

ncpl_reff <- wrap_elements(
  rew_ncpl_plt + eff_ncpl_plt +
    plot_layout(ncol = 2, widths = c(0.55, 0.45)) &
    plot_annotation(
      title = "session 1: non-completers vs. completers",
      tag_levels = list("C")
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = "Open Sans SemiBold", size = 20, colour = "grey30", margin = ggplot2::margin(l = 40, t = 12)
      ),
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28, vjust = -0.8)
    )
)

figure_5abc <- wrap_elements(
  figure_5a_symptoms / mf_ncpl_reff / ncpl_reff +
    plot_layout(nrow = 3, heights = c(0.32, 0.38, 0.3))
)

# SUPPLEMENTARY MATERIALS =======================================================
## sensitivity analyses --------------------------------------------------------
stan_ls_reff_cca <- prep_data(data_long_reff, type = "effort", cca = TRUE)
rf_cca <- "rewEff-multisess-CCA-intervention"
mod_rf_cca <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_cca, ".stan"))

## fit models with cmdstan
fit_rf_cca <- mod_rf_cca$sample(
  data = stan_ls_reff_cca,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
# fit_rf_cca$save_object(file = paste0(pth_t2, "stan-fits/", rf_cca, "-16000.rds"))
# fit_rf_cca <- readRDS(paste0(pth_t2, "stan-fits/", rf_cca, "-16000.rds"))

rf_cca_vars <- c(
  "mu_effSens[1]", "mu_effSens[2]", "mu_rewSens[1]", "mu_rewSens[2]",
  "effSens_int", "rewSens_int", "delta_effSens", "delta_rewSens"
)
draws_rf_cca <- fit_rf_cca$draws(format = "df", variables = rf_cca_vars)
long_draws_rf_cca <- draws_rf_cca |> tidybayes::gather_draws(!!!rlang::syms(rf_cca_vars))
tidy_rf_cca <- extract_draws(draws_rf_cca, long_draws_rf_cca, hdi = 0.95)
tidy_rf_cca$hdi

# rf_cca_i_draws <- fit_rf_cca$draws(format = "df", variables = rf_i_pars)
# rf_cca_hdi <- extract_indiv_pars(rf_cca_i_draws, stan_ls_reff_cca)

## Less strict ITT/non-completion assumptions
rf_strat <- "rewEff-multisess-ITT-intervention-stratified"
mod_rf_strat <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_strat, ".stan"))

## fit models with cmdstan
fit_rf_strat <- mod_rf_strat$sample(
  data = stan_ls_reff,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
# fit_rf_strat$save_object(file = paste0(pth_t2, "stan-fits/", rf_strat, "-16000.rds"))
# fit_rf_strat <- readRDS(paste0(pth_t2, "stan-fits/", rf_strat, "-16000.rds"))
# fit_rf_strat$diagnostic_summary()
fit_rf_strat_loo <- fit_rf_strat$loo(variables = "log_lik")
# saveRDS(fit_rf_strat_loo, file = paste0(pth_t2, "stan-fits/", rf_strat, "-loo.rds"))

draws_rf_strat <- fit_rf_strat$draws(format = "df", variables = rf_vars) |>
  dplyr::rename(
    effSens_int = effSens_int_c, # the estimated effect in completers
    rewSens_int = rewSens_int_c
  )
long_draws_rf_strat <- draws_rf_strat |> tidybayes::gather_draws(!!!rlang::syms(rf_vars))
tidy_rf_strat <- extract_draws(draws_rf_strat, long_draws_rf_strat, hdi = 0.95)
tidy_rf_strat$hdi

## Alternative parameterisation with mean and difference
rf_delta <- "rewEff-multisess-ITT-intervention-delta"
mod_rf_delta <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_delta, ".stan"))

## fit models with cmdstan
fit_rf_delta <- mod_rf_delta$sample(
  data = stan_ls_reff,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
# fit_rf_delta$save_object(file = paste0(pth_t2, "stan-fits/", rf_delta, "-16000.rds"))
# fit_rf_delta$diagnostic_summary()
fit_rf_delta_loo <- fit_rf_delta$loo(variables = "log_lik")
# saveRDS(fit_rf_delta_loo, file = paste0(pth_t2, "stan-fits/", rf_delta, "-loo.rds"))

draws_rf_delta <- fit_rf_delta$draws(format = "df", variables = rf_vars)
long_draws_rf_delta <- draws_rf_delta |> tidybayes::gather_draws(!!!rlang::syms(rf_vars))
tidy_rf_delta <- extract_draws(draws_rf_delta, long_draws_rf_delta, hdi = 0.95)
tidy_rf_delta$hdi

## Worksheet quality
# ratings for in-study completers only (the scrambled file is already filtered to them)
qc_all <- read_stand_data(pth_t2, "worksheet_qc") |>
  dplyr::mutate(
    dplyr::across(tidyselect::where(is.numeric), ~tidyr::replace_na(., 0))
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    qc_score = ifelse(
      intervention == "BA",
      activ_quality + prac_quality + 0.5 * (reliv_quality + smart_quality),
      check_quality + rework_quality + fair_quality + 0.5 * (digg_quality + win_quality)
    )
  ) |>
  dplyr::ungroup()

## Including worksheet quality
stan_ls_reff_qc <- prep_data(data_long_reff, type = "effort", qc = qc_all)
rf_itt_qc <- "rewEff-multisess-ITT-intervention-additive-qc"
mod_rf_itt_qc <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_itt_qc, ".stan"))

## fit models with cmdstan
fit_rf_itt_qc <- mod_rf_itt_qc$sample(
  data = stan_ls_reff_qc,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
# fit_rf_itt_qc$save_object(file = paste0(pth_t2, "stan-fits/", rf_itt_qc, "-16000.rds"))
# fit_rf_itt_qc$diagnostic_summary()
fit_rf_qc_loo <- fit_rf_itt_qc$loo(variables = "log_lik")
# saveRDS(fit_rf_qc_loo, file = paste0(pth_t2, "stan-fits/", rf_itt_qc, "-loo.rds"))

draws_rf_qc <- fit_rf_itt_qc$draws(format = "df", variables = rf_vars)
long_draws_rf_qc <- draws_rf_qc |> tidybayes::gather_draws(!!!rlang::syms(rf_vars))
tidy_rf_qc <- extract_draws(draws_rf_qc, long_draws_rf_qc, hdi = 0.95)
tidy_rf_qc$hdi

## check influence of inattentive responders -----------------------------------
### we have the effort task catch trials in data_long_reff + the attention checks in self reports
### 1. stricter effort task exclusion - any one catch trial failed
stan_ls_reff_strict <- prep_data(data_long_reff, type = "effort", excl_eff_task = "strict")
rf_itt <- "rewEff-multisess-ITT-intervention-additive"
mod_rf_itt <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_itt, ".stan"))

fit_rf_itt_strict <- mod_rf_itt$sample(
  data = stan_ls_reff_strict,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
# fit_rf_itt_strict$save_object(file = paste0(pth_t2, "stan-fits/", rf_itt, "-strict-16000.rds"))
# fit_rf_itt_strict <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt, "-strict-16000.rds"))

draws_rf_itt_strict <- fit_rf_itt_strict$draws(format = "df", variables = rf_vars)
long_draws_rf_itt_strict <- draws_rf_itt_strict |> tidybayes::gather_draws(!!!rlang::syms(rf_vars))
tidy_rf_itt_strict <- extract_draws(draws_rf_itt_strict, long_draws_rf_itt_strict, hdi = 0.95)
tidy_rf_itt_strict$hdi

# 2. strictest - also exclude those who failed any attention checks in questionnaires
stan_ls_reff_strictest <- prep_data(
  data_long_reff, type = "effort", excl_eff_task = "strict", excl_quest = "strict",
  quest_t1 = quest_t1, quest_t2 = quest_t2
)
fit_rf_itt_strictest <- mod_rf_itt$sample(
  data = stan_ls_reff_strictest,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
# fit_rf_itt_strictest$save_object(file = paste0(pth_t2, "stan-fits/", rf_itt, "-strictest-16000.rds"))
# fit_rf_itt_strictest <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt, "-strictest-16000.rds"))

draws_rf_itt_strictest <- fit_rf_itt_strictest$draws(format = "df", variables = rf_vars)
long_draws_rf_itt_strictest <- draws_rf_itt_strictest |> tidybayes::gather_draws(!!!rlang::syms(rf_vars))
tidy_rf_itt_strictest <- extract_draws(draws_rf_itt_strictest, long_draws_rf_itt_strictest, hdi = 0.95)
tidy_rf_itt_strictest$hdi
## individual item level associations between parameters and symptoms ----------
# phq9_items <- qqnrs_t2 |>
#   dplyr::select(
#     subID, group, completed, dplyr::matches("PHQ9_[1-9]_t[12]$")
#   ) |>
#   # make new columns for delta scores for each item by subtracting t1 from t2
#   dplyr::mutate(
#     dplyr::across(
#       tidyselect::matches("PHQ9_[1-9]_t2$"),
#       ~ . - get(stringr::str_replace(dplyr::cur_column(), "_t2$", "_t1")),
#       .names = "{stringr::str_replace(.col, '_t2$', '_delta')}"
#     )
#   ) |>
#   dplyr::select(subID, group, completed, dplyr::matches("PHQ9_[1-9]_delta$"))
# indiv_symptom_corrs <- get_symptom_corrs(
#   data_long = data_long_reff,
#   draws = rf_i_draws,
#   stan_ls = stan_ls_reff,
#   qqnrs = phq9_items,
#   delta_par_pattern = "delta_(.*)_p\\[(.*)\\]",
#   par_labels = c(
#     "rewSens" = "reward sensitivity", "effSens" = "effort sensitivity"
#   )
# )
# print(indiv_symptom_corrs$summary_df, n = 36)
# ba_phq9_heatmap <- plot_correlation_heatmap(
#   cor_summary_df = indiv_symptom_corrs$summary_df,
#   int_group = "BA",
#   par_labels = c(
#     "rewSens" = "reward sensitivity", "effSens" = "effort sensitivity"
#   ),
#   qn_labels = c(
#     "PHQ9_1" = "ΔPHQ-9 item 1", "PHQ9_2" = "ΔPHQ-9 item 2", "PHQ9_3" = "ΔPHQ-9 item 3",
#     "PHQ9_4" = "ΔPHQ-9 item 4", "PHQ9_5" = "ΔPHQ-9 item 5", "PHQ9_6" = "ΔPHQ-9 item 6",
#     "PHQ9_7" = "ΔPHQ-9 item 7", "PHQ9_8" = "ΔPHQ-9 item 8", "PHQ9_9" = "ΔPHQ-9 item 9"
#   ),
#   cor_lims = c(-0.2, 0.2),
#   htmp_title = "goal-setting",
#   title_clr = "#88a0dc"
# ) + ggplot2::theme(legend.position = "none")
# cr_phq9_heatmap <- plot_correlation_heatmap(
#   cor_summary_df = indiv_symptom_corrs$summary_df,
#   int_group = "CR",
#   par_labels = c(
#     "rewSens" = "reward sensitivity", "effSens" = "effort sensitivity"
#   ),
#   qn_labels = c(
#     "PHQ9_1" = "ΔPHQ-9 item 1", "PHQ9_2" = "ΔPHQ-9 item 2", "PHQ9_3" = "ΔPHQ-9 item 3",
#     "PHQ9_4" = "ΔPHQ-9 item 4", "PHQ9_5" = "ΔPHQ-9 item 5", "PHQ9_6" = "ΔPHQ-9 item 6",
#     "PHQ9_7" = "ΔPHQ-9 item 7", "PHQ9_8" = "ΔPHQ-9 item 8", "PHQ9_9" = "ΔPHQ-9 item 9"
#   ),
#   cor_lims = c(-0.2, 0.2),
#   htmp_title = "restructuring",
#   title_clr = "#ed968c"
# ) + ggplot2::theme(axis.text.y = ggplot2::element_blank())

## models including accept bias parameter --------------------------------------
stan_ls_reff_ab <- prep_data(
  data_long_reff, type = "effort", excl_eff_task = "standard",
  quest_t1 = quest_t1, quest_t2 = quest_t2, excl_quest = "standard",
  recode_effort = TRUE # make it so that choice==1 is accepting the high effort option, rather than just route 2
)

rf_itt_ab <- "rewEff-multisess-ITT-intervention-acceptbias"
mod_rf_itt_ab <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_itt_ab, ".stan"))

fit_rf_itt_ab <- mod_rf_itt_ab$sample(
  data = stan_ls_reff_ab,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
fit_rf_itt_ab$save_object(file = paste0(pth_t2, "stan-fits/", rf_itt_ab, "-16000.rds"))
# fit_rf_itt_ab <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt_ab, "-16000.rds"))
fit_rf_itt_ab_loo <- fit_rf_itt_ab$loo(variables = "log_lik", cores = 4)
# saveRDS(fit_rf_itt_ab_loo, file = paste0(pth_t2, "stan-fits/", rf_itt_ab, "-loo.rds"))

# fit_rf_loo <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt, "-loo.rds"))
# fit_rf_qc_loo <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt_qc, "-loo.rds"))
# fit_rf_delta_loo <- readRDS(paste0(pth_t2, "stan-fits/", rf_delta, "-loo.rds"))
# fit_rf_strat_loo <- readRDS(paste0(pth_t2, "stan-fits/", rf_strat, "-loo.rds"))
# fit_rf_itt_ab_loo <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt_ab, "-loo.rds"))

# does it fit better?
mc_rf <- loo::loo_compare(
  list(
    "baseline" = fit_rf_loo,
    "qc" = fit_rf_qc_loo,
    "delta" = fit_rf_delta_loo,
    "stratified" = fit_rf_strat_loo,
    "acceptbias" = fit_rf_itt_ab_loo
  )
)
print(mc_rf, simplify = FALSE)

fit_rf_itt_ab$summary(variables = "effSens_int")

rf_ab_vars <- c(
  "mu_effSens[1]", "mu_effSens[2]", "mu_rewSens[1]", "mu_rewSens[2]", "mu_alpha[1]", "mu_alpha[2]",
  "effSens_int", "rewSens_int", "alpha_int", "effSens_ncpl", "rewSens_ncpl", "alpha_ncpl",
  "delta_effSens", "delta_rewSens", "delta_alpha"
)
draws_rf_itt_ab <- fit_rf_itt_ab$draws(format = "df", variables = rf_ab_vars)
long_draws_rf_ab <- draws_rf_itt_ab |> tidybayes::gather_draws(!!!rlang::syms(rf_ab_vars))
tidy_rf_ab <- extract_draws(draws_rf_itt_ab, long_draws_rf_ab, hdi = 0.95)
tidy_rf_ab$hdi

pd_es_ab <- tidy_rf_ab$df |> dplyr::filter(.variable == "effSens_int") |> dplyr::pull(.value)
pd_ncpl_rs_ab <- tidy_rf_ab$df |> dplyr::filter(.variable == "rewSens_ncpl") |> dplyr::pull(.value)
pd_ncpl_alpha_ab <- tidy_rf_ab$df |> dplyr::filter(.variable == "alpha_ncpl") |> dplyr::pull(.value)
bayestestR::p_direction(pd_es_ab, method = "direct")
bayestestR::p_direction(pd_ncpl_rs_ab, method = "direct")
bayestestR::p_direction(pd_ncpl_alpha_ab, method = "direct")

## association with symptoms
rf_ab_i_pars <- c(
  "effSens", "delta_effSens_p", "rewSens", "delta_rewSens_p", "alpha", "delta_alpha_p"
)
rf_ab_i_draws <- fit_rf_itt_ab$draws(format = "df", variables = rf_ab_i_pars)
rf_ab_hdi <- extract_indiv_pars(rf_ab_i_draws, stan_ls_reff_ab)

# suppl fig 3a - plot of densities

rew_sens_ab_plts <- plot_reff_change(
  draws_rf_itt_ab,
  diff_lab = "behav. activ.\nvs. cogn. restr.",
  flip_order = TRUE,
  param = "rewSens",
  param_nm = "reward sensitivity",
  col_nums = c(1, 4),
  fnt_sz = 1.6
)
eff_sens_ab_plts <- plot_reff_change(
  draws_rf_itt_ab,
  diff_lab = "behav. activ.\nvs. cogn. restr.",
  flip_order = TRUE,
  suppress_labels = TRUE,
  param = "effSens",
  param_nm = "effort sensitivity",
  col_nums = c(1, 4),
  fnt_sz = 1.6
)
accept_bias_ab_plts <- plot_reff_change(
  draws_rf_itt_ab,
  diff_lab = "behav. activ.\nvs. cogn. restr.",
  flip_order = TRUE,
  suppress_labels = TRUE,
  param = "alpha",
  param_nm = "acceptance bias",
  col_nums = c(1, 4),
  fnt_sz = 1.6
)

rew_sens_change_ab <- wrap_elements(
  rew_sens_ab_plts$change_panel &
    plot_annotation(tag_levels = list(c("A", ""))) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
    )
)
eff_sens_change_ab <- wrap_elements(
  eff_sens_ab_plts$change_panel &
    plot_annotation(tag_levels = list("B")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28, hjust = 2)
    )
)
accept_bias_change_ab <- wrap_elements(
  accept_bias_ab_plts$change_panel &
    plot_annotation(tag_levels = list("C")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28, hjust = 2)
    )
)

# figure s2abc - plot of densities and change scores
figure_s3abc <- wrap_elements(
  rew_sens_change_ab + eff_sens_change_ab + accept_bias_change_ab +
    plot_layout(ncol = 3, widths = c(0.4, 0.3, 0.3))
)

# non-completion effects for acceptance bias model
rew_ab_ncpl_plt <- plot_ncpl_diff(
  draws = draws_rf_itt_ab,
  param = "rewSens",
  param_nm = "reward sensitivity",
  diff_lab = "non-completers\nvs. completers",
  flip_order = TRUE,
  col_nums = c(2, 5, 3),
  fnt_sz = 1.6
) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())

eff_ab_ncpl_plt <- plot_ncpl_diff(
  draws = draws_rf_itt_ab,
  param = "effSens",
  param_nm = "effort sensitivity",
  diff_lab = "non-completers\nvs. completers",
  flip_order = TRUE,
  suppress_labels = TRUE,
  col_nums = c(2, 5, 3),
  fnt_sz = 1.6
) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())

alpha_ab_ncpl_plt <- plot_ncpl_diff(
  draws = draws_rf_itt_ab,
  param = "alpha",
  param_nm = "acceptance bias",
  diff_lab = "non-completers\nvs. completers",
  flip_order = TRUE,
  suppress_labels = TRUE,
  col_nums = c(2, 5, 3),
  fnt_sz = 1.6
) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())

figure_s3f <- wrap_elements(
  rew_ab_ncpl_plt + eff_ab_ncpl_plt + alpha_ab_ncpl_plt +
    plot_layout(nrow = 1) &
    plot_annotation(
      title = "session 1: non-completers vs. completers",
      tag_levels = list("F")
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = "Open Sans SemiBold", size = 20, colour = "grey30", margin = ggplot2::margin(l = 40, t = 12)
      ),
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28, vjust = -0.8)
    )
)

symptom_corrs_ab <- get_symptom_corrs(
  data_long = data_long_reff,
  model_type = "effort",
  draws = rf_ab_i_draws,
  stan_ls = stan_ls_reff_ab,
  qqnrs = qqnrs_t2,
  delta_par_pattern = "delta_(.*)_p\\[(.*)\\]",
  par_labels = c(
    "rewSens" = "reward sensitivity",
    "effSens" = "effort sensitivity",
    "alpha" = "acceptance bias"
  ),
  qns_to_invert = c("ERQCR_delta", "BADS_delta")
)

# plot
ab_ba_heatmap <- plot_correlation_heatmap(
  cor_summary_df = symptom_corrs_ab$summary_df,
  int_group = "BA",
  par_labels = c(
    "rewSens" = "reward<br>sensitivity",
    "effSens" = "effort<br>sensitivity",
    "alpha" = "acceptance<br>bias"
  ),
  qn_labels = c(
    "PHQ9" = "ΔPHQ-9", "DAS" = "ΔDAS", "miniSPIN" = "ΔminiSPIN",
    "BADS" = "ΔBADS", "ERQCR" = "ΔERQ-CR", "AMI" = "ΔAMI-BA"
  ),
  cor_lims = c(-0.35, 0.35),
  htmp_title = "behavioural activation",
  label_type = "r_only",
  title_clr = "#88a0dc",
  txt_sz = 6,
  fnt_sz = 1
)
ab_cr_heatmap <- plot_correlation_heatmap(
  cor_summary_df = symptom_corrs_ab$summary_df,
  int_group = "CR",
  par_labels = c(
    "rewSens" = "reward<br>sensitivity",
    "effSens" = "effort<br>sensitivity",
    "alpha" = "acceptance<br>bias"
  ),
  qn_labels = c(
    "PHQ9" = "ΔPHQ-9", "DAS" = "ΔDAS", "miniSPIN" = "ΔminiSPIN",
    "BADS" = "ΔBADS", "ERQCR" = "ΔERQ-CR", "AMI" = "ΔAMI-BA"
  ),
  cor_lims = c(-0.35, 0.35),
  htmp_title = "cognitive restructuring",
  label_type = "r_only",
  title_clr = "#ed968c",
  txt_sz = 6,
  fnt_sz = 1
) + ggplot2::theme(axis.text.y = ggplot2::element_blank())

# parameter x symptom correlation heatmaps
figure_s3_corr <- wrap_elements(
  ab_ba_heatmap + ab_cr_heatmap +
    plot_layout(ncol = 2, guides = "collect") &
    ggplot2::theme(legend.position = "right") &
    plot_annotation(tag_levels = list(c("D", ""))) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
    )
)

# Shapley decompositions
dec_std <- reff_delta_decomp(fit_rf_itt, stan_ls_reff)
dec_ab  <- reff_delta_decomp(fit_rf_itt_ab, stan_ls_reff_ab)

figure_s3_decomp <- wrap_elements(
  plot_delta_decomp_compare(
    dec_std, dec_ab,
    main_label = "primary model",
    ext_label = "extended model",
    base_size = 20,
    include_total = TRUE
  ) +
    plot_annotation(tag_levels = list("E")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
    )
)

figure_s3_de <- wrap_elements(figure_s3_corr + figure_s3_decomp + plot_layout(widths = c(0.54, 0.46)))

# make supplementary figure 3: group-level acceptance-bias results
suppl_fig_3 <- figure_s3abc / figure_s3_de / figure_s3f + plot_layout(nrow = 3, heights = c(0.3, 0.4, 0.3))
suppl_fig_3

## parameter recovery for baseline and acceptance bias models ------------------
# config
n_t1 <- 200
n_t2 <- 120
sd_rew <- 1.5
sd_eff <- 3
sd_alpha <- 1.5

iter_warmup <- 2000
iter_sampling <- 4000
chains <- 4

# seeds for reproducible parameter recovery
seed_base <- 123
seed_task_bsl <- seed_base + 1
seed_task_ab <- seed_base + 2
seed_params_bsl <- seed_base + 3
seed_params_ab <- seed_base + 4
seed_sim_bsl <- seed_base + 5
seed_sim_ab <- seed_base + 6
seed_fit_bsl <- seed_base + 7
seed_fit_ab <- seed_base + 8

# data lists
stan_ls_reff <- prep_data(
  data_long_reff,
  type = "effort",
  excl_eff_task = "standard",
  quest_t1 = quest_t1,
  quest_t2 = quest_t2,
  excl_quest = "standard"
)

stan_ls_reff_ab <- prep_data(
  data_long_reff,
  type = "effort",
  excl_eff_task = "standard",
  quest_t1 = quest_t1,
  quest_t2 = quest_t2,
  excl_quest = "standard",
  recode_effort = TRUE # make 1 = high effort option chosen
)
# saveRDS(stan_ls_reff_ab, file = paste0(pth_t2, "stan-fits/stan_ls_reff_ab_sim.rds"))

# get empirical means of parameters to use as centre of sampling distributions for parameter recovery

rf_itt <- "rewEff-multisess-ITT-intervention-additive"
rf_itt_ab <- "rewEff-multisess-ITT-intervention-acceptbias"
# fit_rf_itt <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt, "-16000.rds"))
# fit_rf_itt_ab <- readRDS(paste0(pth_t2, "stan-fits/", rf_itt_ab, "-16000.rds"))

emp_means_bsl <- get_empirical_means_from_fit(fit_rf_itt, stan_ls_reff, include_alpha = FALSE)
emp_means_ab <- get_empirical_means_from_fit(fit_rf_itt_ab, stan_ls_reff_ab, include_alpha = TRUE)

# parameter recovery: baseline model
task_template_bsl <- sample_task_structure(data_long_reff, n_t1 = n_t1, n_t2 = n_t2, seed = seed_task_bsl)
task_template_ab <- sample_task_structure(data_long_reff, n_t1 = n_t1, n_t2 = n_t2, seed = seed_task_ab)

params_bsl <- sample_param_draws(
  task_template_bsl,
  mean_rew = emp_means_bsl$rew,
  mean_eff = emp_means_bsl$eff,
  sd_rew = sd_rew,
  sd_eff = sd_eff,
  include_alpha = FALSE,
  seed = seed_params_bsl
)

sim_bsl <- simulate_rew_eff_from_template(
  template_df = task_template_bsl,
  params_df = params_bsl,
  include_alpha = FALSE,
  recode_effort = FALSE,
  seed = seed_sim_bsl
)

stan_ls_sim_bsl <- prep_data(
  sim_bsl,
  type = "effort",
  excl_eff_task = "none",
  excl_quest = "none",
  recode_effort = FALSE
)

mod_bsl_sim <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_itt, ".stan"))
fit_bsl_sim <- mod_bsl_sim$sample(
  data = stan_ls_sim_bsl,
  seed = seed_fit_bsl,
  chains = chains,
  parallel_chains = chains,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
# fit_bsl_sim$save_object(file = paste0(pth_t2, "stan-fits/sim-bsl-fit.rds"))
# fit_bsl_sim <- readRDS(paste0(pth_t2, "stan-fits/sim-bsl-fit.rds"))

bsl_recovery <- plot_reff_recovery(
  sim_params = params_bsl,
  sim_data = sim_bsl,
  sim_fit = fit_bsl_sim,
  stan_ls_sim = stan_ls_sim_bsl,
  vars = c("rewSens", "effSens"),
  recode_effort = FALSE
)

## parameter recovery: acceptance-bias model

params_ab <- sample_param_draws(
  task_template_ab,
  mean_rew = emp_means_ab$rew,
  mean_eff = emp_means_ab$eff,
  mean_alpha = emp_means_ab$alpha,
  sd_rew = sd_rew,
  sd_eff = sd_eff,
  sd_alpha = sd_alpha,
  include_alpha = TRUE,
  seed = seed_params_ab
)
# saveRDS(params_ab, file = paste0(pth_t2, "stan-fits/sim_params_ab.rds"))

sim_ab <- simulate_rew_eff_from_template(
  template_df = task_template_ab,
  params_df = params_ab,
  include_alpha = TRUE,
  recode_effort = TRUE,
  seed = seed_sim_ab
)

stan_ls_sim_ab <- prep_data(
  sim_ab,
  type = "effort",
  excl_eff_task = "none",
  excl_quest = "none",
  recode_effort = TRUE
)

# recovery_data <- list(
#   bsl = list(
#     params = params_bsl,
#     sim_data = sim_bsl,
#     stan_ls_sim = stan_ls_sim_bsl
#   ),
#   ab = list(
#     params = params_ab,
#     sim_data = sim_ab,
#     stan_ls_sim = stan_ls_sim_ab
#   )
# )
# saveRDS(recovery_data, file = paste0(pth_t2, "stan-fits/sim_recovery_data.rds"))

mod_ab_sim <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/rew-eff/", rf_itt_ab, ".stan"))
fit_ab_sim <- mod_ab_sim$sample(
  data = stan_ls_sim_ab,
  seed = seed_fit_ab,
  chains = chains,
  parallel_chains = chains,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
# fit_ab_sim$save_object(file = paste0(pth_t2, "stan-fits/sim-ab-fit.rds"))
# fit_ab_sim <- readRDS(paste0(pth_t2, "stan-fits/sim-ab-fit.rds"))

ab_recovery <- plot_reff_recovery(
  sim_params = params_ab,
  sim_data = sim_ab,
  sim_fit = fit_ab_sim,
  stan_ls_sim = stan_ls_sim_ab,
  vars = c("rewSens", "effSens", "alpha"),
  recode_effort = TRUE
)

suppl_fig_1 <- wrap_elements(
  wrap_elements(
    bsl_recovery$scatter + ab_recovery$scatter +
      plot_layout(ncol = 2, widths = c(0.4, 0.6))  &
      plot_annotation(tag_levels = list(c("A", "B"))) &
      ggplot2::theme(
        plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
      )
  ) +
    wrap_elements(
      bsl_recovery$heatmap + ab_recovery$heatmap +
        plot_layout(ncol = 2, widths = c(0.4, 0.6))  &
        plot_annotation(tag_levels = list(c("C", "D"))) &
        ggplot2::theme(
          plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
        )
    ) +
    plot_layout(ncol = 1, heights = c(0.55, 0.45))
)
suppl_fig_1

## post-task ratings -------------------------------------------------------------------

post_ratings_reff <- read_stand_data(pth_t1, "STAND-tasks-1-rew-eff-task-post-task-data-long") |>
  dplyr::mutate(timepoint = 0) |>
  dplyr::bind_rows(
    read_stand_data(pth_t2, "STAND-tasks-2-rew-eff-task-post-task-data-long") |>
      dplyr::mutate(timepoint = 1, intCond = NA)
  ) |>
  dplyr::group_by(subID) |>
  # for each subject, code intervention condition for timepoint 0 from timepoint 1
  dplyr::mutate(
    intCond = dplyr::if_else(
      timepoint == 1,
      dplyr::first(intCond[timepoint == 0]),
      intCond
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    question_type = dplyr::case_when(
      grepl("pleased", question, ignore.case = TRUE) ~ "pleased",
      grepl("achievement", question, ignore.case = TRUE) ~ "achievement",
      TRUE ~ "bored"
    ),
    intCond = factor(intCond, levels = c("CR", "BA")),
    block_no = as.numeric(stringr::str_replace(questionStage, "postBlock", "")),
  ) |>
  dplyr::select(1:2, intCond, timepoint, question_type, block_no) |>
  dplyr::mutate(rating = answer / 100)

# plot distributions of ratings by timepoint and condition
post_ratings_reff |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = rating, colour = interaction(intCond, timepoint), fill = interaction(intCond, timepoint)
    )
  ) +
  ggplot2::geom_density(alpha = 0.1, linewidth = 1) +
  ggplot2::facet_wrap(~question_type, scales = "free_y") +
  ggplot2::scale_fill_manual(
    values = c("CR.0" = "#ed968c", "BA.0" = "#88a0dc", "CR.1" = "#2c5d7f", "BA.1" = "#b8b74d"),
    labels = c("CR.0" = "CR baseline", "BA.0" = "BA baseline", "CR.1" = "CR follow-up", "BA.1" = "BA follow-up")
  ) +
  ggplot2::scale_color_manual(
    values = c("CR.0" = "#ed968c", "BA.0" = "#88a0dc", "CR.1" = "#2c5d7f", "BA.1" = "#b8b74d"),
    labels = c("CR.0" = "CR baseline", "BA.0" = "BA baseline", "CR.1" = "CR follow-up", "BA.1" = "BA follow-up")
  ) +
  ggplot2::labs(x = "Rating", y = "Density", fill = "", color = "") +
  cowplot::theme_half_open(font_family = "Open Sans") +
  ggplot2::theme(
    legend.position = "bottom",
    strip.background = ggplot2::element_rect(fill = "white"),
    strip.text = ggplot2::element_text(face = "plain", size = 14)
  ) +
  ggplot2::ggtitle(
    "reward-effort post-task self-reports"
  )
# realistically these are so skewed that modelling won't really be insightful
# but below is code to do so if desired with zero-one inflated beta regression

# # all three rating types
# results_pleased <- analyze_post_ratings(post_ratings_reff, "pleased")
# results_achievement <- analyze_post_ratings(post_ratings_reff, "achievement")
# results_bored <- analyze_post_ratings(post_ratings_reff, "bored")

# # view interaction effects
# results_pleased$conditional_effects_plot
# results_achievement$conditional_effects_plot
# results_bored$conditional_effects_plot

# # summarize interaction effects
# rating_pdir <- dplyr::bind_rows(
#   bayestestR::p_direction(results_pleased$interaction_effect$`b_intCondBA:timepoint`) |>
#     dplyr::mutate(rating = "pleased"),
#   bayestestR::p_direction(results_achievement$interaction_effect$`b_intCondBA:timepoint`) |>
#     dplyr::mutate(rating = "achievement"),
#   bayestestR::p_direction(results_bored$interaction_effect$`b_intCondBA:timepoint`) |>
#     dplyr::mutate(rating = "bored")
# )
# rating_hdis <- dplyr::bind_rows(
#   bayestestR::hdi(results_pleased$interaction_effect$`b_intCondBA:timepoint`, ci = 0.95) |>
#     dplyr::mutate(rating = "pleased"),
#   bayestestR::hdi(results_achievement$interaction_effect$`b_intCondBA:timepoint`, ci = 0.95) |>
#     dplyr::mutate(rating = "achievement"),
#   bayestestR::hdi(results_bored$interaction_effect$`b_intCondBA:timepoint`, ci = 0.95) |>
#     dplyr::mutate(rating = "bored")
# )

# # summarize inflation differences
# bayestestR::p_direction(results_pleased$one_inflation_diff$difference)
# bayestestR::p_direction(results_achievement$one_inflation_diff$difference)
# bayestestR::p_direction(results_bored$zero_inflation_diff$difference)