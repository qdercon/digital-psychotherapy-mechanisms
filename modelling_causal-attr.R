## MODELLING -------------------------------------------------------------------

library(patchwork)
pth_t1 <- "STAND-tasks-1/analysis/"
pth_t2 <- "STAND-tasks-2/analysis/"
dir.create("outputs", showWarnings = FALSE)
dir.create(paste0(pth_t2, "stan-fits"), showWarnings = FALSE)
source("model_fns.R")

## Load questionnaire data for exclusions
quest_t1 <- read_stand_data(pth_t1, "STAND-tasks-1-self-report-data")
quest_t2 <- read_stand_data(pth_t2, "STAND-tasks-2-self-report-data")
qqnrs_t2 <- read_stand_data(pth_t2, "qqnrs_t2")

# MAIN TEXT: causal attribution ================================================
data_long_ca <- read_stand_data(pth_t2, "stand-tasks-both-causal-attr-task-data-long")

## model-free analysis ---------------------------------------------------------
ca_behav_df <- data_long_ca |>
  dplyr::group_by(subID) |>
  dplyr::mutate(
    non_completed = !any(timepoint == 1),
    intCond = factor(intCond, levels = c("BA", "CR"))
  ) |>
  dplyr::ungroup()

ca_behav_plot_pos <- cattr_behav_plot(
  ca_behav_df,
  plot_type = "positive",
  mode = "randomised",
  plot_style = "slope",
  line_width = 1.5,
  err_width = 0.3,
  dodge_width = 0.45,
  point_size = 6,
  label_colour = "#b75347",
  ylim = c(0, 0.6),
  fnt_sz = 2,
  axis_fnt_sc = 11,
  brew_col = "Archambault",
  col_nums = c(4, 1)
)

ca_behav_plot_neg <- cattr_behav_plot(
  ca_behav_df,
  plot_type = "negative",
  mode = "randomised",
  plot_style = "slope",
  line_width = 1.5,
  err_width = 0.3,
  dodge_width = 0.45,
  point_size = 6,
  label_colour = "#2f70a1",
  ylim = c(0, 0.6),
  fnt_sz = 2,
  axis_fnt_sc = 11,
  brew_col = "Archambault",
  col_nums = c(4, 1)
) + ggplot2::theme(legend.position = "none")

figure_4ab <- wrap_elements(
  ca_behav_plot_pos + ca_behav_plot_neg +
    plot_layout(nrow = 2, heights = c(0.5, 0.5)) &
    plot_annotation(tag_levels = list(c("A", "B"))) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", face = "bold", size = 36)
    )
)

## linear mixed models of behavioural data
ca_pos_trials <- data_long_ca |> dplyr::filter(valence == "positive")
ca_neg_trials <- data_long_ca |> dplyr::filter(valence == "negative")

internalPos_mod <- brms::brm(
  internalChosen ~ intCond * timepoint + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = ca_pos_trials,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(internalPos_mod, rel = "CR")

globalPos_mod <- brms::brm(
  globalChosen ~ intCond * timepoint + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = ca_pos_trials,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(globalPos_mod, rel = "CR")

internalNeg_mod <- brms::brm(
  internalChosen ~ intCond * timepoint + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = ca_neg_trials,
  warmup = 2000,
  iter = 12000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(internalNeg_mod, rel = "CR")

globalNeg_mod <- brms::brm(
  globalChosen ~ intCond * timepoint + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = ca_neg_trials,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(globalNeg_mod, rel = "CR")

## model-based analyses --------------------------------------------------------------

stan_ls_ca <- prep_data(
  data_long_ca, type = "attribution", quest_t1 = quest_t1, quest_t2 = quest_t2,
  excl_quest = "standard"
)
ca_itt <- "causAttr-multisess-ITT-intervention-additive"
mod_ca_itt <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/caus-attr/", ca_itt, ".stan"))

## fit models with cmdstan
fit_ca_itt <- mod_ca_itt$sample(
  data = stan_ls_ca,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
fit_ca_itt$diagnostic_summary()
fit_ca_itt$summary(variables = "theta_int_internal_neg")
fit_ca_itt$save_output_files(
  dir = paste0(pth_t2, "stan-fits/"), timestamp = FALSE, random = FALSE
)

# extract and tidy draws
ca_vars <- c(
  "mu_internal_theta_neg[1]", "mu_internal_theta_neg[2]",
  "mu_internal_theta_pos[1]", "mu_internal_theta_pos[2]",
  "mu_global_theta_neg[1]", "mu_global_theta_neg[2]",
  "mu_global_theta_pos[1]", "mu_global_theta_pos[2]",
  "theta_int_internal_neg", "theta_int_internal_pos",
  "theta_int_global_neg", "theta_int_global_pos",
  "theta_ncompl_internal_neg", "theta_ncompl_internal_pos",
  "theta_ncompl_global_neg", "theta_ncompl_global_pos",
  "delta_internal_neg", "delta_internal_pos",
  "delta_global_neg", "delta_global_pos"
)
draws_ca_itt <- fit_ca_itt$draws(format = "df", variables = ca_vars)
saveRDS(draws_ca_itt, paste0(pth_t2, "stan-fits/", ca_itt, "_draws.rds"))
# draws_ca_itt <- readRDS(paste0(pth_t2, "stan-fits/", ca_itt, "_draws.rds"))
long_draws_ca <- draws_ca_itt |> tidybayes::gather_draws(!!!rlang::syms(ca_vars))

fit_ca_itt_tidy <- extract_draws(draws_ca_itt, long_draws_ca, hdi = 0.95)
print(fit_ca_itt_tidy$hdi, n = 24)
write.csv(fit_ca_itt_tidy$hdi, file = paste0("outputs/", ca_itt, "_hdi.csv"))

pd_in <- fit_ca_itt_tidy$df |> dplyr::filter(.variable == "theta_int_internal_neg") |> dplyr::pull(.value)
bayestestR::p_direction(pd_in, method = "direct")
pd_gp <- fit_ca_itt_tidy$df |> dplyr::filter(.variable == "theta_int_global_pos") |> dplyr::pull(.value)
bayestestR::p_direction(pd_gp, method = "direct")

# ncpl
pd_nc_in <- fit_ca_itt_tidy$df |> dplyr::filter(.variable == "theta_ncompl_internal_neg") |> dplyr::pull(.value)
pd_nc_ip <- fit_ca_itt_tidy$df |> dplyr::filter(.variable == "theta_ncompl_internal_pos") |> dplyr::pull(.value)
pd_nc_gn <- fit_ca_itt_tidy$df |> dplyr::filter(.variable == "theta_ncompl_global_neg") |> dplyr::pull(.value)
pd_nc_gp <- fit_ca_itt_tidy$df |> dplyr::filter(.variable == "theta_ncompl_global_pos") |> dplyr::pull(.value)
bayestestR::p_direction(pd_nc_in, method = "direct")
bayestestR::p_direction(pd_nc_ip, method = "direct")
bayestestR::p_direction(pd_nc_gn, method = "direct")
bayestestR::p_direction(pd_nc_gp, method = "direct")

## combined group-level results: absolute levels (top) + pre→post change (bottom)

# Associations between individuals' changes in parameters and changes in outcomes
ca_i_pars <- c(
  "theta_internal_neg", "theta_internal_pos", "theta_global_neg", "theta_global_pos",
  "delta_internal_neg_p", "delta_internal_pos_p", "delta_global_neg_p", "delta_global_pos_p"
)
ca_i_draws <- fit_ca_itt$draws(format = "df", variables = ca_i_pars)
saveRDS(ca_i_draws, paste0(pth_t2, "stan-fits/", ca_itt, "_indiv_draws.rds"))
# ca_i_draws <- readRDS(paste0(pth_t2, "stan-fits/", ca_itt, "_indiv_draws.rds"))
ca_hdi <- extract_indiv_pars(ca_i_draws, stan_ls_ca, type = "attribution")

## param: used for delta_{param} and theta_int_{param} column names
## mu_prefix: used for {mu_prefix}[1]/{mu_prefix}[2] group-level mean columns
##   (Stan's mu_* naming puts "internal"/"global" before "theta", unlike param's
##   "internal_pos"/"global_pos" ordering, so it can't be derived mechanically)
ca_params <- list(
  list(param = "internal_pos", mu_prefix = "mu_internal_theta_pos", nm = "θ internal-positive"),
  list(param = "global_pos",   mu_prefix = "mu_global_theta_pos",   nm = "θ global-positive"),
  list(param = "internal_neg", mu_prefix = "mu_internal_theta_neg", nm = "θ internal-negative"),
  list(param = "global_neg",   mu_prefix = "mu_global_theta_neg",   nm = "θ global-negative")
)
ca_change_tags <- c("C", "D", "E", "F")

ca_change_plts <- lapply(ca_params, function(p) {
  plot_ca_change(
    draws_ca_itt,
    param     = p$param,
    mu_prefix = p$mu_prefix,
    param_nm  = p$nm,
    diff_lab = "cogn. restr. vs.\n behav. activ.",
    col_nums    = c(1, 4),
    fnt_sz      = 2.3
  )
})
names(ca_change_plts) <- sapply(ca_params, `[[`, "param")

ca_change_panels <- mapply(function(prm, tag) {
  if (tag == "C" || tag == "E") {
    wrap_elements(
      ca_change_plts[[prm]]$change_panel &
        plot_annotation(tag_levels = list(tag)) &
        ggplot2::theme(
          plot.tag = ggplot2::element_text(family = "Open Sans", face = "bold", size = 36),
          plot.subtitle = ggplot2::element_blank()
        )
    )
  } else {
    wrap_elements(
      ca_change_plts[[prm]]$change_panel &
        plot_annotation(tag_levels = list(tag)) &
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          plot.subtitle = ggplot2::element_blank(),
          plot.tag     = ggplot2::element_text(
            family = "Open Sans", face = "bold", size = 36, hjust = 1.5
          )
        )
    )
  }
}, names(ca_change_plts), ca_change_tags, SIMPLIFY = FALSE)

figure_4cdef <- wrap_elements(
  (ca_change_panels[[1]] + ca_change_panels[[2]] + plot_layout(nrow = 1, widths = c(0.31, 0.23))) /
    (ca_change_panels[[3]] + ca_change_panels[[4]] + plot_layout(nrow = 1, widths = c(0.31, 0.23))) +
    plot_layout(nrow = 2)
)

# test the effect of intervention on individual-level parameters
# mod_ca_indiv <- lme4::lmer(
#   mean ~ group * t + (1 | id),
#   weights = 1 / hdi_width,
#   data = ca_hdi[ca_hdi$variable == "theta_global_pos" & ca_hdi$timepoint != "t2-t1", ]
# )
# library(lmerTest)
# broom.mixed::tidy(
#   as(mod_ca_indiv, "merModLmerTest"),
#   effects = "fixed",
#   conf.method = "boot",
#   conf.int = TRUE
# )

## Individual-level parameters

# int_neg <- plot_indiv_pars(
#   ca_hdi, "theta_internal_neg", "\u03B8 internal-negative", col_nums = c(4, 1),
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting")
# ) +
#   ggplot2::labs(
#     x = "mean (± s.d.) timepoint 1 (a.u.)", y = "mean (± s.d.) timepoint 2 (a.u.)",
#     title = "**internal-negative** attributions"
#   )
# int_pos <- plot_indiv_pars(
#   ca_hdi, "theta_internal_pos", "\u03B8 internal-positive", col_nums = c(4, 1),
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting")
# ) +
#   ggplot2::labs(
#     x = "mean (± s.d.) timepoint 1 (a.u.)", y = "mean (± s.d.) timepoint 2 (a.u.)",
#     title = "**internal-positive** attributions"
#   )
# glo_neg <- plot_indiv_pars(
#   ca_hdi, "theta_global_neg", "\u03B8 global-negative", col_nums = c(4, 1),
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting")
# ) +
#   ggplot2::labs(
#     x = "mean (± s.d.) timepoint 1 (a.u.)", y = "mean (± s.d.) timepoint 2 (a.u.)",
#     title = "**global-negative** attributions"
#   )
# glo_pos <- plot_indiv_pars(
#   ca_hdi, "theta_global_pos", "\u03B8 global-positive", col_nums = c(4, 1),
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting")
# ) +
#   ggplot2::labs(
#     x = "mean (± s.d.) timepoint 1 (a.u.)", y = "mean (± s.d.) timepoint 2 (a.u.)",
#     title = "**global-positive** attributions"
#   )

# indiv_ca_point_plt <- wrap_elements(
#   int_pos + int_neg + glo_pos + glo_neg +
#     plot_layout(nrow = 2, axis_title = "collect", guides = "collect") &
#     ggplot2::theme(
#       legend.position = "none",
#       legend.text = ggplot2::element_text(size = 18 * 1.2),
#       legend.key.spacing.x = grid::unit(18 * 1.2, "pt"),
#       axis.title.x = ggplot2::element_text(size = 20),
#       axis.title.y = ggplot2::element_text(size = 20)
#     )
# ) &
#   plot_annotation(tag_levels = list("B")) &
#   ggplot2::theme(
#     plot.tag.position = c(0.025, 0.975),
#     plot.tag = ggplot2::element_text(family = "Open Sans", size = 28, face = "bold")
#   )

# Boxplots of individual-level parameters
# int_pos_box <- plot_indiv_pars(
#   ca_hdi, "theta_internal_pos", "\u03B8 internal-positive", type = "box",
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#   compl = FALSE, indiv_lines = TRUE, fnt_sz = 1.2, col_nums = c(4, 1),
#   sz_scale = c(0.5, 2)
# )
# int_neg_box <- plot_indiv_pars(
#   ca_hdi, "theta_internal_neg", "\u03B8 internal-negative", type = "box",
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#   compl = FALSE, indiv_lines = TRUE, fnt_sz = 1.2, col_nums = c(4, 1),
#   sz_scale = c(0.5, 2)
# )
# glo_pos_box <- plot_indiv_pars(
#   ca_hdi, "theta_global_pos", "\u03B8 global-positive", type = "box",
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#   compl = FALSE, indiv_lines = TRUE, fnt_sz = 1.2, col_nums = c(4, 1),
#   sz_scale = c(0.5, 2)
# )
# glo_neg_box <- plot_indiv_pars(
#   ca_hdi, "theta_global_neg", "\u03B8 global-negative", type = "box",
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#   compl = FALSE, indiv_lines = TRUE, fnt_sz = 1.2, col_nums = c(4, 1),
#   sz_scale = c(0.5, 2)
# )

# ca_box <- wrap_elements(
#   int_pos_box + int_neg_box + glo_pos_box + glo_neg_box +
#     plot_layout(nrow = 2, guides = "collect") &
#     ggplot2::theme(legend.position = "bottom") &
#     plot_annotation(tag_levels = list("C")) &
#     ggplot2::theme(
#       plot.tag = ggplot2::element_text(family = "Open Sans", size = 28),
#       plot.tag.position = c(0, 1.14)
#     )
# )

# Effect of non-completion on timepoint 1 parameter estimates
# ca_ncpl <- plot_indiv_pars(
#   ca_hdi,
#   var = c(
#     "theta_internal_pos", "theta_internal_neg", "theta_global_pos", "theta_global_neg"
#   ),
#   var_nm = c(
#     "\u03B8 **internal-positive**", "\u03B8 **internal-negative**",
#     "\u03B8 **global-positive**", "\u03B8 **global-negative**"
#   ),
#   type = "box", tpt = "t1", sz_scale = c(1, 3),
#   clr = "completed", fnt_sz = 1.2, clr_levs = c(0, 1),
#   clr_labs = c("non-completers", "completers"),
#   compl = FALSE, col_nums = c(3, 6)
# )
# ca_ncpl_box <- wrap_elements(
#   ca_ncpl +
#     plot_layout(nrow = 1, guides = "collect") &
#     ggplot2::theme(legend.position = "bottom") &
#     plot_annotation(tag_levels = list("D")) &
#     ggplot2::theme(
#       plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
#     )
# )

## heatmap of correlations between change in parameters and change in symptoms
ca_par_labels <- c(
  "internal_pos" = "θ internal-<br>positive",
  "global_pos" = "θ global-<br>positive",
  "internal_neg" = "θ internal-<br>negative",
  "global_neg" = "θ global-<br>negative"
)

ca_qn_labels <- c(
  "PHQ9" = "ΔPHQ-9", "DAS" = "ΔDAS", "miniSPIN" = "ΔminiSPIN",
  "BADS" = "ΔBADS", "ERQCR" = "ΔERQ-CR", "AMI" = "ΔAMI-BA"
)

# raw per-timepoint draw columns are named theta_{param}[id,t], not
# {param}[id,t] -- map them onto the same keys used in ca_par_labels
ca_indiv_par_map <- c(
  "theta_internal_pos" = "internal_pos",
  "theta_global_pos" = "global_pos",
  "theta_internal_neg" = "internal_neg",
  "theta_global_neg" = "global_neg"
)

symptom_corrs_ca <- get_symptom_corrs(
  data_long = data_long_ca,
  model_type = "attribution",
  draws = ca_i_draws,
  stan_ls = stan_ls_ca,
  qqnrs = qqnrs_t2,
  delta_par_pattern = "delta_(.*)_p\\[(.*)\\]",
  par_labels = ca_par_labels,
  indiv_par_map = ca_indiv_par_map,
  qns_to_invert = c("ERQCR_delta", "BADS_delta")
)

ba_catr_heatmap <- plot_correlation_heatmap(
  cor_summary_df = symptom_corrs_ca$summary_df,
  int_group = "BA",
  par_labels = ca_par_labels,
  qn_labels = ca_qn_labels,
  cor_lims = c(-0.5, 0.5),
  htmp_title = "behavioural activation",
  title_clr = "#88a0dc",
  label_type = "r_only",
  fnt_sz = 1.4,
  txt_sz = 8.5
)
cr_catr_heatmap <- plot_correlation_heatmap(
  cor_summary_df = symptom_corrs_ca$summary_df,
  int_group = "CR",
  par_labels = ca_par_labels,
  qn_labels = ca_qn_labels,
  cor_lims = c(-0.5, 0.5),
  htmp_title = "cognitive restructuring",
  title_clr = "#ed968c",
  label_type = "r_only",
  fnt_sz = 1.4,
  txt_sz = 8.5
)

# pre-/post-intervention level correlations (parameter value vs. symptom score
# at the same timepoint), e.g.:
# ba_catr_heatmap_t1 <- plot_correlation_heatmap(
#   cor_summary_df = symptom_corrs_ca$summary_df,
#   int_group = "BA", timepoint = "t1",
#   par_labels = ca_par_labels, qn_labels = ca_qn_labels,
#   htmp_title = "behavioural activation (timepoint 1)",
#   title_clr = "#88a0dc", label_type = "r_only", fnt_sz = 1.4, txt_sz = 8.5
# )

figure_4g <- wrap_elements(
  cr_catr_heatmap + ba_catr_heatmap +
    plot_layout(nrow = 2, guides = "collect") &
    ggplot2::theme(legend.position = "right") &
    plot_annotation(tag_levels = list("G")) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(family = "Open Sans", face = "bold", size = 36)
    )
)

# Association between parameters and questionnaires
# ca_hdi_qn <- ca_hdi |>
#   dplyr::left_join(ids, by = c("id", "completed")) |>
#   dplyr::left_join(qqnrs_t2, by = c("group", "completed", "subID"))

# int_neg_phq9 <- plot_indiv_pars(
#   ca_hdi_qn, "internal_neg", "\u03B8 internal-negative", compl = TRUE,
#   qnr = "PHQ9_total_delta", qnr_nm = "PHQ-9", col_nums = c(4, 1),
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#   #cor_xlab = -10, cor_ylab = c(2.2, 1.8), sz_scale = c(1, 3)
# )
# int_pos_phq9 <- plot_indiv_pars(
#   ca_hdi_qn, "internal_pos", "\u03B8 internal-positive", compl = TRUE,
#   qnr = "PHQ9_total_delta", qnr_nm = "PHQ-9", col_nums = c(4, 1),
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#   cor_xlab = -10, cor_ylab = c(3.0, 2.5), sz_scale = c(1, 3)
# )
# glo_neg_phq9 <- plot_indiv_pars(
#   ca_hdi_qn, "global_neg", "\u03B8 global-negative", compl = TRUE,
#   qnr = "PHQ9_total_delta", qnr_nm = "PHQ-9", col_nums = c(4, 1),
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#   cor_xlab = -10, cor_ylab = c(2.2, 1.7), sz_scale = c(1, 3)
# )
# glo_pos_phq9 <- plot_indiv_pars(
#   ca_hdi_qn, "global_pos", "\u03B8 global-positive", compl = TRUE,
#   qnr = "PHQ9_total_delta", qnr_nm = "PHQ-9", col_nums = c(4, 1),
#   clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#   cor_xlab = -10, cor_ylab = c(2.2, 1.8), sz_scale = c(1, 3)
# )

# ca_phq9_cors <- wrap_elements(
#   int_pos_phq9 + int_neg_phq9 + glo_pos_phq9 + glo_neg_phq9 +
#     plot_layout(axis_title = "collect", guides = "collect") &
#     ggplot2::theme(legend.position = "none") &
#     plot_annotation(tag_levels = list("E")) &
#     ggplot2::theme(
#       plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
#     )
# )
# cors_ncpl_all <- wrap_elements(ca_ncpl_box + ca_phq9_cors)

# # theta internal negative predicts all questionnaires
# qrns <- colnames(qqnrs_t2)[grep("total_delta", colnames(qqnrs_t2))][-1]
# qnrnms <- c("miniSPIN", "DAS", "ERQ-CR", "BADS", "AMI-CA")
# int_neg_ls <- list()
# for (q in qrns) {
#   min_score <- min(qqnrs_t2[[q]], na.rm = TRUE)
#   int_neg_ls[[q]] <- plot_indiv_pars(
#     ca_hdi_qn, "internal_neg", "\u03B8 internal-negative", compl = TRUE,
#     qnr = q, qnr_nm = qnrnms[which(qrns == q)], col_nums = c(4, 1),
#     clr_levs = c("CR", "BA"), clr_labs = c("restructuring", "goal-setting"),
#     cor_xlab = min_score, cor_ylab = c(2.2, 1.8), sz_scale = c(1, 3)
#   )
# }

# int_neg_qrn <- wrap_elements(
#   int_neg_ls$miniSPIN_total_delta +
#     (int_neg_ls$DAS_total_delta + ggplot2::ylab("")) +
#     (int_neg_ls$ERQCR_total_delta + ggplot2::ylab("")) +
#     (int_neg_ls$BADS_total_delta + ggplot2::ylab("")) +
#     (int_neg_ls$AMI_total_delta + ggplot2::ylab("")) +
#     plot_layout(nrow = 1) &
#     ggplot2::theme(legend.position = "none") &
#     plot_annotation(tag_levels = list("F")) &
#     ggplot2::theme(
#       plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
#     )
# )

# pars_point_plts_ca + cors_ncpl_all + int_neg_qrn +
#   plot_layout(heights = c(0.55, 0.3, 0.15))

## mediation analysis: group * time * parameter on PHQ9 -------------------
# prepare data
## z_sd_ref = "t2" puts PHQ-9 in units of the COMPLETE-CASE sd(t2), i.e. the same b_std the
## questionnaire ANCOVA forest is drawn on (qnr_analysis_ancova.R), so this model's direct_arm and
## total_arm can be checked against it below. NB not a pure relabel -- the .stan priors are fixed, so
## the divisor changes how informative they are, and fits on the old "pooled" scaling are not
## comparable. See prep_data()'s qnr_df block.
stan_ls_symptoms_ca <- prep_data(
  data_long_ca, type = "attribution", qnr_df = qqnrs_t2,
  quest_t1 = quest_t1, quest_t2 = quest_t2, excl_quest = "standard", z_sd_ref = "t2"
)

## joint mediation model: all four attribution-style parameters adjusted for each other, with an
## ANCOVA symptom submodel. the previous model (...-joint-decomp, no suffix) put a participant random
## intercept on BOTH PHQ-9 rows, which imposes compound symmetry -- it forces Var(t1) == Var(t2), and
## PHQ-9 violates that here (sd 3.29 -> 3.86, variance ratio 1.37, because enrolment screened on
## baseline PHQ-9 but not on follow-up). The ANCOVA version conditions on the observed baseline
## instead. The ATTRIBUTION submodel is byte-identical between the two; see the new file's header for
## what changes in the reportable quantities (one direct_arm instead of a per-arm pair, one
## total_arm, no prop_med_*).
ca_joint_med_mod_name <- "causAttr-multisess-ITT-symptoms-joint-decomp-ancova"
mod_ca_joint_med <- cmdstanr::cmdstan_model(paste0(pth_t2, "/stan-models/caus-attr/", ca_joint_med_mod_name, ".stan"))

fit_ca_joint_med <- mod_ca_joint_med$sample(
  data = stan_ls_symptoms_ca,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 4000,
  refresh = 1000
)
fit_ca_joint_med$save_object(file = paste0(pth_t2, "stan-fits/", ca_joint_med_mod_name, "-mediation-16000.rds"))
# ca_joint_med_rds <- paste0(pth_t2, "stan-fits/", ca_joint_med_mod_name, "-mediation-16000.rds")
# if (!file.exists(ca_joint_med_rds)) {
#   stop(
#     "no cached fit at ", ca_joint_med_rds, ".\n",
#     "The symptom submodel was reformulated as an ANCOVA, so the old ...-joint-decomp fit is stale ",
#     "and is deliberately NOT reused (its draws have direct_ba/direct_cr, which the ANCOVA does not ",
#     "define). Uncomment the $sample() block above and re-run.",
#     call. = FALSE
#   )
# }
# fit_ca_joint_med <- readRDS(ca_joint_med_rds)

ca_med_vars <- c(
  "alpha_phq", "beta_t1_phq", "beta_phq_group", "sigma_phq",
  "a_int_neg_cr", "a_int_neg_ba", "a_int_pos_cr", "a_int_pos_ba",
  "a_glob_neg_cr", "a_glob_neg_ba", "a_glob_pos_cr", "a_glob_pos_ba",
  # observed-t2-only companions to the ITT a-paths: only 96 of 188 return, so the ITT a-paths lean
  # heavily on the model-imputed t2 thetas, and the gap between the two is worth reporting
  "a_int_neg_cr_obs", "a_int_neg_ba_obs", "a_int_pos_cr_obs", "a_int_pos_ba_obs",
  "a_glob_neg_cr_obs", "a_glob_neg_ba_obs", "a_glob_pos_cr_obs", "a_glob_pos_ba_obs",
  "b_int_neg_cr", "b_int_neg_ba", "b_int_pos_cr", "b_int_pos_ba",
  "b_glob_neg_cr", "b_glob_neg_ba", "b_glob_pos_cr", "b_glob_pos_ba",
  # ANCOVA: ONE direct effect and ONE total, not a per-arm pair -- the outcome is the t2 LEVEL, so
  # there is no within-arm change term to split per arm, and prop_med_* no longer exists
  "direct_arm", "total_arm",
  "indirect_int_neg_cr", "indirect_int_neg_ba", "indirect_int_pos_cr", "indirect_int_pos_ba",
  "indirect_glob_neg_cr", "indirect_glob_neg_ba", "indirect_glob_pos_cr", "indirect_glob_pos_ba",
  "index_mod_med_int_neg", "index_mod_med_int_pos", "index_mod_med_glob_neg", "index_mod_med_glob_pos"
)
# confirm the fit actually carries every quantity read below, BEFORE any of it is used
ca_symp_decomp_vars_names <- function() {
  c(
    paste0("beta_phq_", rep(c("int_neg", "int_pos", "glob_neg", "glob_pos"), each = 4),
           c("B", "B_group", "W", "W_group")),
    paste0("delta_", c("internal_neg", "internal_pos", "global_neg", "global_pos"), "_p")
  )
}
ca_missing_vars <- setdiff(
  sub("\\[.*", "", c(ca_med_vars, ca_symp_decomp_vars_names())),
  fit_ca_joint_med$metadata()$stan_variables
)
if (length(ca_missing_vars)) {
  stop("fit is missing generated quantities: ", paste(ca_missing_vars, collapse = ", "), call. = FALSE)
}
ca_symp_decomp_vars <- c(
  "beta_phq_int_negB", "beta_phq_int_negB_group", "beta_phq_int_negW", "beta_phq_int_negW_group",
  "beta_phq_int_posB", "beta_phq_int_posB_group", "beta_phq_int_posW", "beta_phq_int_posW_group",
  "beta_phq_glob_negB", "beta_phq_glob_negB_group", "beta_phq_glob_negW", "beta_phq_glob_negW_group",
  "beta_phq_glob_posB", "beta_phq_glob_posB_group", "beta_phq_glob_posW", "beta_phq_glob_posW_group",
  "delta_internal_neg_p", "delta_internal_pos_p", "delta_global_neg_p", "delta_global_pos_p"
)
ca_joint_med_vars <- c(ca_med_vars, ca_symp_decomp_vars)
# fit_ca_joint_med$summary(variables = ca_joint_med_vars)
draws_ca_joint_med <- fit_ca_joint_med$draws(format = "df", variables = ca_joint_med_vars)

## between-/within-person PHQ-9 <-> attribution-parameter associations
## from the prefix, so this mapping has to be exact.
ca_param_short_long <- c(
  int_neg = "internal_neg", int_pos = "internal_pos",
  glob_neg = "global_neg", glob_pos = "global_pos"
)
for (short in names(ca_param_short_long)) {
  long <- ca_param_short_long[[short]]
  draws_ca_joint_med[[paste0("phq_", long, "_t1_group_ba")]] <- draws_ca_joint_med[[paste0("beta_phq_", short, "B")]]
  draws_ca_joint_med[[paste0("phq_", long, "_t1_group_cr")]] <-
    draws_ca_joint_med[[paste0("beta_phq_", short, "B")]] + draws_ca_joint_med[[paste0("beta_phq_", short, "B_group")]]
  draws_ca_joint_med[[paste0("phq_", long, "_t1_diff")]]  <- draws_ca_joint_med[[paste0("beta_phq_", short, "B_group")]]
  draws_ca_joint_med[[paste0("phq_", long, "_t2_change_ba")]] <- draws_ca_joint_med[[paste0("beta_phq_", short, "W")]]
  draws_ca_joint_med[[paste0("phq_", long, "_t2_change_cr")]] <-
    draws_ca_joint_med[[paste0("beta_phq_", short, "W")]] + draws_ca_joint_med[[paste0("beta_phq_", short, "W_group")]]
  draws_ca_joint_med[[paste0("phq_", long, "_t2_change_diff")]] <-
    draws_ca_joint_med[[paste0("beta_phq_", short, "W_group")]]
}

# difference in b paths
get_hdi_pd(draws_ca_joint_med, "beta_phq_int_negW_group")
get_hdi_pd(draws_ca_joint_med, "beta_phq_int_posW_group")
get_hdi_pd(draws_ca_joint_med, "beta_phq_glob_negW_group")
get_hdi_pd(draws_ca_joint_med, "beta_phq_glob_posW_group")

## does symptom change still mediate for each attribution dimension, adjusting for the other three?
get_hdi_pd(draws_ca_joint_med, "indirect_int_pos_cr")
get_hdi_pd(draws_ca_joint_med, "indirect_int_pos_ba")
get_hdi_pd(draws_ca_joint_med, "index_mod_med_int_pos")

get_hdi_pd(draws_ca_joint_med, "indirect_glob_pos_cr")
get_hdi_pd(draws_ca_joint_med, "indirect_glob_pos_ba")
get_hdi_pd(draws_ca_joint_med, "index_mod_med_glob_pos")

## direct / indirect / total decomposition, with the model's own total_arm checked against the
## mediator-free questionnaire ANCOVA. condition == 1 = CR in the attribution models (the OPPOSITE
## of the reward-effort ones), so direct_arm and total_arm are CR - BA while the questionnaire
## forest is BA - CR; mediation_direct_total() flips the reference value to match rather than
## comparing them on trust.
ca_med_decomp <- mediation_direct_total(
  draws_ca_joint_med,
  suffixes = c("int_neg", "int_pos", "glob_neg", "glob_pos"),
  mediator_nms = c(
    int_neg = "theta internal-negative", int_pos = "theta internal-positive",
    glob_neg = "theta global-negative", glob_pos = "theta global-positive"
  ),
  arm_dir = "CR - BA",
  out_path = "outputs/mediation_direct_total_ca.csv"
)
print(ca_med_decomp, n = Inf)

## symptom decomposition panels (raw scatter + within-person association)
ca_decomp_param_nms <- c(
  internal_pos = "θ internal-positive",
  global_pos   = "θ global-positive",
  internal_neg = "θ internal-negative",
  global_neg   = "θ global-negative"
)
ca_symp_decomp_plts <- lapply(seq_along(ca_decomp_param_nms), function(i) {
  long <- names(ca_decomp_param_nms)[i]
  plot_sympt_decomp(
    draws_ca_joint_med,
    stan_ls  = stan_ls_symptoms_ca,
    prefix   = paste0("phq_", long),
    param_nm = ca_decomp_param_nms[[long]],
    type     = "attribution",
    diff_lab = "cogn. restr. vs.\nbehav. activ.",
    flip_order = FALSE,
    col_nums = c(1, 4), fnt_sz = 1.6
  )
})
names(ca_symp_decomp_plts) <- names(ca_decomp_param_nms)

ca_symp_change_tags <- c("H", "I", "J", "K")
ca_symp_change_panels <- Map(function(long, tag) {
  wrap_elements(
    wrap_elements(ca_symp_decomp_plts[[long]]$scatter_panel) +
      wrap_elements(ca_symp_decomp_plts[[long]]$within_panel) +
      plot_layout(widths = c(0.4, 0.6)) &
      plot_annotation(tag_levels = list(c(tag, ""))) &
      ggplot2::theme(
        plot.tag = ggplot2::element_text(family = "Open Sans", face = "bold", size = 36)
      )
  )
}, names(ca_decomp_param_nms), ca_symp_change_tags)

figure_4hijk <- wrap_elements(
  (ca_symp_change_panels[[1]] / ca_symp_change_panels[[2]] + plot_layout(ncol = 2))  /
    (ca_symp_change_panels[[3]] / ca_symp_change_panels[[4]] + plot_layout(ncol = 2)) +
    plot_layout(nrow = 2)
)

### mediation path diagrams
ca_mediator_cols <- stats::setNames(
  MetBrewer::met.brewer("Hokusai1")[c(3, 4, 6, 7)],
  c("int_neg", "glob_neg", "int_pos", "glob_pos")
)

## title colours match the valence colouring already used for the model-free plot
## (cattr_behav_plot()'s label_colour = c("#b75347", "#2f70a1") -- positive/negative)
ca_neg_med_plts <- plot_joint_mediation(
  draws_ca_joint_med,
  mediator1_suffix = "int_neg", mediator2_suffix = "glob_neg",
  mediator1_nm = "internal attribution", mediator2_nm = "global attribution",
  mediator_cols = ca_mediator_cols[c(1:2)],
  diff_lab = "cogn. restr. vs.\nbehav. activ.",
  title = "negative events", title_col = "#2f70a1",
  col_nums = c(1, 4), fnt_sz = 2.0
)
ca_pos_med_plts <- plot_joint_mediation(
  draws_ca_joint_med,
  mediator1_suffix = "int_pos", mediator2_suffix = "glob_pos",
  mediator1_nm = "internal attribution", mediator2_nm = "global attribution",
  mediator_cols = ca_mediator_cols[c(3:4)],
  diff_lab = "cogn. restr. vs.\nbehav. activ.",
  title = "positive events", title_col = "#b75347",
  brew_col = "Archambault",
  col_nums = c(1, 4), fnt_sz = 2.0
)

## inset the a1/b1, a2/b2 and residual (c') path plots onto the diagram, then stack the
## combined indirect-effect panel underneath -- same layout as figure_3g in
## modelling_reward-effort.R, factored out here since it's needed twice (neg/pos valence)
build_ca_med_figure <- function(plts, tag) {
  diagram <- wrap_elements(
    plts$diagram +
      patchwork::inset_element(
        plts$a1_path, left = 0.065, bottom = 0.65, right = 0.265, top = 0.975, align_to = "panel"
      ) +
      patchwork::inset_element(
        plts$b1_path, left = 0.735, bottom = 0.65, right = 0.935, top = 0.975, align_to = "panel"
      ) +
      patchwork::inset_element(
        plts$a2_path, left = 0.065, bottom = 0.025, right = 0.265, top = 0.35, align_to = "panel"
      ) +
      patchwork::inset_element(
        plts$b2_path, left = 0.735, bottom = 0.025, right = 0.935, top = 0.35, align_to = "panel"
      ) +
      patchwork::inset_element(
        plts$direct_path, left = 0.34, bottom = 0.43, right = 0.6, top = 0.72, align_to = "panel"
      ) +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  )
  wrap_elements(
    diagram + wrap_elements(plts$indirect) +
      plot_layout(widths = c(0.64, 0.36)) +
      plot_annotation(tag_levels = list(tag)) &
      ggplot2::theme(
        plot.tag = ggplot2::element_text(family = "Open Sans", face = "bold", size = 36)
      )
  )
}

figure_4l <- build_ca_med_figure(ca_pos_med_plts, "L")
figure_4m <- build_ca_med_figure(ca_neg_med_plts, "M")

## figure 4 all together --------------------------------------------------------

figure_4abcdef <- figure_4ab + figure_4cdef + plot_layout(nrow = 1, widths = c(0.44, 0.56))
figure_4ghijk <- wrap_elements(figure_4g + figure_4hijk + plot_layout(widths = c(0.35, 0.65)))
figure_4lm    <- figure_4l / figure_4m + plot_layout(nrow = 2)

fig_4 <- figure_4abcdef / figure_4ghijk / figure_4lm + plot_layout(nrow = 3, heights = c(0.28, 0.28, 0.44))
fig_4

## Simple baseline differences in non-completers -------------------------------

bsl_data <- prep_data(data_long_ca, type = "attribution", ret_df = TRUE) |>
  dplyr::filter(timepoint == 0) |>
  dplyr::mutate(non_completer = completed == 0)

## as above, fit each attribution style within its own valence
bsl_pos <- bsl_data |> dplyr::filter(valence == "positive")
bsl_neg <- bsl_data |> dplyr::filter(valence == "negative")

internalPos_ncpl <- brms::brm(
  internalChosen ~ non_completer + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = bsl_pos,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(internalPos_ncpl, rel = NULL, var = "non_completer")

globalPos_ncpl <- brms::brm(
  globalChosen ~ non_completer + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = bsl_pos,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(globalPos_ncpl, rel = NULL, var = "non_completer")

internalNeg_ncpl <- brms::brm(
  internalChosen ~ non_completer + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = bsl_neg,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(internalNeg_ncpl, rel = NULL, var = "non_completer")

globalNeg_ncpl <- brms::brm(
  globalChosen ~ non_completer + trialNo + (1 | subID),
  family = brms::bernoulli(link = "logit"),
  data = bsl_neg,
  warmup = 2000,
  iter = 12000,
  refresh = 1000,
  cores = 4,
  backend = "cmdstanr"
)
get_hdi_pd_brms(globalNeg_ncpl, rel = NULL, var = "non_completer")

## non-completer vs completer plots (CA) ----------------------------------------

mf_ncpl_ca <- cattr_behav_plot(
  ca_behav_df,
  plot_type = "both",
  mode = "completer_status",
  plot_style = "bar",
  err_width = 0.2,
  bar_width = 0.65,
  line_width = 1.5,
  ylim = c(0, 0.6),
  fnt_sz = 1.5,
  brew_col = "Archambault",
  col_nums = c(2, 5),
  label_colour = c("#b75347", "#2f70a1")
)

mf_ncpl_ca_plots <- wrap_elements(
  mf_ncpl_ca +
    plot_annotation(tag_levels = list("D")) &
    ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28)
    )
)

ca_ncpl_params <- list(
  list(
    prm = "theta_internal_pos", nm = "θ internal-positive",
    mu_col = "mu_internal_theta_pos[1]", ncpl_col = "theta_ncompl_internal_pos",
    suppress_labels = FALSE
  ),
  list(
    prm = "theta_global_pos", nm = "θ global-positive",
    mu_col = "mu_global_theta_pos[1]", ncpl_col = "theta_ncompl_global_pos",
    suppress_labels = TRUE
  ),
  list(
    prm = "theta_internal_neg", nm = "θ internal-negative",
    mu_col = "mu_internal_theta_neg[1]", ncpl_col = "theta_ncompl_internal_neg",
    suppress_labels = FALSE
  ),
  list(
    prm = "theta_global_neg", nm = "θ global-negative",
    mu_col = "mu_global_theta_neg[1]", ncpl_col = "theta_ncompl_global_neg",
    suppress_labels = TRUE
  )
)

ca_ncpl_plts <- lapply(ca_ncpl_params, function(p) {
  plot_ncpl_diff(
    draws_ca_itt,
    param    = p$prm,
    param_nm = p$nm,
    mu_col   = p$mu_col,
    ncpl_col = p$ncpl_col,
    suppress_labels = p$suppress_labels,
    diff_lab = "completers vs.\nnon-completers",
    col_nums = c(2, 5),
    fnt_sz   = 1.6
  ) +
    ggplot2::theme(plot.subtitle = ggplot2::element_blank())
})

ncpl_ca <- wrap_elements(
  (ca_ncpl_plts[[1]] + ca_ncpl_plts[[2]]) /
    (ca_ncpl_plts[[3]] + ca_ncpl_plts[[4]]) +
    plot_layout(nrow = 2) &
    plot_annotation(
      title = "session 1: non-completers vs. completers",
      tag_levels = list("E")
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = "Open Sans SemiBold", size = 20, colour = "grey30", margin = ggplot2::margin(l = 40, t = 10)
      ),
      plot.tag = ggplot2::element_text(family = "Open Sans", size = 28, vjust = -0.8)
    )
)

# Combined non-completer figure (rows: reff | ca)
# figure_5ab defined in modelling_reward-effort.R:
#   mf_ncpl_reff + ncpl_reff + plot_layout(nrow = 1, widths = c(0.4, 0.6))
figure_5de <- mf_ncpl_ca_plots + ncpl_ca + plot_layout(nrow = 1, widths = c(0.35, 0.65))
fig_5 <- figure_5abc / figure_5de + plot_layout(heights = c(0.67, 0.33))
fig_5
