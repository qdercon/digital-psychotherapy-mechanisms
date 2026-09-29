################################################################################
######## Figure 2 -- questionnaire analyses, baseline-adjusted (ANCOVA) #########
################################################################################

library(patchwork)
source("model_fns.R")

pth_t1 <- "STAND-tasks-1/analysis/"
pth_t2 <- "STAND-tasks-2/analysis/"
out_dir <- "outputs"
ancova_cache_dir <- paste0(pth_t2, "stan-fits/qnr")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ancova_cache_dir, showWarnings = FALSE, recursive = TRUE)

# arm direction: comp - ref = BA - CR, matching the earlier change-score analysis
# (compute_interaction(ref = "CR", comp = "BA")) and the sign of its forest panel
ARM_REF <- "CR"
ARM_INT <- "BA"

# fixed seed on every brm() below. without it these fits are not reproducible run-to-run and the
# lm/gls cross-checks become intermittent rather than deterministic.
ANCOVA_SEED <- 20260927

## 1. measure metadata --------------------------------------------------------------
# identical to the earlier change-score analysis, so the two cover exactly the same outcomes.
# the dagger marks measures where a HIGHER score is better (BADS, ERQ-CR); no sign flipping is done
# for the model panels -- the raw direction is kept and the caveat is carried as a note, exactly as
# the change-score figure does.
measure_meta <- tibble::tribble(
  ~measure,   ~prefix,    ~title,               ~subtitle,                                                  ~ylab,
  "PHQ-9",    "PHQ9",     "PHQ-9",              "patient health questionnaire",                             "PHQ-9 total (/27)", #nolint
  "DAS",      "DAS",      "DAS",                "dysfunctional attitudes scale",                            "DAS total (/36)", #nolint
  "miniSPIN", "miniSPIN", "miniSPIN",           "mini social phobia inventory",                             "miniSPIN total (/12)", #nolint
  "BADS",     "BADS",     "BADS<sup>†</sup>",   "behavioural activation for depression scale",              "BADS<sup>†</sup> total (/54)", #nolint
  "AMI",      "AMI",      "AMI-BA",             "apathy and motivation index\n(behavioural activation)",     "AMI-BA mean (/4)", #nolint
  "ERQ-CR",   "ERQCR",    "ERQ-CR<sup>†</sup>", "emotion regulation questionnaire\n(cognitive reappraisal)", "ERQ-CR<sup>†</sup> total (/36)" #nolint
)
main_measures <- measure_meta$measure
total_cols <- as.vector(outer(paste0(measure_meta$prefix, "_total_"), c("t1", "t2"), paste0))
prefix_to_measure <- stats::setNames(measure_meta$measure, measure_meta$prefix)

## 2. data --------------------------------------------------------------------------
# same rebuild as the earlier change-score analysis with ONE difference: no drop_na(PHQ9_total_delta).
# the cLDA in section 5c needs all 188 randomised, and the complete-case restriction is applied per
# measure by make_ancova_df() instead of globally here.
data_quest_wide_t2 <- read_stand_data(pth_t2, "STAND-tasks-2-self-report-data")
data_quest_wide_t1 <- read_stand_data(pth_t1, "STAND-tasks-1-self-report-data")
interventions <- read_stand_data(pth_t2, "stand-tasks-both-rew-eff-task-data-long") |>
  dplyr::distinct(subID, group)

qqnrs_t2 <- data_quest_wide_t2 |>
  dplyr::right_join(data_quest_wide_t1, by = "subID", suffix = c("_t2", "_t1")) |>
  dplyr::left_join(interventions, by = "subID") |>
  dplyr::rename(AMI_total_t1 = AMI_behavActiv_t1, AMI_total_t2 = AMI_behavActiv_t2) |>
  dplyr::select(subID, group, tidyselect::all_of(total_cols)) |>
  dplyr::mutate(
    dplyr::across(
      tidyselect::all_of(paste0(measure_meta$prefix, "_total_t2")),
      ~ . - get(sub("_t2$", "_t1", dplyr::cur_column())),
      .names = "{sub('_total_t2', '_total_delta', .col)}"
    ),
    group = factor(group, levels = c(ARM_REF, ARM_INT)),
    completed = !is.na(PHQ9_total_t2)
  )

stopifnot(
  nrow(qqnrs_t2) == 188,                       # all randomised
  sum(qqnrs_t2$completed) == 96,               # with a t2 PHQ-9
  !anyNA(qqnrs_t2$group),
  # the delta columns must reconstruct exactly -- guards the across() renaming above
  all(abs(qqnrs_t2$PHQ9_total_delta - (qqnrs_t2$PHQ9_total_t2 - qqnrs_t2$PHQ9_total_t1)) < 1e-9,
      na.rm = TRUE),
  all(abs(qqnrs_t2$ERQCR_total_delta - (qqnrs_t2$ERQCR_total_t2 - qqnrs_t2$ERQCR_total_t1)) < 1e-9,
      na.rm = TRUE)
)

# differential attrition is the specific threat to a complete-case analysis, so the retention rates
# are reported alongside the test rather than only the test
retention <- table(qqnrs_t2$group, qqnrs_t2$completed)
ret_test <- stats::chisq.test(retention)

## 3. shared aesthetics -------------------------------------------------------------
# unchanged from the earlier change-score figure -- same palette, shapes, fonts
grp_cols <- c(CR = "#ed968c", BA = "#88a0dc")[c(ARM_REF, ARM_INT)]
grp_shapes <- c(CR = 17, BA = 16)[c(ARM_REF, ARM_INT)]
grp_labs <- c(CR = "cognitive restructuring", BA = "behavioural activation")[c(ARM_REF, ARM_INT)]
# DISPLAY order for every panel: BA first, then CR. This is presentation only -- the MODEL keeps CR
# as its reference level (make_ancova_df() and make_clda_df() re-level to c(ref, comp) internally
# from ARM_REF/ARM_INT), so nothing here touches the sign of the contrast, which stays BA - CR.
# Colours, shapes and labels are matched BY NAME, so they follow their arm rather than a position.
arm_plot_levels <- c(ARM_INT, ARM_REF)
measure_cols <- stats::setNames(MetBrewer::met.brewer("Cassatt2", length(main_measures)), main_measures)
ancova_labels <- stats::setNames(measure_meta$title, measure_meta$measure)
# y axis is the measure itself, not a change score, so no delta prefix here
forest_labels <- c(
  "PHQ-9" = "PHQ-9", "DAS" = "DAS", "miniSPIN" = "miniSPIN",
  "BADS" = "BADS<sup>†</sup>", "AMI" = "AMI-BA", "ERQ-CR" = "ERQ-CR<sup>†</sup>"
)
tag_std <- ggplot2::element_text(family = "Open Sans", face = "bold", size = 28, colour = "black")

## 4. per-measure constants ---------------------------------------------------------
# computed ONCE on the primary frame and reused unchanged for every other frame, so all frames sit
# on one reference grid and their estimates are directly comparable
ancova_const <- stats::setNames(
  lapply(measure_meta$prefix, function(p) ancova_constants(qqnrs_t2, p)),
  measure_meta$measure
)
# sd_av must sit strictly between the two sds it averages, for every measure -- a cheap guard against
# a swapped or missing field silently rescaling panel A
stopifnot(vapply(ancova_const, function(c1) {
  c1$sd_av > min(c1$sd_t1, c1$sd_t2) && c1$sd_av < max(c1$sd_t1, c1$sd_t2) && c1$n == 96
}, logical(1)))

## 5a. primary ANCOVA: final score ~ arm + baseline ---------------------------------
fit_ancova <- function(df, pfx, frame_key, formula = t2 ~ group + t1_c, suffix = "") {
  brms::brm(
    formula = stats::as.formula(formula), data = df, family = stats::gaussian(),
    chains = 4, cores = 4, iter = 6000, warmup = 2000,
    backend = "cmdstanr", refresh = 0, seed = ANCOVA_SEED,
    file = file.path(ancova_cache_dir, paste0("ancova_", pfx, "_", frame_key, suffix)),
    file_refit = "on_change"
  )
}

# cross-check against the frequentist fit -- guards against a factor-level or centring slip. brms
# defaults leave class b FLAT, so there is no prior shrinkage to allow for and the two should agree
# to Monte Carlo error. the tolerance is ABSOLUTE in sd(t2) units rather than a ratio, since a ratio
# explodes on the near-null coefficients these analyses produce.
check_vs_lm <- function(v, df, formula, sd_t2, tol = 0.02) {
  b_lm <- stats::coef(stats::lm(stats::as.formula(formula), data = df))[[paste0("group", ARM_INT)]]
  stopifnot(abs(mean(v) - b_lm) <= tol * sd_t2)
  b_lm
}

ancova_fits <- list()
ancova_draws <- list()
ancova_summ <- list()

for (i in seq_len(nrow(measure_meta))) {
  meas <- measure_meta$measure[i]
  pfx <- measure_meta$prefix[i]
  cst <- ancova_const[[meas]]

  df <- make_ancova_df(qqnrs_t2, pfx, cst$t1_mean, ref = ARM_REF, comp = ARM_INT)
  fit <- fit_ancova(df, pfx, "primary")
  cd <- ancova_contrast_draws(fit, comp = ARM_INT)
  b_lm <- check_vs_lm(cd$contrast$.value, df, t2 ~ group + t1_c, cst$sd_t2)

  ancova_fits[[meas]] <- fit
  ancova_draws[[meas]] <- cd
  ancova_summ[[meas]] <- summarise_ancova_draws(fit, cd$contrast, cst, meas, "primary", nrow(df)) |>
    dplyr::mutate(b_lm = b_lm, .after = n)
}

ancova_tbl <- dplyr::bind_rows(ancova_summ) |>
  dplyr::mutate(measure = factor(measure, levels = main_measures)) |>
  dplyr::arrange(measure)
stopifnot(all(ancova_tbl$n == 96))
readr::write_csv(ancova_tbl, file.path(out_dir, "ancova_primary.csv"))
readr::write_csv(
  ancova_tbl |> dplyr::select(measure, frame, n, rhat_max, ess_bulk_min, ess_tail_min, divergent),
  file.path(out_dir, "ancova_diagnostics.csv")
)

# long draws for the forest, standardised by the primary-frame sd(t2) -- a fixed constant, so this
# is a pure relabelling of the posterior, not a refit
ancova_draws_long <- dplyr::bind_rows(lapply(main_measures, function(m) {
  tibble::tibble(
    measure = m, .draw = ancova_draws[[m]]$contrast$.draw,
    .value = ancova_draws[[m]]$contrast$.value / ancova_const[[m]]$sd_t2
  )
})) |>
  dplyr::mutate(measure = factor(measure, levels = rev(main_measures)))

## 5c. sensitivity: cLDA on all 188 randomised --------------------------------------
# the primary ANCOVA is necessarily complete-case (see the frame note at the top of this file), so
# 92 of the 188 randomised contribute nothing to it -- a far larger share than in the sibling study
# (69 of 361). this refits the arm effect in a model every randomised participant enters (everyone
# has a baseline score) to test whether the restriction changes the conclusion, rather than only
# acknowledging that it might.
#
# constrained longitudinal data analysis: score ~ post + post:group, no group main effect, so the
# arms are constrained to a shared baseline mean (justified by randomisation) and post:group{BA} is
# the arm difference in change.
#
# TWO THINGS THIS DOES NOT DO, both for the write-up:
#  1. it does not rescue a broken analysis. complete-case ANCOVA is ALREADY valid under MAR when
#     missingness depends on baseline, because it conditions on t1. agreement is the expected
#     result, not a lucky one.
#  2. with two timepoints, the 92 baseline-only participants inform the baseline mean and the
#     between-subject variance far more than the follow-up contrast, which is still carried by the
#     96. "includes all 188" is true but thin.
# neither model addresses MNAR.
clda_fits <- list()
clda_draws <- list()
clda_gls <- list()

for (i in seq_len(nrow(measure_meta))) {
  meas <- measure_meta$measure[i]
  pfx <- measure_meta$prefix[i]
  cst <- ancova_const[[meas]]
  dl <- make_clda_df(qqnrs_t2, pfx, ref = ARM_REF, comp = ARM_INT)

  # the constraint is ASSERTED, never assumed: the wrong (factor) formula fits happily and returns a
  # plausible number, so a silent regression here would be invisible downstream
  mm_c <- stats::model.matrix(~ post + post:group, data = dl)
  b_grp <- paste0("post:group", ARM_INT)
  stopifnot(
    !(paste0("group", ARM_INT) %in% colnames(mm_c)),  # no free baseline arm difference
    b_grp %in% colnames(mm_c),
    all(mm_c[dl$post == 0L, b_grp] == 0),             # group column is 0 at baseline
    dplyr::n_distinct(dl$subID) == nrow(qqnrs_t2),    # every randomised participant present
    sum(dl$post == 0L) == 188
  )

  # UNSTRUCTURED covariance + per-timepoint residual sd, NOT a random intercept. This is
  # load-bearing, not a refinement: `(1 | subID)` imposes compound symmetry, i.e. forces
  # Var(t1) == Var(t2), which is the very misspecification the mediation models are being
  # reformulated to remove (PHQ-9 goes sd 3.29 -> 3.86 here, a variance ratio of 1.37, because
  # enrolment screened on baseline PHQ-9 but not on follow-up). Using it here would leave the paper
  # internally inconsistent. The other five measures' variances barely move, so they agree either way.
  fit <- brms::brm(
    formula = brms::bf(score ~ post + post:group + unstr(time = tf, gr = subID), sigma ~ 0 + tf),
    data = dl, family = stats::gaussian(),
    chains = 4, cores = 4, iter = 6000, warmup = 2000,
    backend = "cmdstanr", refresh = 0, seed = ANCOVA_SEED,
    file = file.path(ancova_cache_dir, paste0("clda_", pfx)), file_refit = "on_change"
  )
  v <- coef_draws(fit, paste0("b_", b_grp))

  # frequentist cross-check on an independent implementation. gls with corSymm + varIdent is the
  # classic MMRM/cLDA fit and MATCHES THIS MODEL'S COVARIANCE -- lme4::lmer would not, being
  # random-intercept, so it would silently endorse the very specification the brm() above avoids.
  # tolerance is ABSOLUTE in sd(t2) units, not a ratio: these effects are near zero, where a
  # multiplicative tolerance explodes on a difference of no consequence.
  b_gls <- tryCatch(
    stats::coef(nlme::gls(
      score ~ post + post:group, data = dl,
      correlation = nlme::corSymm(form = ~ post + 1 | subID),
      weights = nlme::varIdent(form = ~ 1 | tf), method = "ML"
    ))[[b_grp]],
    error = function(e) {
      warning("cLDA gls cross-check failed to fit for ", meas, ": ", conditionMessage(e),
              " -- brms estimate is UNVERIFIED for this measure", call. = FALSE)
      NA_real_
    }
  )
  if (!is.na(b_gls)) {
    stopifnot(abs(mean(v) - b_gls) <= 0.05 * cst$sd_t2)
    # sign agreement only where the effect is big enough for a sign to be meaningful
    if (abs(b_gls) > 0.05 * cst$sd_t2) stopifnot(sign(mean(v)) == sign(b_gls))
  }

  clda_fits[[meas]] <- fit
  clda_draws[[meas]] <- v
  clda_gls[[meas]] <- b_gls
}

# side-by-side with the primary ANCOVA, standardised by the SAME sd(t2) constant so the two are
# directly comparable
clda_tbl <- dplyr::bind_rows(lapply(main_measures, function(meas) {
  cst <- ancova_const[[meas]]
  a <- ancova_summ[[meas]]
  dplyr::bind_rows(
    summarise_effect_draws(ancova_draws[[meas]]$contrast$.value, cst$sd_t2) |>
      dplyr::mutate(
        model = "ANCOVA (with outcome data)", n_ppts = 96L, b_freq = a$b_lm,
        rhat_max = a$rhat_max, ess_bulk_min = a$ess_bulk_min,
        ess_tail_min = a$ess_tail_min, divergent = a$divergent
      ),
    dplyr::bind_cols(
      summarise_effect_draws(clda_draws[[meas]], cst$sd_t2) |>
        dplyr::mutate(model = "cLDA (all randomised)", n_ppts = 188L, b_freq = clda_gls[[meas]]),
      brms_diag(clda_fits[[meas]])
    )
  ) |>
    dplyr::mutate(measure = meas, .before = 1)
}))
readr::write_csv(clda_tbl, file.path(out_dir, "ancova_clda_sensitivity.csv"))

## 6. Figure 2 panels ----------------------------------------------------------------

### panel A: per-participant standardised change, by measure and arm ----------------
# replaces the six per-measure pre->post slope panels (they move to the supplementary figure below):
# the same descriptive result on one effect-size axis, freeing the main figure for the model-based
# panels either side of it.
#
# UNITS: d_av -- (t2 - t1) / ((sd_t1 + sd_t2)/2) -- NOT the sd(t2) the model panels use, so the
# annotated numbers are a conventional effect size comparable with other literature. Two consequences
# for the caption:
#   1. panels A and B are on DIFFERENT standardisers. immaterial for five of six measures (their sds
#      barely move over two weeks, ratios 0.98-1.02) but PHQ-9's sd goes 3.29 -> 3.86 (ratio 1.17),
#      so pulling sd_t1 into the denominator shrinks it. that is an artefact of the study's own
#      screening criterion truncating baseline, not a discrepancy between the panels.
#   2. the denominator is computed ONCE per measure across both arms (see ancova_constants), so the
#      arms sit on a common scale -- a per-arm d_av would put them on different ones, and this panel
#      exists to be compared across arms by eye.
qqnrs_cc <- qqnrs_t2[qqnrs_t2$completed, ]

change_draws_long <- dplyr::bind_rows(lapply(seq_len(nrow(measure_meta)), function(i) {
  mm <- measure_meta[i, ]
  cst <- ancova_const[[mm$measure]]
  d <- make_ancova_df(qqnrs_t2, mm$prefix, cst$t1_mean, ref = ARM_REF, comp = ARM_INT)
  tibble::tibble(measure = mm$measure, group = d$group, .value = (d$t2 - d$t1) / cst$sd_av)
})) |>
  dplyr::mutate(
    measure = factor(measure, levels = main_measures),
    group = factor(as.character(group), levels = arm_plot_levels)
  )

# pooled (both arms) headline per measure: raw-score change AND d_av, so a reader can judge clinical
# magnitude on the measure's own scale alongside the effect size. deliberately NOT per arm -- that
# would double the 12 value labels already on this panel. d_av is computed exactly as in the earlier
# change-score figure, so the two figures cannot drift apart.
#
# NO paired t-test / significance star, unlike that figure: this figure set is presented entirely
# in Bayesian terms, so `panel_stats` carries no `p_star_timepoint` column and the slope panels in
# section 8 are drawn with show_stars = FALSE. (plot_qnr_measure() only reads that column inside its
# star branch, so its absence is safe -- but it is why the two must be changed together.)
compute_panel_stats <- function(meas, col_prefix) {
  t1 <- qqnrs_cc[[paste0(col_prefix, "_total_t1")]]
  t2 <- qqnrs_cc[[paste0(col_prefix, "_total_t2")]]
  keep <- !is.na(t1) & !is.na(t2)
  t1 <- t1[keep]
  t2 <- t2[keep]
  tibble::tibble(
    measure = meas,
    d_av = effectsize::repeated_measures_d(t2, t1, method = "av", adjust = FALSE)[["d_av"]],
    raw_change = mean(t2 - t1)
  )
}
panel_stats <- dplyr::bind_rows(
  lapply(seq_len(nrow(measure_meta)), function(i) {
    compute_panel_stats(measure_meta$measure[i], measure_meta$prefix[i])
  })
) |>
  dplyr::mutate(measure = factor(measure, levels = main_measures))

panel_raw_change <- panel_stats |>
  dplyr::transmute(
    measure, label = sprintf("raw &Delta; = %+.1f<br>*d*<sub>av</sub> = %+.2f", raw_change, d_av)
  )

# no dagger note on this panel: panel C carries the one explanation of it for the whole figure, and
# repeating it inside A only competes with the twelve value labels already there.
#
# the two label layers are sized and positioned together. `raw_change_nudge` is measured from the
# TOP of the per-arm value labels (plot_change_distributions() adds value_nudge before it), so it is
# the clearance between the two, and it has to grow with the fonts or the pooled raw/d_av block sits
# on top of the per-arm d it is meant to sit above.
change_plt <- plot_change_distributions(
  change_draws_long, labels = ancova_labels, base_size = 22, tag = "A",
  ylab = "pre → post change (*d*<sub>av</sub>)",
  grp_cols = grp_cols, grp_labs = grp_labs, grp_shapes = grp_shapes,
  show_values = TRUE, value_nudge = 0.5, value_size = 7, dodge_w = 0.8,
  raw_change = panel_raw_change, raw_change_nudge = 2, raw_change_size = 6,
  top_expand = 0.18
) +
  ggplot2::theme(plot.tag = tag_std)

### panel B: the pre-registered arm effect (forest) --------------------------------
# the group x questionnaire forest. axis is *b*<sub>std</sub>: the ANCOVA coefficient divided by that
# measure's primary-frame sd(t2). replaces the earlier frequentist change-score
# Cohen's d forest with significance stars.
forest_xlab <- paste0(
  "group difference in baseline-adjusted final score (*b*<sub>std</sub>)<br>",
  "<span style='font-family: \"Open Sans SemiBold\"; color:", grp_cols[[ARM_INT]], "'>",
  "(behavioural activation</span><span> − </span>",
  "<span style='font-family: \"Open Sans SemiBold\"; color:", grp_cols[[ARM_REF]], "'>",
  "cognitive restructuring)</span>"
)

forest_plt <- plot_ancova_forest(
  ancova_draws_long, labels = forest_labels, measure_cols = measure_cols,
  tag = "B", xlab = forest_xlab
) +
  ggplot2::theme(plot.tag = tag_std)

### panel C: correlations between change scores (heatmap) ------------------------
# same panel as the earlier change-score figure with the significance stars removed -- the cell
# label is now the correlation alone. computed on the 96 complete cases, and the two "higher = better"
# measures are sign-flipped so every delta points the same way.
pretty_nms <- list(
  "PHQ9_delta" = "ΔPHQ-9", "DAS_delta" = "ΔDAS", "miniSPIN_delta" = "ΔminiSPIN",
  "BADS_delta" = "ΔBADS<sup>†</sup>", "AMI_delta" = "ΔAMI-BA",
  "ERQCR_delta" = "ΔERQ-CR<sup>†</sup>"
)
stopifnot(identical(sub("_delta$", "", names(pretty_nms)), measure_meta$prefix))

cor_input_data <- qqnrs_cc |>
  dplyr::select(tidyselect::all_of(paste0(measure_meta$prefix, "_total_delta"))) |>
  stats::setNames(names(pretty_nms)) |>
  dplyr::mutate(ERQCR_delta = -ERQCR_delta, BADS_delta = -BADS_delta)

cor_results <- psych::corr.test(
  cor_input_data, method = "pearson", adjust = "none", ci = TRUE, minlength = 20
)
cor_df <- cor_results$ci |>
  tibble::as_tibble(rownames = "pair") |>
  tidyr::separate(col = "pair", into = c("qqnr1", "qqnr2"), sep = "-") |>
  dplyr::mutate(
    cor_coeff = cor_results$r[cbind(qqnr1, qqnr2)],
    label = sprintf("%.2f", cor_coeff),  # no significance stars
    qqnr1 = factor(qqnr1, levels = names(pretty_nms)),
    qqnr2 = factor(qqnr2, levels = rev(names(pretty_nms)))
  )

htmp <- cor_df |>
  ggplot2::ggplot(ggplot2::aes(x = qqnr1, y = qqnr2)) +
  ggplot2::geom_tile(ggplot2::aes(fill = cor_coeff)) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 4.8, colour = "white", family = "Open Sans") +
  ggtext::geom_richtext(
    data = data.frame(
      x = length(main_measures) - 1.9, y = length(main_measures) - 1.2,
      label = "<sup>†</sup> higher score = better;<br>signs flipped for correlation<br>analysis only"
    ),
    ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE, fill = NA,
    label.color = "slateblue4", label.padding = grid::unit(rep(0.4, 4), "lines"),
    label.r = grid::unit(0.5, "lines"), colour = "grey20", lineheight = 1.25,
    size = 5.2, family = "Open Sans"
  ) +
  ggplot2::scale_x_discrete(labels = unlist(pretty_nms)) +
  ggplot2::scale_y_discrete(labels = unlist(pretty_nms)) +
  ggplot2::scale_fill_gradientn(
    name = "*r*", colours = MetBrewer::met.brewer("Hokusai2", type = "continuous")
  ) +
  ggplot2::labs(x = "", y = "", tag = "C") +
  cowplot::theme_half_open(font_size = 18, font_family = "Open Sans") +
  ggplot2::theme(
    axis.text.y = ggtext::element_markdown(angle = 0, hjust = 1),
    axis.text.x = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(size = 14),
    plot.tag = tag_std
  )

## 7. Figure 2 assembly -------------------------------------------------------------
# A spans the full width (the descriptive change, all six measures); B and C sit beneath it.
# The cLDA was the other candidate for panel C and is NOT plotted: on the numbers (see
# outputs/ancova_clda_sensitivity.csv) it agrees with the complete-case ANCOVA to within 0.005 SD on
# every one of the six measures, which is reassurance that the complete-case restriction costs
# nothing rather than a finding in its own right. It stays in the analysis, the CSV and the write-up.

# printing is guarded only so that a graphics device which cannot resolve Open Sans (the default pdf
# device under a bare `Rscript` run, for instance) does not abort the script before the write-up
# section below. In an interactive session this behaves exactly like printing the object.
show_fig <- function(p) {
  tryCatch(print(p), error = function(e) {
    message("could not render on the current device (", conditionMessage(e),
            ") -- the object is still assigned")
  })
}
fig_2 <- wrap_elements(change_plt) /
  (wrap_elements(forest_plt) + wrap_elements(htmp) + plot_layout(widths = c(0.48, 0.52))) +
  plot_layout(heights = c(0.5, 0.5))
show_fig(fig_2)

## 8. Supplementary figure -----------------------------------------------------------

### the six per-measure pre->post slope panels ------------------------------------
# these used to be Figure 2A-F. NB plot_qnr_measure() reads `slope_summ`, `slope_long`,
# `panel_stats`, `grp_cols`, `grp_labs`, `grp_shapes` and `grp_lty` from the GLOBAL env rather than
# taking them as arguments (model_fns.R:314-317, 348-350), so all seven must exist before it is
# called. `grp_lty` is defined here rather than in section 3 to keep that dependency visible at the
# point of use.
# show_stars = FALSE: no significance bracket and no star, since nothing else in this figure set is
# presented frequentist-ly. What remains is the descriptive pooled d_av (both arms together), which
# is NOT the arm contrast -- that is Figure 2B, and a bracket spanning the two timepoints in a
# two-arm figure was easy to misread as one.
grp_lty <- c(CR = "32", BA = "solid")[c(ARM_REF, ARM_INT)]

slope_long <- qqnrs_cc |>
  dplyr::select(subID, group, tidyselect::all_of(total_cols)) |>
  tidyr::pivot_longer(
    cols = -c(subID, group), names_to = c("measure", "timepoint"),
    names_pattern = "(.*)_total_(t[12])"
  ) |>
  dplyr::mutate(
    measure = factor(prefix_to_measure[measure], levels = main_measures),
    timepoint = factor(timepoint, levels = c("t1", "t2"), labels = c("baseline", "post-intervention")),
    group = factor(as.character(group), levels = arm_plot_levels),
    subID = factor(subID)
  ) |>
  tidyr::drop_na(value)

slope_summ <- slope_long |>
  dplyr::group_by(measure, group, timepoint) |>
  dplyr::summarise(
    mean_score = mean(value), n = dplyr::n(),
    se = stats::sd(value) / sqrt(dplyr::n()),
    ci = stats::qt(0.975, dplyr::n() - 1) * (stats::sd(value) / sqrt(dplyr::n())),
    .groups = "drop"
  )

slope_panels <- lapply(seq_len(nrow(measure_meta)), function(i) {
  mm <- measure_meta[i, ]
  # this repo's plot_qnr_measure() themes plot.tag itself and takes no tag_theme argument, so the
  # shared tag style is appended as an override
  plot_qnr_measure(mm$measure, mm$title, mm$subtitle, mm$ylab, tag = LETTERS[i], dodge_w = 0.45,
                   show_stars = FALSE) +
    ggplot2::theme(plot.tag = tag_std)
})
suppl_slopes <- patchwork::wrap_plots(slope_panels, nrow = 2, guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
show_fig(suppl_slopes)

## 9. write-up statistics ------------------------------------------------------------
cat("\n========== QUESTIONNAIRE ANCOVA: t2 ~ arm + t1 (", ARM_INT, " − ", ARM_REF, ") ==========\n",
    sep = "")
cat(paste0(
  "the primary frame is the ", sum(qqnrs_t2$completed), " of ", nrow(qqnrs_t2), " randomised with BOTH\n",
  "timepoints. this is complete-case BY CONSTRUCTION (an ANCOVA whose outcome is the t2 score\n",
  "cannot use a participant without one), so it must be reported as 'with outcome data', never as\n",
  "intention-to-treat. section 5c refits every randomised participant into a cLDA.\n\n"
))

cat(sprintf(
  "retention: %s %d/%d (%.1f%%), %s %d/%d (%.1f%%); chi2(%d) = %.2f, p = %.3f\n\n",
  ARM_REF, retention[ARM_REF, "TRUE"], sum(retention[ARM_REF, ]),
  100 * retention[ARM_REF, "TRUE"] / sum(retention[ARM_REF, ]),
  ARM_INT, retention[ARM_INT, "TRUE"], sum(retention[ARM_INT, ]),
  100 * retention[ARM_INT, "TRUE"] / sum(retention[ARM_INT, ]),
  ret_test$parameter, ret_test$statistic, ret_test$p.value
))

cat("per-measure test-retest r (the quantity that makes ANCOVA worth the change from change scores:\n")
cat("change-score residual variance is 2(1-r), ANCOVA's is (1-r^2)):\n")
for (meas in main_measures) {
  cst <- ancova_const[[meas]]
  cat(sprintf(
    "  %-9s r = %.2f   sd_t1 = %5.2f  sd_t2 = %5.2f (ratio %.2f)   var ratio chg/ANCOVA = %.2f\n",
    meas, cst$r_t1t2, cst$sd_t1, cst$sd_t2, cst$sd_t2 / cst$sd_t1,
    2 * (1 - cst$r_t1t2) / (1 - cst$r_t1t2^2)
  ))
}

cat("\n==== primary ANCOVA (raw points; b/SD is b divided by that measure's sd(t2)) ====\n\n")
for (j in seq_len(nrow(ancova_tbl))) {
  r <- ancova_tbl[j, ]
  cat(sprintf(
    "%-9s n=%3d  b=%+.3f  95%%HDI=[%+.3f, %+.3f]  b/SD=%+.3f  b/sigma=%+.3f  pd=%.1f%%  (lm b=%+.3f)\n",
    r$measure, r$n, r$est, r$lo95, r$hi95, r$d_sd_t2, r$d_sigma, 100 * r$pd, r$b_lm
  ))
}

cat("\n==== diagnostics (worst across b_* and sigma per fit) ====\n\n")
print(
  ancova_tbl |>
    dplyr::select(measure, frame, rhat_max, ess_bulk_min, ess_tail_min, divergent) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 3))),
  n = Inf
)
if (any(ancova_tbl$rhat_max > 1.01) || any(ancova_tbl$divergent > 0)) {
  warning("ANCOVA fits: R-hat > 1.01 or divergent transitions present -- inspect before reporting")
} else {
  cat("\nall primary fits: R-hat <= 1.01, zero divergent transitions\n")
}

cat("\n\n========== MISSING DATA: complete-case ANCOVA vs cLDA on all randomised ==========\n")
cat(paste0(
  "complete-case ANCOVA is ALREADY MAR-valid given baseline, so agreement is the expected result,\n",
  "not a lucky one -- and neither model addresses MNAR. with two timepoints the 92 baseline-only\n",
  "participants inform the baseline mean and between-subject variance far more than the follow-up\n",
  "contrast, which is still carried by the 96.\n\n"
))
for (meas in main_measures) {
  rows <- clda_tbl |> dplyr::filter(measure == meas)
  for (j in seq_len(nrow(rows))) {
    r <- rows[j, ]
    cat(sprintf(
      "%-9s %-28s n=%3d  b=%+.3f  95%%HDI=[%+.3f, %+.3f]  b/SD=%+.3f  pd=%.1f%%  rhat=%.3f%s\n",
      if (j == 1) meas else "", r$model, r$n_ppts, r$est, r$lo95, r$hi95, r$d_sd_t2,
      100 * r$pd, r$rhat_max,
      if (is.na(r$b_freq)) "  [freq check UNAVAILABLE]" else sprintf("  (freq b=%+.3f)", r$b_freq)
    ))
  }
  d <- rows$d_sd_t2
  cat(sprintf("%-9s %-28s %+.3f SD\n\n", "", "-> cLDA - ANCOVA:", d[2] - d[1]))
}

cat("outputs written to ", out_dir, "/: ",
    paste(c("ancova_primary.csv", "ancova_diagnostics.csv", "ancova_clda_sensitivity.csv"),
          collapse = ", "), "\n", sep = "")


## Generate tables of demographic data -----------------------------------------

# helper functions

format_digits <- function(x, digit = 2, small = 2) {
  format(round(x, digits = digit), nsmall = small, scientific = FALSE)
}

mean_sd_range <- function(data, questionnaire, d = 2, s = 2) {
  mean_val <- mean(data[[questionnaire]], na.rm = TRUE)
  sd_val <- sd(data[[questionnaire]], na.rm = TRUE)
  range_val <- range(data[[questionnaire]], na.rm = TRUE)

  range_str <-
    if (is.integer(range_val[1]) && is.integer(range_val[2])) {
      paste0(range_val[1], "-", range_val[2])
    } else {
      paste0(
        format_digits(range_val[1], d, s), "-",
        format_digits(range_val[2], d, s)
      )
    }

  result <- paste0(
    format_digits(mean_val, d, s), " (",
    format_digits(sd_val, d, s), "; ",
    range_str, ")"
  )

  result
}

# 1. simple table comparing baseline scores by exclusion
# BEST (Kruschke, 2013, J Exp Psychol Gen 142:573) in place of the Student t-test. each group gets
# its own mean AND sd under a Student-t likelihood with a shared normality parameter nu, so it is
# robust to outliers and drops the equal-variance assumption the old t.test(var.equal = TRUE) made.
# the BEST package has left CRAN, so the model is rebuilt in brms with Kruschke's own broad priors,
# scaled to the pooled data (y = both groups together):
#   mu_g    ~ normal(mean(y), 1000 * sd(y))
#   sigma_g ~ uniform(sd(y) / 1000, sd(y) * 1000)  -- identity sigma link, so this is on sigma itself
#   nu      ~ exponential(1/29), bounded below at 1 by brms. by memorylessness that IS Kruschke's
#             nu - 1 ~ exponential(1/29), so no shift is needed.
# the table reports p_direction of mu_completed - mu_dropped where the p-value used to be.
qnr_names <- c(
  "PHQ-9" = "PHQ9_total", "DAS" = "DAS_total", "BADS" = "BADS_total",
  "ERQ-CR" = "ERQCR_total", "AMI-BA" = "AMI_behavActiv", "miniSPIN" = "miniSPIN_total"
)

complete <- data_quest_wide_t1 |>
  dplyr::filter(subID %in% data_quest_wide_t2$subID)

drop_out <- data_quest_wide_t1 |>
  dplyr::filter(!subID %in% data_quest_wide_t2$subID)

stopifnot(nrow(complete) == 96, nrow(drop_out) == 92, all(qnr_names %in% names(data_quest_wide_t1)))

fit_best <- function(col) {
  df <- data.frame(
    y = c(complete[[col]], drop_out[[col]]),
    g = factor(rep(c("completed", "dropped"), c(nrow(complete), nrow(drop_out))))
  )
  df <- df[!is.na(df$y), ]
  m <- mean(df$y)
  s <- stats::sd(df$y)
  # the sigma bounds go into the Stan code twice (parameter declaration AND uniform prior), so they
  # are formatted ONCE and reused: if the two disagree in the last digit, values at the edge of the
  # declared support fall outside the uniform's and the sampler cannot initialise
  lo <- format(signif(s / 1000, 4), scientific = FALSE)
  hi <- format(signif(s * 1000, 4), scientific = FALSE)
  brms::brm(
    formula = brms::bf(y ~ 0 + g, sigma ~ 0 + g),
    data = df, family = brms::student(link_sigma = "identity"),
    prior = c(
      brms::set_prior(sprintf("normal(%.4f, %.4f)", m, 1000 * s), class = "b"),
      brms::set_prior(
        sprintf("uniform(%s, %s)", lo, hi), class = "b", dpar = "sigma", lb = lo, ub = hi
      ),
      brms::set_prior("exponential(1.0 / 29)", class = "nu")
    ),
    chains = 4, cores = 4, iter = 6000, warmup = 2000,
    backend = "cmdstanr", refresh = 0, seed = ANCOVA_SEED,
    file = file.path(ancova_cache_dir, paste0("best_", col)), file_refit = "on_change"
  )
}

best_fits <- lapply(qnr_names, fit_best)
best_diag <- dplyr::bind_rows(lapply(best_fits, brms_diag, vars = c("b_", "nu")), .id = "measure")
if (any(best_diag$rhat_max > 1.01) || any(best_diag$divergent > 0)) {
  warning("BEST fits: R-hat > 1.01 or divergent transitions present -- inspect before reporting")
}

# completed - dropped, so pd is the posterior mass on the sign of the completer mean's deviation
best_diff <- lapply(best_fits, function(fit) {
  coef_draws(fit, "b_gcompleted") - coef_draws(fit, "b_gdropped")
})

bsl_compare <- tibble::tibble(
  "Questionnaire" = names(qnr_names),
  "Completed" = vapply(qnr_names, function(col) mean_sd_range(complete, col, 1, 1), character(1)),
  "Dropped out" = vapply(qnr_names, function(col) mean_sd_range(drop_out, col, 1, 1), character(1)),
  "P(direction) for difference" = vapply(best_diff, function(v) {
    format_digits(as.numeric(bayestestR::p_direction(v, method = "direct")), 3, 3)
  }, character(1))
)

# 2. sample characteristics by intervention in complete cases

# helper function

compl_t2 <- data_quest_wide_t2 |>
  dplyr::select(
    subID, tidyselect::starts_with("demogs"), tidyselect::contains("catch"), PHQ9_total
  ) |>
  dplyr::rename_all(~ stringr::str_replace(., "demogs_", "")) |>
  dplyr::mutate(
    group = ifelse(
      subID %in% interventions[interventions$group == "BA", ]$subID, "BA", "CR"
    )
  ) |>
  # get baseline PHQ-9 score from t1 data
  dplyr::left_join(
    data_quest_wide_t1 |>
      dplyr::select(subID, PHQ9_total) |>
      dplyr::rename(PHQ9_total_t1 = PHQ9_total),
    by = "subID"
  )

demog_rows <- c(
  "Cohort size",
  "Age, mean (SD; range)",
  "Gender, number (%)",
  "  Man",
  "  Woman",
  "Employment status, number (%)",
  "  Employed",
  "  Unemployed",
  "  Not seeking employment",
  "Financial status, number (%)",
  "  Doing okay",
  "  Getting by",
  "  Struggling",
  "Housing status, number (%)",
  "  Homeowner",
  "  Renting",
  "  Other",
  "Neurodivergent?, number (%)",
  "  Yes",
  "  No",
  "  Prefer not to say",
  "Disability that affects any of the below, number (%)",
  "  Concentrate for extended periods",
  "  Perform physically effortful activities",
  "  Read, write, or do maths",
  "  Deal with people you do not know",
  "  Other form of impact not listed",
  "  None of the above",
  "  Prefer not to say",
  "Previous treatment for mental ill-health, number (%)",
  "  Yes - talking therapy",
  "  Yes - medication",
  "  Yes - self-guided",
  "  Yes - other",
  "  No",
  "  Prefer not to say",
  "Baseline PHQ-9 score, mean (SD; range)"
)

compl_ba <- compl_t2[compl_t2$group == "BA", ]
n_ba <- nrow(compl_ba)
compl_cr <- compl_t2[compl_t2$group == "CR", ]
n_cr <- nrow(compl_cr)

demog_table <- tibble::tibble(
  "Characteristic" = demog_rows,
  "Behavioural activation" = c(
    n_ba,
    mean_sd_range(compl_ba, "age", 1, 1),
    NA,
    paste0(
      sum(compl_ba$gender == "man"), " (",
      format_digits(round(sum(compl_ba$gender == "man") * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(compl_ba$gender == "woman"), " (",
      format_digits(round(sum(compl_ba$gender == "woman") * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("full-time", compl_ba$employment)), " (",
      format_digits(round(sum(grepl("full-time", compl_ba$employment)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("seekers", compl_ba$employment)), " (",
      format_digits(round(sum(grepl("seekers", compl_ba$employment)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("not", compl_ba$employment)), " (",
      format_digits(round(sum(grepl("not", compl_ba$employment)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("okay", compl_ba$financial)), " (",
      format_digits(round(sum(grepl("okay", compl_ba$financial)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("getting", compl_ba$financial)), " (",
      format_digits(round(sum(grepl("getting", compl_ba$financial)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("struggling", compl_ba$financial)), " (",
      format_digits(round(sum(grepl("struggling", compl_ba$financial)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("homeowner", compl_ba$housing)), " (",
      format_digits(round(sum(grepl("homeowner", compl_ba$housing)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("tenant", compl_ba$housing)), " (",
      format_digits(round(sum(grepl("tenant", compl_ba$housing)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("other", compl_ba$housing)), " (",
      format_digits(round(sum(grepl("other", compl_ba$housing)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("yes", compl_ba$neurodiv)), " (",
      format_digits(round(sum(grepl("yes", compl_ba$neurodiv)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("no", compl_ba$neurodiv)), " (",
      format_digits(round(sum(grepl("no", compl_ba$neurodiv)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("prefer", compl_ba$neurodiv)), " (",
      format_digits(round(sum(grepl("prefer", compl_ba$neurodiv)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("concentrate", compl_ba$disability)), " (",
      format_digits(round(sum(grepl("concentrate", compl_ba$disability)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("perform", compl_ba$disability)), " (",
      format_digits(round(sum(grepl("perform", compl_ba$disability)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("read", compl_ba$disability)), " (",
      format_digits(round(sum(grepl("read", compl_ba$disability)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("deal", compl_ba$disability)), " (",
      format_digits(round(sum(grepl("deal", compl_ba$disability)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("other", compl_ba$disability)), " (",
      format_digits(round(sum(grepl("other", compl_ba$disability)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("none", compl_ba$disability)), " (",
      format_digits(round(sum(grepl("none", compl_ba$disability)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("prefer", compl_ba$disability)), " (",
      format_digits(round(sum(grepl("prefer", compl_ba$disability)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("talking", compl_ba$tx)), " (",
      format_digits(round(sum(grepl("talking", compl_ba$tx)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("medication", compl_ba$tx)), " (",
      format_digits(round(sum(grepl("medication", compl_ba$tx)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("self", compl_ba$tx)), " (",
      format_digits(round(sum(grepl("self", compl_ba$tx)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("other", compl_ba$tx)), " (",
      format_digits(round(sum(grepl("other", compl_ba$tx)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("no", compl_ba$tx)), " (",
      format_digits(round(sum(grepl("no", compl_ba$tx)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("prefer", compl_ba$tx)), " (",
      format_digits(round(sum(grepl("prefer", compl_ba$tx)) * 100 / n_ba, 3), 1, 1),
      ")"
    ),
    mean_sd_range(compl_ba, "PHQ9_total_t1", 1, 1)
  ),
  "Cognitive restructuring" = c(
    n_cr,
    mean_sd_range(compl_cr, "age", 1, 1),
    NA,
    paste0(
      sum(compl_cr$gender == "man"), " (",
      format_digits(round(sum(compl_cr$gender == "man") * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(compl_cr$gender == "woman"), " (",
      format_digits(round(sum(compl_cr$gender == "woman") * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("full-time", compl_cr$employment)), " (",
      format_digits(round(sum(grepl("full-time", compl_cr$employment)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("seekers", compl_cr$employment)), " (",
      format_digits(round(sum(grepl("seekers", compl_cr$employment)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("not", compl_cr$employment)), " (",
      format_digits(round(sum(grepl("not", compl_cr$employment)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("okay", compl_cr$financial)), " (",
      format_digits(round(sum(grepl("okay", compl_cr$financial)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("getting", compl_cr$financial)), " (",
      format_digits(round(sum(grepl("getting", compl_cr$financial)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("struggling", compl_cr$financial)), " (",
      format_digits(round(sum(grepl("struggling", compl_cr$financial)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("homeowner", compl_cr$housing)), " (",
      format_digits(round(sum(grepl("homeowner", compl_cr$housing)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("tenant", compl_cr$housing)), " (",
      format_digits(round(sum(grepl("tenant", compl_cr$housing)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("other", compl_cr$housing)), " (",
      format_digits(round(sum(grepl("other", compl_cr$housing)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("yes", compl_cr$neurodiv)), " (",
      format_digits(round(sum(grepl("yes", compl_cr$neurodiv)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("no", compl_cr$neurodiv)), " (",
      format_digits(round(sum(grepl("no", compl_cr$neurodiv)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("prefer", compl_cr$neurodiv)), " (",
      format_digits(round(sum(grepl("prefer", compl_cr$neurodiv)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("concentrate", compl_cr$disability)), " (",
      format_digits(round(sum(grepl("concentrate", compl_cr$disability)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("perform", compl_cr$disability)), " (",
      format_digits(round(sum(grepl("perform", compl_cr$disability)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("read", compl_cr$disability)), " (",
      format_digits(round(sum(grepl("read", compl_cr$disability)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("deal", compl_cr$disability)), " (",
      format_digits(round(sum(grepl("deal", compl_cr$disability)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("other", compl_cr$disability)), " (",
      format_digits(round(sum(grepl("other", compl_cr$disability)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("none", compl_cr$disability)), " (",
      format_digits(round(sum(grepl("none", compl_cr$disability)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("prefer", compl_cr$disability)), " (",
      format_digits(round(sum(grepl("prefer", compl_cr$disability)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    NA,
    paste0(
      sum(grepl("talking", compl_cr$tx)), " (",
      format_digits(round(sum(grepl("talking", compl_cr$tx)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("medication", compl_cr$tx)), " (",
      format_digits(round(sum(grepl("medication", compl_cr$tx)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("self", compl_cr$tx)), " (",
      format_digits(round(sum(grepl("self", compl_cr$tx)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("other", compl_cr$tx)), " (",
      format_digits(round(sum(grepl("other", compl_cr$tx)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("no", compl_cr$tx)), " (",
      format_digits(round(sum(grepl("no", compl_cr$tx)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    paste0(
      sum(grepl("prefer", compl_cr$tx)), " (",
      format_digits(round(sum(grepl("prefer", compl_cr$tx)) * 100 / n_cr, 3), 1, 1),
      ")"
    ),
    mean_sd_range(compl_cr, "PHQ9_total_t1", 1, 1)
  )
)

write.csv(bsl_compare, file = file.path(out_dir, "bsl_compare.csv"))
write.csv(demog_table, file = file.path(out_dir, "demog_table.csv"))
