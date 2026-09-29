## helper functions for plotting + modelling
`%>%` <- dplyr::`%>%`

# loads the public, id-scrambled copy of a data file (see scramble_ids.R), e.g.
# read_stand_data("STAND-tasks-2/analysis/", "qqnrs_t2"). scrambling moved subID to column 1 and
# kept write.csv()'s old row index as "...1", so the index is dropped by name -- a positional [-1]
# would silently drop subID instead.
read_stand_data <- function(dir, name) {
  df <- utils::read.csv(paste0(dir, name, "_scrambled.csv"))
  if (!"subID" %in% names(df)) stop("no subID column in ", name, "_scrambled.csv")
  df[setdiff(names(df), "...1")]
}

prep_data <- function(
    df,
    type,
    excl_eff_task = "none",
    excl_quest = "none",
    quest_t1 = NULL,
    quest_t2 = NULL,
    ret_df = FALSE,
    ret_ids = FALSE,
    cca = FALSE,
    qc = NULL,
    qnr_df = NULL,
    z_sd_ref = c("pooled", "t2"),
    recode_effort = FALSE) {
  # create session variable from timepoint
  df <- df |> dplyr::mutate(session = ifelse(timepoint == 0, 1, 2))

  # `z_sd_ref` picks how PHQ-9 is standardised when `qnr_df` is supplied (see the qnr_df block
  # at the end of this function). NB this is NOT a pure reparameterisation -- the symptom priors
  # in the .stan files are fixed, so changing the divisor changes how informative they are
  # relative to the data, and fits on the two settings are not comparable.
  z_sd_ref <- match.arg(z_sd_ref)
  excl_quest <- match.arg(excl_quest, choices = c("none", "standard", "strict"))
  excl_eff_task <- match.arg(excl_eff_task, choices = c("none", "standard", "strict"))

  excl_ids <- list()
  excl_ids$t1 <- NULL
  excl_ids$t2 <- NULL
  qn_excl_t1 <- 0
  qn_excl_t2 <- 0

  # define exclusion criteria for questionnaires
  if (excl_quest != "none") {
    if (is.null(quest_t1) || is.null(quest_t2)) {
      stop("Questionnaire data must be provided for questionnaire-based exclusion.")
    }
    # T1 exclusions
    quest_t1 <- quest_t1 |>
      dplyr::mutate(
        catchQ1 = ifelse(catch_1 >= 2, 0, 1),
        catchQ2 = ifelse(catch_2 == 0, 0, 1),
        fail_sum = catchQ1 + catchQ2
      )
    # T2 exclusions
    quest_t2 <- quest_t2 |>
      dplyr::mutate(
        catchQ1 = ifelse(catch_1 <= 1, 0, 1),
        catchQ2 = ifelse(catch_2 >= 2, 0, 1),
        fail_sum = catchQ1 + catchQ2
      )

    if (excl_quest == "standard") {
      excl_ids$t1 <- quest_t1 |>
        dplyr::filter(fail_sum == 2) |>
        dplyr::pull(subID)
      excl_ids$t2 <- quest_t2 |>
        dplyr::filter(fail_sum == 2) |>
        dplyr::pull(subID)
    } else if (excl_quest == "strict") {
      excl_ids$t1 <- quest_t1 |>
        dplyr::filter(fail_sum >= 1) |>
        dplyr::pull(subID)
      excl_ids$t2 <- quest_t2 |>
        dplyr::filter(fail_sum >= 1) |>
        dplyr::pull(subID)
    }
    qn_excl_t1 <- length(unique(excl_ids$t1))
    qn_excl_t2 <- length(unique(excl_ids$t2))
    both_excl <- intersect(excl_ids$t1, excl_ids$t2)

    # print how many excluded at each timepoint
    message(
      paste0(
        "Excluding ", qn_excl_t1, " participants at timepoint 1 and ",
        qn_excl_t2, " participants at timepoint 2, based on questionnaire catch trials."
      )
    )
  }

  if (type == "effort") {
    df <- df |>
      dplyr::mutate(
        catch_trial = ifelse(trialReward1 == 8 | trialReward2 == 8, 1, 0),
        catch_fail = ifelse(catch_trial == 1 & higherRewChosen == 0, 1, 0)
      ) |>
      dplyr::mutate(
        choice01 = dplyr::recode(choice, "route 1" = 0, "route 2" = 1),
        condition01 = dplyr::recode(group, "BA" = 1, "CR" = 0)
      )
    if (recode_effort) {
      df <- df |>
        dplyr::mutate(
          rew1_tmp = trialReward1,
          rew2_tmp = trialReward2,
          eff1_tmp = trialEffortPropMax1,
          eff2_tmp = trialEffortPropMax2,
          swap_effort = eff2_tmp < eff1_tmp,
          trialReward1 = dplyr::if_else(swap_effort, rew2_tmp, rew1_tmp),
          trialReward2 = dplyr::if_else(swap_effort, rew1_tmp, rew2_tmp),
          trialEffortPropMax1 = dplyr::if_else(swap_effort, eff2_tmp, eff1_tmp),
          trialEffortPropMax2 = dplyr::if_else(swap_effort, eff1_tmp, eff2_tmp),
          choice01 = dplyr::if_else(
            swap_effort,
            dplyr::if_else(choice01 == -1, -1, 1 - choice01),
            choice01
          )
        ) |>
        dplyr::select(-rew1_tmp, -rew2_tmp, -eff1_tmp, -eff2_tmp, -swap_effort)
    }

    if (excl_eff_task != "none") {
      # standard exclusion: remove participants who fail both catch trials
      # in EITHER session
      eff_excl_df <- df |>
        dplyr::group_by(subID, session) |>
        dplyr::mutate(catch_fail_no = sum(catch_fail, na.rm = TRUE)) |>
        dplyr::distinct(subID, session, catch_fail_no) |>
        dplyr::ungroup()

      if (excl_eff_task == "standard") eff_excl_df <- eff_excl_df |> dplyr::filter(catch_fail_no == 2)
      else if (excl_eff_task == "strict") eff_excl_df <- eff_excl_df |> dplyr::filter(catch_fail_no >= 1)
      excl_ids$t1 <- c(
        excl_ids$t1,
        eff_excl_df |>
          dplyr::filter(session == 1) |>
          dplyr::pull(subID)
      )
      excl_ids$t2 <- c(
        excl_ids$t2,
        eff_excl_df |>
          dplyr::filter(session == 2) |>
          dplyr::pull(subID)
      )
      both_excl <- intersect(excl_ids$t1, excl_ids$t2)

      # print how many excluded at each timepoint
      message(
        paste0(
          "Excluding ", length(unique(excl_ids$t1)) - qn_excl_t1, " participants at timepoint 1 and ",
          length(unique(excl_ids$t2)) - qn_excl_t2, " participants at timepoint 2 based on effort task catch trials.",
          ifelse(length(both_excl) > 0,
            paste0(" (", length(both_excl), " completers excluded at both timepoints)"),
            ""
          )
        )
      )
    }
  } else if (type == "attribution") {
    df <- df |>
      dplyr::mutate(
        neg_pos = ifelse(valence == "positive", 1, 0),
        condition01 = dplyr::recode(condition, "BA" = 0, "CR" = 1)
      )
  }

  # Exclude from dataframe
  df <- df |>
    dplyr::filter(
      !(session == 1 & subID %in% excl_ids$t1) & !(session == 2 & subID %in% excl_ids$t2)
    )

  # create numeric participant IDs
  df <- df |>
    dplyr::group_by(subID) |>
    dplyr::mutate(
      id = dplyr::cur_group_id(),
      completed = any(timepoint == 1),
      non_completed = !completed
    ) |>
    dplyr::ungroup()

  # get intervention condition for each participant
  all_ids <- df |>
    dplyr::distinct(id, subID, completed, condition01) |>
    dplyr::arrange(id)

  if (ret_ids) return(all_ids)

  if (cca) {
    df <- df[df$completed, ]
    all_ids <- all_ids[all_ids$completed, ]
  }

  if (ret_df) return(df)

  # get number of participants
  nPpts <- length(unique(df$id))
  # get number of trials for each participant
  # first, for session 1
  s1_trials <- df |>
    dplyr::filter(session == 1) |>
    dplyr::group_by(id) |>
    dplyr::summarise(n_trials = dplyr::n())
  # then for session 2
  s2_trials <- df |>
    dplyr::filter(session == 2) |>
    dplyr::group_by(id) |>
    dplyr::summarise(n_trials = dplyr::n())

  # get max number of trials
  nTrials_max <- max(c(s1_trials$n_trials, s2_trials$n_trials), na.rm = TRUE)

  # create arrays for number of trials
  nT_ppts <- array(0, dim = c(nPpts, 2))
  # fill with 0s for non-completers' session 2
  for (p in 1:nPpts) {
    ppt_id <- all_ids$id[p]
    s1_dat <- s1_trials[s1_trials$id == ppt_id, ]
    s2_dat <- s2_trials[s2_trials$id == ppt_id, ]
    nT_ppts[p, 1] <- ifelse(nrow(s1_dat) > 0, s1_dat$n_trials, 0)
    nT_ppts[p, 2] <- ifelse(nrow(s2_dat) > 0, s2_dat$n_trials, 0)
  }

  # create arrays for data
  if (type == "effort") {
    rew1 <- array(-1, dim = c(nPpts, 2, nTrials_max))
    rew2 <- array(-1, dim = c(nPpts, 2, nTrials_max))
    eff1 <- array(-1, dim = c(nPpts, 2, nTrials_max))
    eff2 <- array(-1, dim = c(nPpts, 2, nTrials_max))
    choice01 <- array(-1, dim = c(nPpts, 2, nTrials_max))
  } else if (type == "attribution") {
    internal_neg <- array(-1, dim = c(nPpts, 2, nTrials_max))
    internal_pos <- array(-1, dim = c(nPpts, 2, nTrials_max))
    global_neg <- array(-1, dim = c(nPpts, 2, nTrials_max))
    global_pos <- array(-1, dim = c(nPpts, 2, nTrials_max))
  }

  # loop over participants and sessions to fill arrays
  for (p in 1:nPpts) {
    ppt_id <- all_ids$id[p]
    for (t in 1:2) {
      ppt_sess_data <- df[df$id == ppt_id & df$session == t, ]
      if (nrow(ppt_sess_data) > 0) {
        if (type == "effort") {
          rew1[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$trialReward1
          rew2[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$trialReward2
          eff1[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$trialEffortPropMax1
          eff2[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$trialEffortPropMax2
          choice01[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$choice01
        } else if (type == "attribution") {
          internal_neg[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$internalChosen[ppt_sess_data$neg_pos == 0]
          internal_pos[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$internalChosen[ppt_sess_data$neg_pos == 1]
          global_neg[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$globalChosen[ppt_sess_data$neg_pos == 0]
          global_pos[p, t, 1:nT_ppts[p, t]] <- ppt_sess_data$globalChosen[ppt_sess_data$neg_pos == 1]
        }
      }
    }
  }

  # create list for stan
  stan_ls <- list(
    nTimes = 2,
    nPpts = nPpts,
    condition = all_ids$condition01,
    completed = all_ids$completed |> as.integer(),
    nTrials_max = nTrials_max,
    nT_ppts = nT_ppts
  )

  if (type == "effort") {
    stan_ls$rew1 <- rew1
    stan_ls$eff1 <- eff1
    stan_ls$rew2 <- rew2
    stan_ls$eff2 <- eff2
    stan_ls$choice01 <- choice01
  } else if (type == "attribution") {
    stan_ls$internal_neg <- internal_neg
    stan_ls$internal_pos <- internal_pos
    stan_ls$global_neg <- global_neg
    stan_ls$global_pos <- global_pos
  }

  if (!is.null(qc)) {
    # join qc data to all_ids
    qc_dat <- all_ids |>
      dplyr::left_join(qc, by = "subID") |>
      dplyr::select(subID, qc_score) |>
      # replace NAs with 0
      tidyr::replace_na(list(qc_score = 0))
    stan_ls$quality <- qc_dat$qc_score
  }

  if (!is.null(qnr_df)) {
    # join questionnaire data to all_ids
    qnr_dat <- all_ids |>
      dplyr::left_join(qnr_df, by = "subID") |>
      dplyr::select(subID, PHQ9_total_t1, PHQ9_total_t2)
    # Standardise PHQ-9 so the regression is on ~unit scale and default N(0,1) coefficient priors
    # are weakly informative rather than shrinking associations toward 0. Missing entries are
    # flagged with -999 *after* standardising. `z_sd_ref` picks the reference:
    #   "pooled" -- mean and sd over ALL observed person-timepoints. Historical behaviour, kept as
    #               the default so sourcing this file does not silently invalidate cached fits of
    #               the compound-symmetry symptom models, which were fitted on this scaling.
    #   "t2"     -- centre on the COMPLETE-CASE t1 mean, divide by the COMPLETE-CASE sd(t2). This is
    #               the scaling the ANCOVA symptom submodels expect: it puts every symptom
    #               coefficient in units of sd(t2), i.e. the same b_std the questionnaire ANCOVA
    #               forest is drawn on, so direct_arm / total_arm can be checked against it.
    #               sd(t1) would be the wrong divisor -- baseline is range-restricted by the study's
    #               own screening criterion (PHQ-9 sd 3.29 at t1 vs 3.86 at t2) -- and a pooled sd
    #               is not a quantity the questionnaire analysis ever forms.
    # Both timepoints are shifted and scaled by the SAME constants either way, so beta_t1_phq stays
    # a dimensionless same-instrument slope and only the intercept moves between the two settings.
    phq_cc <- !is.na(qnr_dat$PHQ9_total_t1) & !is.na(qnr_dat$PHQ9_total_t2)
    if (z_sd_ref == "pooled") {
      phq_vals <- c(qnr_dat$PHQ9_total_t1, qnr_dat$PHQ9_total_t2)
      phq_mean <- mean(phq_vals, na.rm = TRUE)
      phq_sd   <- stats::sd(phq_vals, na.rm = TRUE)
    } else {
      if (sum(phq_cc) < 2) stop("prep_data(): z_sd_ref = 't2' needs at least 2 complete PHQ-9 cases.")
      phq_mean <- mean(qnr_dat$PHQ9_total_t1[phq_cc])
      phq_sd   <- stats::sd(qnr_dat$PHQ9_total_t2[phq_cc])
      stopifnot(is.finite(phq_sd), phq_sd > 0)
    }
    qnr_dat <- qnr_dat |>
      dplyr::mutate(
        PHQ9_total_t1 = ifelse(is.na(PHQ9_total_t1), -999, (PHQ9_total_t1 - phq_mean) / phq_sd),
        PHQ9_total_t2 = ifelse(is.na(PHQ9_total_t2), -999, (PHQ9_total_t2 - phq_mean) / phq_sd)
      )
    stan_ls$PHQ9 <- as.matrix(qnr_dat[, -1])
    # keep the standardisation constants so coefficients can be back-transformed
    # to raw PHQ-9 points (multiply by phq_sd) if desired
    stan_ls$PHQ9_mean <- phq_mean
    stan_ls$PHQ9_sd   <- phq_sd
    # carried so a caller can assert the divisor matches ancova_constants()$sd_t2 rather than
    # trusting that two separately-computed constants happen to agree.
    # NB the scaling LABEL goes on as an attribute, not as a list element: cmdstanr rejects a
    # character entry in the data list outright ("Variable 'x' is of invalid type"), even for a
    # variable the model never declares. PHQ9_n_complete is numeric, so it travels as data like
    # PHQ9_mean/PHQ9_sd already do and is simply ignored by the model.
    attr(stan_ls, "PHQ9_z_sd_ref") <- z_sd_ref
    stan_ls$PHQ9_n_complete <- sum(phq_cc)  # NB captured BEFORE the -999 substitution above
  }

  stan_ls
}

# plot function for questionnaire measures
plot_qnr_measure <- function(meas, title, subtitle, ylab, tag, dodge_w,
                             err            = c("ci", "se"),
                             show_points    = TRUE,
                             annotate_stats = TRUE,
                             # the pooled pre-post paired-test bracket + star. TRUE keeps the
                             # existing behaviour (the earlier frequentist change-score figure);
                             # FALSE drops the bracket and star while keeping the d_av label, for a
                             # figure presented entirely in Bayesian terms. When FALSE, `panel_stats`
                             # need not carry a `p_star_timepoint` column at all -- it is only ever
                             # read inside the branch below.
                             show_stars     = TRUE,
                             base_size      = 18) {
  err <- match.arg(err)
  sm  <- slope_summ[slope_summ$measure == meas, ]
  sm$err <- if (err == "ci") sm$ci else sm$se
  iv  <- slope_long[slope_long$measure == meas, ]
  st  <- panel_stats[panel_stats$measure == meas, ]

  pd <- ggplot2::position_dodge(width = dodge_w)

  p <- ggplot2::ggplot(
    sm,
    ggplot2::aes(x = timepoint, y = mean_score, colour = group, group = group)
  )

  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = iv,
      ggplot2::aes(x = timepoint, y = value, colour = group, shape = group),
      position = ggbeeswarm::position_quasirandom(
        width = 0.12, dodge.width = dodge_w
      ),
      alpha = 0.34, size = 1.5, show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::geom_line(
      ggplot2::aes(linetype = group), position = pd, linewidth = 1.1
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_score - err, ymax = mean_score + err),
      position = pd, width = 0.12, linewidth = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = group), position = pd, size = 3
    ) +
    ggplot2::scale_colour_manual(values = grp_cols, labels = grp_labs, name = NULL) +
    ggplot2::scale_shape_manual(values = grp_shapes, labels = grp_labs, name = NULL) +
    ggplot2::scale_linetype_manual(values = grp_lty, labels = grp_labs, name = NULL) +
    ggplot2::labs(x = NULL, y = ylab, title = title, subtitle = subtitle, tag = tag) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = base_size) +
    cowplot::background_grid(major = "y", minor = "none") +
    ggplot2::theme(
      # title/axis.title.y rendered as markdown -- harmless for measures with
      # plain-text titles, and lets BADS/ERQ-CR render their superscript
      # dagger (set in measure_meta) via <sup>...</sup>.
      plot.title    = ggtext::element_markdown(size = base_size + 3, hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size - 5, hjust = 0.5, colour = "grey30"),
      plot.tag        = ggplot2::element_text(face = "bold", size = base_size + 7, colour = "black"),
      axis.title.y    = ggtext::element_markdown(size = base_size - 4),
      legend.position = "bottom",
      legend.justification.bottom = "center",
      legend.background = ggplot2::element_rect(
        fill = "white", colour = "#cacaca", linewidth = 0.5
      ),
      legend.margin = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm")
    )

  if (annotate_stats) {
    p <- p +
      ggtext::geom_richtext(
        data = data.frame(
          x = 1.5, y = Inf,
          label = paste0("*d*<sub>av</sub> = ", formatC(st$d_av, format = "f", digits = 2))
        ),
        ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
        vjust = 1.3, family = "Open Sans", size = (base_size - 6) / 2,
        fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt")
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.22)))

    if (show_stars) {
      data_range <- range(
        c(sm$mean_score - sm$err, sm$mean_score + sm$err, iv$value), na.rm = TRUE
      )
      bracket_y <- max(sm$mean_score + sm$err) + diff(data_range) * 0.10

      p <- p +
        ggplot2::annotate(
          "segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
          linewidth = 0.7, color = "grey70"
        ) +
        ggplot2::annotate(
          "text", x = 1.5, y = bracket_y, label = st$p_star_timepoint,
          vjust = -0.35, size = 4.3, family = "Open Sans", colour = "grey40"
        )
    }
  }
  p
}

# plot function for per-measure pre-post change (delta), split by group --
# companion to plot_qnr_measure(), for the supplementary subgroup-change figure.
# annotates the group x time interaction p-value (via compute_interaction()) and,
# per group, a one-sample d_s (mean delta / sd delta) labelling the within-group
# standardised change.
plot_qnr_change <- function(meas, title, subtitle, ylab, tag,
                            err            = c("ci", "se"),
                            show_points    = TRUE,
                            annotate_stats = TRUE,
                            base_size      = 18) {
  err <- match.arg(err)
  col <- delta_cols[[meas]]
  dd  <- qqnrs_t2[, c("group", col)]
  names(dd) <- c("group", "delta")
  dd  <- dd[!is.na(dd$delta), ]

  dsumm <- dd |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      mean_delta = mean(delta),
      n          = dplyr::n(),
      sd_delta   = stats::sd(delta),
      se         = sd_delta / sqrt(n),
      ci         = stats::qt(0.975, n - 1) * se,
      d_s        = mean_delta / sd_delta,
      .groups    = "drop"
    )
  dsumm$err <- if (err == "ci") dsumm$ci else dsumm$se

  intr <- compute_interaction(meas)

  p <- ggplot2::ggplot(dsumm, ggplot2::aes(x = group, y = mean_delta, colour = group))

  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = dd,
      ggplot2::aes(x = group, y = delta, colour = group, shape = group),
      position = ggbeeswarm::position_quasirandom(width = 0.25),
      alpha = 0.2, size = 1.5, show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "32", colour = "grey60") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_delta - err, ymax = mean_delta + err),
      width = 0.12, linewidth = 0.9
    ) +
    ggplot2::geom_point(ggplot2::aes(shape = group), size = 3) +
    ggplot2::scale_colour_manual(values = grp_cols, labels = grp_labs, name = NULL) +
    ggplot2::scale_shape_manual(values = grp_shapes, labels = grp_labs, name = NULL) +
    ggplot2::scale_x_discrete(labels = grp_labs_br) +
    ggplot2::labs(x = NULL, y = ylab, title = title, subtitle = subtitle, tag = tag) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = base_size) +
    cowplot::background_grid(major = "y", minor = "none") +
    ggplot2::theme(
      plot.title      = ggtext::element_markdown(size = base_size + 3, hjust = 0.5, face = "bold"),
      plot.subtitle   = ggplot2::element_text(size = base_size - 5, hjust = 0.5, colour = "grey30"),
      plot.tag        = ggplot2::element_text(face = "bold", size = base_size + 7, colour = "black"),
      axis.title.y    = ggtext::element_markdown(size = base_size - 4),
      legend.position = "bottom",
      legend.justification.bottom = "center",
      legend.background = ggplot2::element_rect(
        fill = "white", colour = "#cacaca", linewidth = 0.5
      ),
      legend.margin = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm")
    )

  if (annotate_stats) {
    data_range <- range(
      c(dsumm$mean_delta - dsumm$err, dsumm$mean_delta + dsumm$err, dd$delta), na.rm = TRUE
    )
    bracket_y <- max(dsumm$mean_delta + dsumm$err) + diff(data_range) * 0.10

    d_lab <- paste0(
      "<span style='color:", grp_cols[["BA"]], "'>*d*<sub>s</sub> = ",
      formatC(dsumm$d_s[dsumm$group == "BA"], format = "f", digits = 2), "</span> &middot; ",
      "<span style='color:", grp_cols[["CR"]], "'>*d*<sub>s</sub> = ",
      formatC(dsumm$d_s[dsumm$group == "CR"], format = "f", digits = 2), "</span>"
    )

    p <- p +
      ggtext::geom_richtext(
        data = data.frame(x = 1.5, y = Inf, label = d_lab),
        ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
        vjust = 1.3, family = "Open Sans", size = (base_size - 6) / 2,
        fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt")
      ) +
      ggplot2::annotate(
        "segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.7, color = "grey70"
      ) +
      ggplot2::annotate(
        "text", x = 1.5, y = bracket_y, label = intr$star,
        vjust = -0.35, size = 4.3, family = "Open Sans", colour = "grey40"
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.22)))
  }
  p
}

# Truncated normal sampler (vectorized)
rtruncnorm <- function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  if (sd <= 0) stop("sd must be > 0")
  p_low <- stats::pnorm(lower, mean = mean, sd = sd)
  p_high <- stats::pnorm(upper, mean = mean, sd = sd)
  if (any(p_low >= p_high)) stop("Invalid truncation bounds")
  u <- stats::runif(n, min = p_low, max = p_high)
  stats::qnorm(u, mean = mean, sd = sd)
}

# Sample task structures from observed data
sample_task_structure <- function(data_long, n_t1, n_t2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (n_t2 > n_t1) stop("n_t2 must be <= n_t1")

  df <- data_long |>
    dplyr::mutate(session = ifelse(timepoint == 0, 1, 2))

  t1_templates <- df |>
    dplyr::filter(session == 1) |>
    dplyr::distinct(subID)

  t2_templates <- df |>
    dplyr::filter(session == 2) |>
    dplyr::distinct(subID)

  if (nrow(t1_templates) == 0) stop("No timepoint 1 data found")
  if (n_t2 > 0 && nrow(t2_templates) == 0) stop("No timepoint 2 data found")

  t1_ids <- sample(t1_templates$subID, size = n_t1, replace = TRUE)
  sim_ids <- paste0("sim_", sprintf("%03d", seq_len(n_t1)))

  t1_out <- dplyr::bind_rows(lapply(seq_along(t1_ids), function(i) {
    df |>
      dplyr::filter(subID == t1_ids[i], session == 1) |>
      dplyr::mutate(subID = sim_ids[i], timepoint = 0)
  }))

  if (n_t2 == 0) return(t1_out)

  t2_ids <- sample(t2_templates$subID, size = n_t2, replace = TRUE)
  t2_sim_ids <- sample(sim_ids, size = n_t2, replace = FALSE)

  t2_out <- dplyr::bind_rows(lapply(seq_along(t2_ids), function(i) {
    df |>
      dplyr::filter(subID == t2_ids[i], session == 2) |>
      dplyr::mutate(subID = t2_sim_ids[i], timepoint = 1)
  }))

  # Ensure intervention condition is consistent across sessions
  t1_groups <- t1_out |>
    dplyr::distinct(subID, group)

  t2_out <- t2_out |>
    dplyr::select(-group) |>
    dplyr::left_join(t1_groups, by = "subID")

  dplyr::bind_rows(t1_out, t2_out)
}

# Sample participant parameters from truncated distributions
sample_param_draws <- function(ids_df,
                               mean_rew = c(0, 0),
                               mean_eff = c(0, 0),
                               mean_alpha = c(0, 0),
                               sd_rew = 1.5,
                               sd_eff = 3,
                               sd_alpha = 1.5,
                               bounds_rew = c(-1, 3),
                               bounds_eff = c(-1, 9),
                               include_alpha = FALSE,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  ids <- ids_df |>
    dplyr::distinct(subID, timepoint) |>
    dplyr::mutate(session = ifelse(timepoint == 0, 1, 2))

  ids <- ids |>
    dplyr::rowwise() |>
    dplyr::mutate(
      rewSens = rtruncnorm(1, mean = mean_rew[session], sd = sd_rew,
                           lower = bounds_rew[1], upper = bounds_rew[2]),
      effSens = rtruncnorm(1, mean = mean_eff[session], sd = sd_eff,
                           lower = bounds_eff[1], upper = bounds_eff[2]),
      alpha = if (include_alpha) stats::rnorm(1, mean = mean_alpha[session], sd = sd_alpha) else 0
    ) |>
    dplyr::ungroup() |>
    dplyr::select(subID, timepoint, rewSens, effSens, alpha)

  ids
}

# Simulate choices using provided task structure and parameters
simulate_rew_eff_from_template <- function(template_df,
                                           params_df,
                                           include_alpha = FALSE,
                                           recode_effort = FALSE,
                                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  df_sim <- template_df |>
    dplyr::mutate(session = ifelse(timepoint == 0, 1, 2))

  if (recode_effort) {
    df_sim <- df_sim |>
      dplyr::mutate(
        rew1_tmp = trialReward1,
        rew2_tmp = trialReward2,
        eff1_tmp = trialEffortPropMax1,
        eff2_tmp = trialEffortPropMax2,
        swap_effort = eff2_tmp < eff1_tmp,
        trialReward1 = dplyr::if_else(swap_effort, rew2_tmp, rew1_tmp),
        trialReward2 = dplyr::if_else(swap_effort, rew1_tmp, rew2_tmp),
        trialEffortPropMax1 = dplyr::if_else(swap_effort, eff2_tmp, eff1_tmp),
        trialEffortPropMax2 = dplyr::if_else(swap_effort, eff1_tmp, eff2_tmp)
      ) |>
      dplyr::select(-rew1_tmp, -rew2_tmp, -eff1_tmp, -eff2_tmp, -swap_effort)
  }

  df_sim <- df_sim |>
    dplyr::left_join(params_df, by = c("subID", "timepoint")) |>
    dplyr::mutate(
      v1 = rewSens * trialReward1 - effSens * trialEffortPropMax1,
      v2 = rewSens * trialReward2 - effSens * trialEffortPropMax2,
      p_choose2 = if (include_alpha) stats::plogis(alpha + v2 - v1) else stats::plogis(v2 - v1),
      choice01 = stats::rbinom(dplyr::n(), size = 1, prob = p_choose2),
      choice = dplyr::if_else(choice01 == 1, "route 2", "route 1")
    )

  df_sim
}

# Extract draws and return tidy factor data.frame

extract_draws <- function(draws, gathered, hdi = 0.95) {
  mod_type <- NULL
  if (any(grepl("delta_rew|delta_eff|delta_alpha", gathered[[".variable"]]))) {
    int_eff_entries <- list()
    if (all(c("delta_rewSens", "rewSens_int") %in% names(draws))) {
      draws <- draws |>
        dplyr::mutate(
          delta_rewSens_int = delta_rewSens + rewSens_int
        )
      int_eff_entries[[length(int_eff_entries) + 1]] <- tibble::tibble(
        .variable = rep("delta_rewSens_int", each = nrow(draws)),
        .value = draws$delta_rewSens_int
      )
    }
    if (all(c("delta_effSens", "effSens_int") %in% names(draws))) {
      draws <- draws |>
        dplyr::mutate(
          delta_effSens_int = delta_effSens + effSens_int
        )
      int_eff_entries[[length(int_eff_entries) + 1]] <- tibble::tibble(
        .variable = rep("delta_effSens_int", each = nrow(draws)),
        .value = draws$delta_effSens_int
      )
    }
    if (all(c("delta_alpha", "alpha_int") %in% names(draws))) {
      draws <- draws |>
        dplyr::mutate(
          delta_alpha_int = delta_alpha + alpha_int
        )
      int_eff_entries[[length(int_eff_entries) + 1]] <- tibble::tibble(
        .variable = rep("delta_alpha_int", each = nrow(draws)),
        .value = draws$delta_alpha_int
      )
    }
    int_eff_df <- dplyr::bind_rows(int_eff_entries)
    mod_type <- "reff"
  } else if (any(grepl("delta_int", gathered[[".variable"]]))) {
    draws <- draws |>
      dplyr::mutate(
        delta_internal_neg_int = delta_internal_neg + theta_int_internal_neg,
        delta_internal_pos_int = delta_internal_pos + theta_int_internal_pos,
        delta_global_neg_int = delta_global_neg + theta_int_global_neg,
        delta_global_pos_int = delta_global_pos + theta_int_global_pos
      )
    int_eff_df <- tibble::tibble(
      .variable = rep(
        c(
          "delta_internal_neg_int", "delta_internal_pos_int",
          "delta_global_neg_int", "delta_global_pos_int"
        ),
        each = nrow(draws)
      ),
      .value = c(
        draws$delta_internal_neg_int,
        draws$delta_internal_pos_int,
        draws$delta_global_neg_int,
        draws$delta_global_pos_int
      )
    )
    mod_type <- "cattr"
  } else {
    int_eff_df <- tibble::tibble()
    mod_type <- ifelse(
      any(grepl("rew|eff|alpha", gathered[[".variable"]])), "reff", "cattr"
    )
  }
  lvl_map <- c(
    "mu_effSens[1]" = "mean **effort sensitivity**:<br>timepoint 1",
    "mu_effSens[2]" = "mean **effort sensitivity**:<br>timepoint 2",
    "mu_rewSens[1]" = "mean **reward sensitivity**:<br>timepoint 1",
    "mu_rewSens[2]" = "mean **reward sensitivity**:<br>timepoint 2",
    "mu_alpha[1]" = "mean **acceptance bias**:<br>timepoint 1",
    "mu_alpha[2]" = "mean **acceptance bias**:<br>timepoint 2",
    "effSens_int" = "behavioural activation vs. cognitive restructuring:<br>effect on **effort sensitivity**",
    "rewSens_int" = "behavioural activation vs. cognitive restructuring:<br>effect on **reward sensitivity**",
    "alpha_int" = "behavioural activation vs. cognitive restructuring:<br>effect on **acceptance bias**",
    "effSens_ncpl" = "effect of non-completion:<br>**effort sensitivity**",
    "rewSens_ncpl" = "effect of non-completion:<br>**reward sensitivity**",
    "alpha_ncpl" = "effect of non-completion:<br>**acceptance bias**",
    "delta_effSens" = "change in **effort sensitivity**", #:<br>restructuring",
    "delta_effSens_int" = "change in **effort sensitivity**", #:<br>goal-setting",
    "delta_rewSens" = "change in **reward sensitivity**", #:<br>restructuring",
    "delta_rewSens_int" = "change in **reward sensitivity**", #:<br>goal-setting",
    "delta_alpha" = "change in **acceptance bias**",
    "delta_alpha_int" = "change in **acceptance bias**",
    "phq_rewSens_t1_group_cr" = "association between **reward sensitivity**<br>and PHQ-9 score, timepoint 1",
    "phq_rewSens_t1_group_ba" = "association between **reward sensitivity**<br>and PHQ-9 score, timepoint 1",
    "phq_rewSens_t1_diff" =
      "behavioural activation vs. cognitive restructuring: association between<br>**reward sensitivity** & PHQ-9 score, timepoint 1", #nolint
    "phq_rewSens_t2_group_cr" = "association between **reward sensitivity**<br>and PHQ-9 score, timepoint 2",
    "phq_rewSens_t2_group_ba" = "association between **reward sensitivity**<br>and PHQ-9 score, timepoint 2",
    "phq_rewSens_t2_diff" =
      "behavioural activation vs. cognitive restructuring: association between<br>**reward sensitivity** & PHQ-9 score, timepoint 2", #nolint
    "phq_rewSens_t2_change_cr" = "pre-post change in association between<br>**reward sensitivity** & PHQ-9 score",
    "phq_rewSens_t2_change_ba" = "pre-post change in association between<br>**reward sensitivity** & PHQ-9 score",
    "phq_rewSens_t2_change_diff" =
      "behavioural activation vs. cognitive restructuring:<br>pre-post change in association between<br>**reward sensitivity** & PHQ-9 score", #nolint
    "phq_effSens_t1_group_cr" = "association between **effort sensitivity**<br>and PHQ-9 score, timepoint 1",
    "phq_effSens_t1_group_ba" = "association between **effort sensitivity**<br>and PHQ-9 score, timepoint 1",
    "phq_effSens_t1_diff" =
      "behavioural activation vs. cognitive restructuring: association between<br>**effort sensitivity** & PHQ-9 score, timepoint 1", #nolint
    "phq_effSens_t2_group_cr" = "association between **effort sensitivity**<br>and PHQ-9 score, timepoint 2",
    "phq_effSens_t2_group_ba" = "association between **effort sensitivity**<br>and PHQ-9 score, timepoint 2",
    "phq_effSens_t2_diff" =
      "behavioural activation vs. cognitive restructuring: association between<br>**effort sensitivity** & PHQ-9 score, timepoint 2", #nolint
    "phq_effSens_t2_change_cr" = "pre-post change in association between<br>**effort sensitivity** & PHQ-9 score",
    "phq_effSens_t2_change_ba" = "pre-post change in association between<br>**effort sensitivity** & PHQ-9 score",
    "phq_effSens_t2_change_diff" =
      "behavioural activation vs. cognitive restructuring:<br>pre-post change in association between<br>**effort sensitivity** & PHQ-9 score", #nolint
    "phq_alpha_t1_group_cr" = "association between **acceptance bias**<br>and PHQ-9 score, timepoint 1",
    "phq_alpha_t1_group_ba" = "association between **acceptance bias**<br>and PHQ-9 score, timepoint 1",
    "phq_alpha_t1_diff" =
      "behavioural activation vs. cognitive restructuring: association between<br>**acceptance bias** & PHQ-9 score, timepoint 1", #nolint
    "phq_alpha_t2_group_cr" = "association between **acceptance bias**<br>and PHQ-9 score, timepoint 2",
    "phq_alpha_t2_group_ba" = "association between **acceptance bias**<br>and PHQ-9 score, timepoint 2",
    "phq_alpha_t2_diff" =
      "behavioural activation vs. cognitive restructuring: association between<br>**acceptance bias** & PHQ-9 score, timepoint 2", #nolint
    "phq_alpha_t2_change_cr" = "pre-post change in association between<br>**acceptance bias** & PHQ-9 score",
    "phq_alpha_t2_change_ba" = "pre-post change in association between<br>**acceptance bias** & PHQ-9 score",
    "phq_alpha_t2_change_diff" =
      "behavioural activation vs. cognitive restructuring:<br>pre-post change in association between<br>**acceptance bias** & PHQ-9 score", #nolint
    "mu_internal_theta_neg[1]" = "mean \u03B8 **internal-negative**:<br>timepoint 1",
    "mu_internal_theta_neg[2]" = "mean \u03B8 **internal-negative**:<br>timepoint 2",
    "mu_internal_theta_pos[1]" = "mean \u03B8 **internal-positive**:<br>timepoint 1",
    "mu_internal_theta_pos[2]" = "mean \u03B8 **internal-positive**:<br>timepoint 2",
    "mu_global_theta_neg[1]" = "mean \u03B8 **global-negative**:<br>timepoint 1",
    "mu_global_theta_neg[2]" = "mean \u03B8 **global-negative**:<br>timepoint 2",
    "mu_global_theta_pos[1]" = "mean \u03B8 **global-positive**:<br>timepoint 1",
    "mu_global_theta_pos[2]" = "mean \u03B8 **global-positive**:<br>timepoint 2",
    "theta_int_internal_neg" = "cognitive restructuring vs. behavioural activation:<br> \u03B8 **internal-negative**",
    "theta_int_internal_pos" = "cognitive restructuring vs. behavioural activation:<br> \u03B8 **internal-positive**",
    "theta_int_global_neg" = "cognitive restructuring vs. behavioural activation:<br> \u03B8 **global-negative**",
    "theta_int_global_pos" = "cognitive restructuring vs. behavioural activation:<br> \u03B8 **global-positive**",
    "theta_ncompl_internal_neg" = "effect of non-completion:<br> \u03B8 **internal-negative**",
    "theta_ncompl_internal_pos" = "effect of non-completion:<br> \u03B8 **internal-positive**",
    "theta_ncompl_global_neg" = "effect of non-completion:<br> \u03B8 **global-negative**",
    "theta_ncompl_global_pos" = "effect of non-completion:<br> \u03B8 **global-positive**",
    "delta_internal_neg" = "change in \u03B8 **internal-negative**", #:<br>goal-setting",
    "delta_internal_neg_int" = "change in \u03B8 **internal-negative**", #:<br>restructuring",
    "delta_internal_pos" = "change in \u03B8 **internal-positive**", #:<br>goal-setting",
    "delta_internal_pos_int" = "change in \u03B8 **internal-positive**", #:<br>restructuring",
    "delta_global_neg" = "change in \u03B8 **global-negative**", #:<br>goal-setting",
    "delta_global_neg_int" = "change in \u03B8 **global-negative**", #:<br>restructuring",
    "delta_global_pos" = "change in \u03B8 **global-positive**", #:<br>goal-setting",
    "delta_global_pos_int" = "change in \u03B8 **global-positive**", #:<br>restructuring"
    "phq_int_neg_t1_group_cr" = "association between \u03B8 **internal-negative**<br>and PHQ-9 score, timepoint 1",
    "phq_int_neg_t1_group_ba" = "association between \u03B8 **internal-negative**<br>and PHQ-9 score, timepoint 1",
    "phq_int_neg_t1_diff" =
      "cognitive restructuring vs. behavioural activation: association between<br>\u03B8 **internal-negative** & PHQ-9 score, timepoint 1", #nolint
    "phq_int_neg_t2_group_cr" = "association between \u03B8 **internal-negative**<br>and PHQ-9 score, timepoint 2",
    "phq_int_neg_t2_group_ba" = "association between \u03B8 **internal-negative**<br>and PHQ-9 score, timepoint 2",
    "phq_int_neg_t2_diff" =
      "cognitive restructuring vs. behavioural activation: association between<br>\u03B8 **internal-negative** & PHQ-9 score, timepoint 2", #nolint
    "phq_int_neg_t2_change_cr" = "pre-post change in association between<br>\u03B8 **internal-negative** & PHQ-9 score",
    "phq_int_neg_t2_change_ba" = "pre-post change in association between<br>\u03B8 **internal-negative** & PHQ-9 score",
    "phq_int_neg_t2_change_diff" =
      "cognitive restructuring vs. behavioural activation:<br>pre-post change in association between<br>\u03B8 **internal-negative** & PHQ-9 score", #nolint
    "phq_int_pos_t1_group_cr" = "association between \u03B8 **internal-positive**<br>and PHQ-9 score, timepoint 1",
    "phq_int_pos_t1_group_ba" = "association between \u03B8 **internal-positive**<br>and PHQ-9 score, timepoint 1",
    "phq_int_pos_t1_diff" =
      "cognitive restructuring vs. behavioural activation: association between<br>\u03B8 **internal-positive** & PHQ-9 score, timepoint 1", #nolint
    "phq_int_pos_t2_group_cr" = "association between \u03B8 **internal-positive**<br>and PHQ-9 score, timepoint 2",
    "phq_int_pos_t2_group_ba" = "association between \u03B8 **internal-positive**<br>and PHQ-9 score, timepoint 2",
    "phq_int_pos_t2_diff" =
      "cognitive restructuring vs. behavioural activation: association between<br>\u03B8 **internal-positive** & PHQ-9 score, timepoint 2", #nolint
    "phq_int_pos_t2_change_cr" = "pre-post change in association between<br>\u03B8 **internal-positive** & PHQ-9 score",
    "phq_int_pos_t2_change_ba" = "pre-post change in association between<br>\u03B8 **internal-positive** & PHQ-9 score",
    "phq_int_pos_t2_change_diff" =
      "cognitive restructuring vs. behavioural activation:<br>pre-post change in association between<br>\u03B8 **internal-positive** & PHQ-9 score", #nolint
    "phq_glob_neg_t1_group_cr" = "association between \u03B8 **global-negative**<br>and PHQ-9 score, timepoint 1",
    "phq_glob_neg_t1_group_ba" = "association between \u03B8 **global-negative**<br>and PHQ-9 score, timepoint 1",
    "phq_glob_neg_t1_diff" =
      "cognitive restructuring vs. behavioural activation: association between<br>\u03B8 **global-negative** & PHQ-9 score, timepoint 1", #nolint
    "phq_glob_neg_t2_group_cr" = "association between \u03B8 **global-negative**<br>and PHQ-9 score, timepoint 2",
    "phq_glob_neg_t2_group_ba" = "association between \u03B8 **global-negative**<br>and PHQ-9 score, timepoint 2",
    "phq_glob_neg_t2_diff" =
      "cognitive restructuring vs. behavioural activation: association between<br>\u03B8 **global-negative** & PHQ-9 score, timepoint 2", #nolint
    "phq_glob_neg_t2_change_cr" = "pre-post change in association between<br>\u03B8 **global-negative** & PHQ-9 score",
    "phq_glob_neg_t2_change_ba" = "pre-post change in association between<br>\u03B8 **global-negative** & PHQ-9 score",
    "phq_glob_neg_t2_change_diff" =
      "cognitive restructuring vs. behavioural activation:<br>pre-post change in association between<br>\u03B8 **global-negative** & PHQ-9 score", #nolint
    "phq_glob_pos_t1_group_cr" = "association between \u03B8 **global-positive**<br>and PHQ-9 score, timepoint 1",
    "phq_glob_pos_t1_group_ba" = "association between \u03B8 **global-positive**<br>and PHQ-9 score, timepoint 1",
    "phq_glob_pos_t1_diff" =
      "cognitive restructuring vs. behavioural activation: association between<br>\u03B8 **global-positive** & PHQ-9 score, timepoint 1", #nolint
    "phq_glob_pos_t2_group_cr" = "association between \u03B8 **global-positive**<br>and PHQ-9 score, timepoint 2",
    "phq_glob_pos_t2_group_ba" = "association between \u03B8 **global-positive**<br>and PHQ-9 score, timepoint 2",
    "phq_glob_pos_t2_diff" =
      "cognitive restructuring vs. behavioural activation: association between<br>\u03B8 **global-positive** & PHQ-9 score, timepoint 2", #nolint
    "phq_glob_pos_t2_change_cr" = "pre-post change in association between<br>\u03B8 **global-positive** & PHQ-9 score",
    "phq_glob_pos_t2_change_ba" = "pre-post change in association between<br>\u03B8 **global-positive** & PHQ-9 score",
    "phq_glob_pos_t2_change_diff" =
      "cognitive restructuring vs. behavioural activation:<br>pre-post change in association between<br>\u03B8 **global-positive** & PHQ-9 score" #nolint
  )
  df_ret <- gathered |>
    dplyr::bind_rows(int_eff_df) |>
    dplyr::mutate(
      var_type = ifelse(
        grepl("^delta.*_int$", .variable) | grepl("^theta_int_", .variable),
        "intervention effect", "group mean"
      )
    )

  if (mod_type == "reff") {
    df_ret <- df_ret |> dplyr::mutate(
      eff_rew = dplyr::case_when(
        grepl("eff", .variable) ~ "effort",
        grepl("rew", .variable) ~ "reward",
        grepl("alpha", .variable) ~ "acceptance bias",
        TRUE ~ "other"
      ),
      var_group = factor(
        dplyr::case_when(
          grepl("_int$|_ba$", .variable) ~ "behavioural activation",
          grepl("^delta|_cr$", .variable) ~ "cognitive restructuring",
          grepl("ncpl|ncompl", .variable) ~ "non-completion",
          grepl("diff$", .variable) ~ "group difference",
          TRUE ~ "group mean"
        ),
        levels = c("behavioural activation", "cognitive restructuring", "non-completion", "group mean")
      ),
      prim_int = ifelse(var_group == "behavioural activation", TRUE, FALSE)
    )
  } else if (mod_type == "cattr") {
    df_ret <- df_ret |> dplyr::mutate(
      pos_neg = ifelse(grepl("pos", .variable), "positive", "negative"),
      int_glob = ifelse(grepl("internal", .variable), "internal", "global"),
      var_group = factor(
        dplyr::case_when(
          grepl("theta_int_|_int$|_cr$", .variable) ~ "cognitive restructuring",
          grepl("^delta|_ba$", .variable) ~ "behavioural activation",
          grepl("ncpl|ncompl", .variable) ~ "non-completion",
          grepl("diff$", .variable) ~ "group difference",
          TRUE ~ "group mean"
        ),
        levels = c("cognitive restructuring", "behavioural activation", "non-completion", "group mean")
      ),
      prim_int = ifelse(var_group == "cognitive restructuring", TRUE, FALSE),
    )
  }
  # Subset the lvl_map based on the column names in "gathered"
  lvls <- intersect(names(lvl_map), colnames(draws))

  # Use the subsetted lvl_map to match levels to labels
  df_ret$variable <- factor(
    df_ret[[".variable"]],
    levels = lvls,
    labels = lvl_map[lvls]
  )

  hdi_ret <- df_ret |>
    dplyr::group_by(variable, var_group) |>
    dplyr::mutate(
      mean = mean(.value),
      hdi_low = signif(bayestestR::hdi(.value, ci = 0.95)[[2]], digits = 3),
      hdi_high = signif(bayestestR::hdi(.value, ci = 0.95)[[3]], digits = 3),
      pd = signif(bayestestR::p_direction(.value, method = "direct")[[2]], digits = 3)
    ) |>
    dplyr::distinct(variable, mean, hdi_low, hdi_high, pd, var_group)

  list("df" = df_ret, "hdi" = hdi_ret)
}

get_empirical_means_from_fit <- function(fit, stan_ls, include_alpha = FALSE) {
  vars <- c("rewSens", "effSens")
  if (include_alpha) vars <- c(vars, "alpha")

  draws <- fit$draws(format = "df", variables = vars)
  pars <- extract_indiv_pars(draws, stan_ls)

  means <- pars |>
    dplyr::filter(timepoint %in% c("t1", "t2"), variable %in% vars) |>
    dplyr::group_by(variable, timepoint) |>
    dplyr::summarise(mean_val = mean(mean), .groups = "drop")

  get_mean <- function(var, tpt) {
    means |>
      dplyr::filter(variable == var, timepoint == tpt) |>
      dplyr::pull(mean_val) |>
      dplyr::first()
  }

  list(
    rew = c(get_mean("rewSens", "t1"), get_mean("rewSens", "t2")),
    eff = c(get_mean("effSens", "t1"), get_mean("effSens", "t2")),
    alpha = if (include_alpha) c(get_mean("alpha", "t1"), get_mean("alpha", "t2")) else c(0, 0)
  )
}

# Decompose the pre->post change in model-predicted acceptance into additive
# (Shapley) contributions from the acceptance-bias, reward-sensitivity and
# effort-sensitivity changes -- all on the single interpretable scale of
# P(accept the higher-effort option). Also returns natural-scale group-level
# deltas, which fixes the latent-probit (delta_rewSens/effSens) vs natural-scale
# (delta_alpha) mismatch in the Stan generated quantities.
#
# Works for the accept-bias model (alpha present) and the standard model (no
# alpha -> its contribution is identically 0). Run it on BOTH fits and compare
# the `reward` contribution: the amount reward loses when the intercept is added
# is the part of the "reward-sensitivity" change that is really acceptance bias.
#
# `stan_ls` MUST be the data list the fit was estimated on -- it supplies the
# trial design used as the fixed reference set (so the change reflects
# parameters, not trial composition). Predictions are for the "typical"
# participant (individual offsets = 0), i.e. the group-level effect. The two arms
# are CR (baseline: mu[2] only) and BA (mu[2] + intervention effect).
reff_delta_decomp <- function(fit, stan_ls, hdi = 0.95) {
  # Stan's Phi_approx, so natural-scale values match the fitted model
  phi_approx <- function(x) stats::plogis(0.07056 * x^3 + 1.5976 * x)
  nat_rew <- function(x) -1 + phi_approx(x) * 4
  nat_eff <- function(x) -1 + phi_approx(x) * 10

  avail <- fit$metadata()$stan_variables
  if (!all(c("mu_rewSens", "mu_effSens") %in% avail)) {
    stop("reff_delta_decomp(): fit is missing mu_rewSens / mu_effSens.")
  }
  has_alpha <- "mu_alpha" %in% avail
  has_int   <- all(c("rewSens_int", "effSens_int") %in% avail)

  # ---- reference trial design: pooled (dRew, dEff) over all valid trials ------
  dRew <- numeric(0)
  dEff <- numeric(0)
  for (p in seq_len(stan_ls$nPpts)) {
    for (t in 1:2) {
      nt <- stan_ls$nT_ppts[p, t]
      if (nt > 0) {
        dRew <- c(dRew, stan_ls$rew2[p, t, 1:nt] - stan_ls$rew1[p, t, 1:nt])
        dEff <- c(dEff, stan_ls$eff2[p, t, 1:nt] - stan_ls$eff1[p, t, 1:nt])
      }
    }
  }
  if (length(dRew) == 0) stop("reff_delta_decomp(): no valid trials found in stan_ls.")

  # Canonicalise every reference trial to the higher-effort framing, so
  # "acceptance" = choosing the higher-effort option for BOTH models. rewSens /
  # effSens are coding-invariant (relabelling the two options and flipping the
  # choice leaves the likelihood unchanged), so this is valid whether or not
  # `stan_ls` was built with recode_effort = TRUE -- and it is REQUIRED for a
  # non-recoded stan_ls, where option 2 is an arbitrary route: there dRew/dEff
  # are symmetric about 0, mean P(choose option 2) is ~0.5 at both timepoints,
  # and every contribution collapses to ~0. After flipping, dEff >= 0 and dRew is
  # the reward premium of the higher-effort option, so the sensitivities act on
  # the acceptance rate. (Effort ties carry no higher-effort option and are dropped.)
  keep0 <- dEff != 0
  dRew <- dRew[keep0]
  dEff <- dEff[keep0]
  flip <- dEff < 0
  dRew[flip] <- -dRew[flip]
  dEff[flip] <- -dEff[flip]
  if (length(dRew) == 0) stop("reff_delta_decomp(): all reference trials were effort ties.")

  # ---- posterior draws of the group-level parameters -------------------------
  vars <- intersect(
    c("mu_rewSens", "mu_effSens", "mu_alpha", "rewSens_int", "effSens_int", "alpha_int"),
    avail
  )
  dr <- fit$draws(variables = vars, format = "df")
  nd <- nrow(dr)
  z0 <- rep(0, nd)

  # t1 (baseline) natural-scale parameters -- shared by both arms
  r1 <- nat_rew(dr[["mu_rewSens[1]"]])
  e1 <- nat_eff(dr[["mu_effSens[1]"]])
  a1 <- if (has_alpha) dr[["mu_alpha[1]"]] else z0

  # mean predicted acceptance over the reference trials for one parameter config
  # (eta = alpha + rewSens*dRew - effSens*dEff, matching the Stan likelihood)
  vconf <- function(a, r, e) {
    eta <- matrix(a, nrow = nd, ncol = length(dRew)) + outer(r, dRew) - outer(e, dEff)
    rowMeans(stats::plogis(eta))
  }

  decomp_arm <- function(arm) {
    is_ba <- identical(arm, "BA")
    r2 <- nat_rew(dr[["mu_rewSens[2]"]] + (if (is_ba && has_int) dr[["rewSens_int"]] else z0))
    e2 <- nat_eff(dr[["mu_effSens[2]"]] + (if (is_ba && has_int) dr[["effSens_int"]] else z0))
    a2 <- if (has_alpha) {
      dr[["mu_alpha[2]"]] + (if (is_ba && "alpha_int" %in% avail) dr[["alpha_int"]] else z0)
    } else {
      z0
    }

    # Predicted acceptance with any subset of the three parameters moved to their
    # t2 value and the rest held at t1. Each requested combination is computed
    # once and cached, because the Shapley formulas below reuse them.
    acc_cache <- new.env(parent = emptyenv())
    p_accept <- function(alpha_t2, reward_t2, effort_t2) {
      key <- paste0(alpha_t2 + 0L, reward_t2 + 0L, effort_t2 + 0L)
      if (!exists(key, envir = acc_cache, inherits = FALSE)) {
        assign(key, vconf(
          if (alpha_t2) a2 else a1,
          if (reward_t2) r2 else r1,
          if (effort_t2) e2 else e1
        ), envir = acc_cache)
      }
      get(key, envir = acc_cache, inherits = FALSE)
    }

    # A parameter's Shapley value is the average of its marginal effect on
    # acceptance over every order the three parameters could be switched from t1
    # to t2. With three parameters the ordering weights are 1/3 when it switches
    # first or last and 1/6 in each of the two middle positions (they sum to 1).
    # Each line contrasts "this parameter at t2" with "this parameter at t1",
    # holding the already-switched parameters fixed.
    shap_alpha <-
      1 / 3 * (p_accept(TRUE, FALSE, FALSE) - p_accept(FALSE, FALSE, FALSE)) + # switched first
      1 / 6 * (p_accept(TRUE, TRUE, FALSE) - p_accept(FALSE, TRUE, FALSE)) + # after reward
      1 / 6 * (p_accept(TRUE, FALSE, TRUE) - p_accept(FALSE, FALSE, TRUE)) + # after effort
      1 / 3 * (p_accept(TRUE, TRUE, TRUE) - p_accept(FALSE, TRUE, TRUE)) # switched last
    shap_reward <-
      1 / 3 * (p_accept(FALSE, TRUE, FALSE) - p_accept(FALSE, FALSE, FALSE)) + # switched first
      1 / 6 * (p_accept(TRUE, TRUE, FALSE) - p_accept(TRUE, FALSE, FALSE)) + # after acceptance bias
      1 / 6 * (p_accept(FALSE, TRUE, TRUE) - p_accept(FALSE, FALSE, TRUE)) + # after effort
      1 / 3 * (p_accept(TRUE, TRUE, TRUE) - p_accept(TRUE, FALSE, TRUE)) # switched last
    shap_effort <-
      1 / 3 * (p_accept(FALSE, FALSE, TRUE) - p_accept(FALSE, FALSE, FALSE)) + # switched first
      1 / 6 * (p_accept(TRUE, FALSE, TRUE) - p_accept(TRUE, FALSE, FALSE)) + # after acceptance bias
      1 / 6 * (p_accept(FALSE, TRUE, TRUE) - p_accept(FALSE, TRUE, FALSE)) + # after reward
      1 / 3 * (p_accept(TRUE, TRUE, TRUE) - p_accept(TRUE, TRUE, FALSE)) # switched last

    total <- p_accept(TRUE, TRUE, TRUE) - p_accept(FALSE, FALSE, FALSE)

    # Shapley contributions are exactly additive; guard against a coding slip
    if (max(abs(shap_alpha + shap_reward + shap_effort - total)) > 1e-8) {
      stop("reff_delta_decomp(): Shapley contributions do not sum to the total change.")
    }

    tibble::tibble(
      arm = arm,
      alpha = shap_alpha, reward = shap_reward, effort = shap_effort, total = total,
      # natural-scale group deltas -- note these keep parameter-specific units
      # (per-coin, per-effort, logit), unlike the P(accept) contributions above
      nat_alpha = if (has_alpha) a2 - a1 else z0,
      nat_reward = r2 - r1, nat_effort = e2 - e1
    )
  }

  arms <- if (has_int) c("CR", "BA") else "CR"
  raw <- dplyr::bind_rows(lapply(arms, decomp_arm))

  summ <- function(df, cols) {
    df |>
      tidyr::pivot_longer(dplyr::all_of(cols), names_to = "component", values_to = "val") |>
      dplyr::group_by(arm, component) |>
      dplyr::summarise(
        mean = mean(val),
        hdi_low = bayestestR::hdi(val, ci = hdi)[[2]],
        hdi_high = bayestestR::hdi(val, ci = hdi)[[3]],
        pd = bayestestR::p_direction(val)[[2]],
        .groups = "drop"
      )
  }

  list(
    decomp = summ(raw, c("alpha", "reward", "effort", "total")),        # scale: ΔP(accept)
    nat_delta = summ(raw, c("nat_alpha", "nat_reward", "nat_effort")),  # natural param units
    draws = raw,
    has_alpha = has_alpha
  )
}

# Overlay two reff_delta_decomp() results (e.g. the main additive model and the
# extended acceptance-bias model) on ONE panel, so the reallocation of the
# reward-sensitivity contribution into acceptance bias is directly visible as a
# shift between adjacent marks, rather than requiring readers to compare two
# separate x-axes side by side (as plotting each decomposition separately would).
# Colour still encodes arm (BA/CR), matching every other panel in this figure;
# model is encoded by shape + errorbar linetype instead (merged into one legend,
# since ggplot combines guides when two aesthetics map the same variable with
# identical labels). Acceptance bias is only drawn for whichever model(s)
# actually have it (`has_alpha`); a model without it contributes no row rather
# than a spurious point at zero.
plot_delta_decomp_compare <- function(dec_main, dec_ext,
                                      main_label = "main model",
                                      ext_label  = "extended model",
                                      base_size  = 14,
                                      arm_cols   = NULL,
                                      model_shapes = c(7, 10),
                                      model_linetypes = c("solid", "dashed"),
                                      include_total = TRUE) {
  full_levels <- c("alpha", "reward", "effort", "total")
  full_labels <- c("acceptance bias", "reward sensitivity", "effort sensitivity", "total change")

  prep_one <- function(decomp_obj, model_lab) {
    keep <- c(if (isTRUE(decomp_obj$has_alpha)) "alpha", "reward", "effort", if (include_total) "total")
    decomp_obj$decomp |>
      dplyr::filter(component %in% keep) |>
      dplyr::mutate(model = model_lab)
  }
  d <- dplyr::bind_rows(prep_one(dec_main, main_label), prep_one(dec_ext, ext_label))

  comp_levels <- full_levels[full_levels %in% unique(d$component)]
  comp_labels <- full_labels[full_levels %in% unique(d$component)]
  arm_lab <- c(CR = "cognitive restructuring", BA = "behavioural activation")

  d <- d |>
    dplyr::mutate(
      component = factor(component, levels = rev(comp_levels), labels = rev(comp_labels)),
      arm = factor(unname(arm_lab[arm]), levels = unname(arm_lab)),
      model = factor(model, levels = c(main_label, ext_label)),
      dplyr::across(c(mean, hdi_low, hdi_high), ~ .x * 100)
    ) |>
    droplevels()

  # explicit dodge order (arm-major, model-minor) so each arm's main/extended
  # pair sits adjacent -- the direct comparison the user is after -- rather than
  # whatever order ggplot would infer from combining colour/shape/linetype
  grp_levels <- as.vector(t(outer(levels(d$arm), levels(d$model), paste, sep = " | ")))
  d$dodge_grp <- factor(paste(d$arm, d$model, sep = " | "), levels = grp_levels)

  if (is.null(arm_cols)) {
    arm_cols <- setNames(MetBrewer::met.brewer("Archambault")[c(4, 1)], unname(arm_lab))
  }
  shape_vals <- setNames(model_shapes, levels(d$model))
  lty_vals   <- setNames(model_linetypes, levels(d$model))
  pd <- ggplot2::position_dodge(width = 0.75)

  ggplot2::ggplot(d, ggplot2::aes(x = mean, y = component, group = dodge_grp, colour = arm)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    # separate the total (always the bottom row) from its components
    (if (include_total) {
      ggplot2::geom_hline(yintercept = 1.5, linetype = "dotted", colour = "grey75", linewidth = 0.5)
    }) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = hdi_low, xmax = hdi_high, linetype = model),
      orientation = "y", width = 0.15, position = pd, linewidth = 0.9
    ) +
    ggplot2::geom_point(ggplot2::aes(shape = model), size = 3, position = pd) +
    ggplot2::scale_colour_manual(values = arm_cols, name = NULL) +
    ggplot2::scale_shape_manual(values = shape_vals, name = NULL) +
    ggplot2::scale_linetype_manual(values = lty_vals, name = NULL) +
    ggplot2::labs(
      x = "contribution to pre → post change in acceptance (%)",
      y = NULL
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(reverse = TRUE, order = 1, position = "bottom"),
      shape = ggplot2::guide_legend(position = "inside"),
      linetype = ggplot2::guide_legend(position = "inside")
    ) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = base_size) +
    cowplot::background_grid(major = "x", minor = "none") +
    ggplot2::theme(
      legend.position.inside = c(0.025, 0.925),
      legend.justification.bottom = "centre",
      legend.background = ggplot2::element_rect(
        fill = "white", colour = "#cacaca", linewidth = 0.5
      ),
      legend.margin = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm")
    )
}

# simple plot of acceptance % vs. difference in reward and effort between routes, by group and timepoint
reff_plot <- function(
    data,
    plot_type = c("reward", "effort", "accept_reward", "accept_effort"),
    mode = c("both", "randomised", "completer_status"),
    error_type = c("sd", "se"),
    label_colour = NULL,
    ln_width = 0.8,
    err_width = 0.08,
    dodge_width = 0.28,
    point_size = 3.3,
    shape_types = c(15, 17),
    n_bins = 5,
    legend_pos = "inside",
    legend_pos_in = c(0.82, 0.12),
    legend_fnt_sc = 9,
    fnt_sz = 1,
    brew_col = "Archambault",
    col_nums = c(1, 4, 3),
    custom_pal = NULL) {

  plot_type  <- match.arg(plot_type)
  mode       <- match.arg(mode)
  error_type <- match.arg(error_type)

  to_plt <- data |>
    dplyr::mutate(
      group_plot = dplyr::case_when(
        as.logical(dplyr::coalesce(.data[["non_completed"]], FALSE)) ~ "non-completers",
        .data[["group"]] == "BA" ~ "behavioural activation",
        .data[["group"]] == "CR" ~ "cognitive restructuring",
        TRUE ~ NA_character_
      ),
      group_plot = factor(
        group_plot, levels = c("behavioural activation", "cognitive restructuring", "non-completers")
      ),
      timepoint_plot = factor(
        as.character(.data[["timepoint"]]), levels = c("0", "1"), labels = c("pre-intervention", "post-intervention")
      )
    )

  if (mode == "randomised") {
    # both timepoints; non-completers merged into their randomised group
    # (at post they simply have no data, so the mean reflects completers only)
    to_plt <- to_plt |>
      dplyr::mutate(
        group_plot = dplyr::case_when(
          .data[["group"]] == "BA" ~ "behavioural activation",
          .data[["group"]] == "CR" ~ "cognitive restructuring",
          TRUE ~ NA_character_
        ),
        group_plot = factor(group_plot, levels = c("behavioural activation", "cognitive restructuring"))
      )
  } else if (mode == "completer_status") {
    # pre-intervention only; completers vs. non-completers regardless of group
    to_plt <- to_plt |>
      dplyr::filter(.data[["timepoint"]] == 0) |>
      dplyr::mutate(
        group_plot = dplyr::if_else(
          group_plot == "non-completers", "non-completers", "completers"
        ),
        group_plot = factor(group_plot, levels = c("completers", "non-completers"))
      )
  }

  full_pal <- MetBrewer::met.brewer(brew_col)
  if (!is.null(custom_pal)) {
    pal <- custom_pal
  } else if (is.null(col_nums)) {
    pal <- sample(full_pal, length(unique(to_plt$group_plot)))
  } else {
    pal <- full_pal[col_nums]
  }

  grp_cols <- if (mode == "both") {
    setNames(pal[1:3], c("behavioural activation", "cognitive restructuring", "non-completers"))
  } else if (mode == "randomised") {
    setNames(pal[1:2], c("behavioural activation", "cognitive restructuring"))
  } else {
    setNames(pal[1:2], c("completers", "non-completers"))
  }

  if (is.null(label_colour)) {
    label_colour <- if (plot_type == "reward") full_pal[2] else full_pal[3]
  }

  if (plot_type == "reward") {
    plot_df <- to_plt |>
      dplyr::transmute(
        id_plot = as.character(.data[["subID"]]),
        group_plot,
        timepoint_plot,
        predictor_value = as.numeric(.data[["highRewardOption"]]),
        outcome_val = as.numeric(.data[["higherRewChosen"]])
      )
  } else {
    plot_df <- to_plt |>
      dplyr::transmute(
        id_plot = as.character(.data[["subID"]]),
        group_plot,
        timepoint_plot,
        predictor_value = 100 * as.numeric(.data[["highEffortOption"]]),
        outcome_val = as.numeric(.data[["higherEffChosen"]])
      )
  }

  plot_df <- plot_df |>
    dplyr::filter(!is.na(id_plot), !is.na(outcome_val), !is.na(predictor_value), !is.na(group_plot))

  subj_df <- plot_df |>
    dplyr::group_by(id_plot, group_plot, timepoint_plot, predictor_value) |>
    dplyr::summarise(outcome_val = mean(outcome_val), .groups = "drop")

  sum_df <- subj_df |>
    dplyr::group_by(group_plot, timepoint_plot, predictor_value) |>
    dplyr::summarise(
      outcome_mean = mean(outcome_val),
      outcome_sd = stats::sd(outcome_val),
      outcome_se = outcome_sd / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    dplyr::mutate(err = if (error_type == "sd") outcome_sd else outcome_se)

  pd <- ggplot2::position_dodge(width = dodge_width)
  legend_key_h <- grid::unit(max(6, legend_fnt_sc * fnt_sz), "pt")
  legend_spacing_y <- grid::unit(max(0, 2 * fnt_sz), "pt")
  x_lab <- if (plot_type == "reward") {
    paste0("higher <span style='color:", label_colour, ";'>reward</span> route (number of coins)")
  } else {
    paste0("higher <span style='color:", label_colour, ";'>effort</span> route (% of max effort)")
  }
  title_lab <- paste0("<span style='color:", label_colour, ";'>", plot_type, "</span>")

  p <- ggplot2::ggplot()

  p <- p +
    ggplot2::geom_errorbar(
      data = sum_df,
      ggplot2::aes(
        x = predictor_value,
        y = outcome_mean,
        ymin = pmax(outcome_mean - err, 0),
        ymax = pmin(outcome_mean + err, 1),
        colour = group_plot,
        group = group_plot
      ),
      width = err_width,
      linewidth = ln_width * 0.6,
      position = pd,
      alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = sum_df,
      ggplot2::aes(
        x = predictor_value,
        y = outcome_mean,
        colour = group_plot,
        fill = group_plot,
        shape = group_plot,
        group = group_plot
      ),
      size = point_size,
      stroke = 0.4,
      position = pd
    ) + # add line through points
    ggplot2::geom_line(
      data = sum_df,
      ggplot2::aes(
        x = predictor_value,
        y = outcome_mean,
        colour = group_plot,
        group = group_plot
      ),
      linewidth = ln_width,
      position = pd
    )

  p +
    ggplot2::geom_hline(
      yintercept = 0.5, linetype = "dashed", colour = "#cacaca", alpha = 0.8
    ) +
    (if (mode != "completer_status") ggplot2::facet_wrap(~timepoint_plot, nrow = 1)) +
    ggplot2::scale_colour_manual(name = "", values = grp_cols) +
    ggplot2::scale_fill_manual(name = "", values = grp_cols) +
    ggplot2::scale_shape_manual(name = "", values = shape_types) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(),
      fill = ggplot2::guide_legend(),
      shape = ggplot2::guide_legend()
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = title_lab,
      x = x_lab,
      y = paste0("prop. accept (± ", ifelse(error_type == "sd", "s.d.", "s.e."), ")")
    ) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = 14 * fnt_sz) +
    ggplot2::theme(
      legend.key.height = legend_key_h,
      legend.spacing.y = legend_spacing_y,
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "italic", size = 11 * fnt_sz),
      plot.title = ggtext::element_markdown(
        family = "Open Sans SemiBold",
        size = 14 * fnt_sz,
        margin = ggplot2::margin(t = 0, r = 0, b = 0.01, l = 0, unit = "npc")
      ),
      axis.title.y = ggplot2::element_text(size = 12 * fnt_sz),
      axis.title.x = ggtext::element_markdown(size = 14 * fnt_sz),
      axis.text = ggplot2::element_text(size = 11 * fnt_sz),
      legend.text = ggplot2::element_text(size = legend_fnt_sc * fnt_sz),
      legend.key.size = ggplot2::unit(2 * legend_fnt_sc, "pt")
    ) +
    (
      if (legend_pos == "inside") {
        ggplot2::theme(legend.position = "inside", legend.position.inside = legend_pos_in)
      } else if (legend_pos == "bottom") {
        ggplot2::theme(
          legend.position = legend_pos,
          legend.justification.bottom = "center",
          legend.background = ggplot2::element_rect(
            fill = "white", colour = "#cacaca", linewidth = 0.5
          ),
          legend.margin = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm")
        )
      } else {
        ggplot2::theme(legend.position = legend_pos)
      }
    )
}

# simple plot of attribution category choices by valence, group, and timepoint
cattr_behav_plot <- function(
    data,
    plot_type = c("positive", "negative", "both"),
    mode = c("both", "randomised", "completer_status"),
    plot_style = c("bar", "slope"),
    error_type = c("sd", "se"),
    label_colour = NULL,
    bar_width = 0.65,
    err_width = 0.18,
    dodge_width = 0.75,
    alpha_bar = 0.8,
    line_width = 0.9,
    point_size = 3,
    shape_types = c(17, 15),
    time_sep = 0.16,
    ylim = NULL,
    fnt_sz = 1,
    axis_fnt_sc = 11,
    brew_col = "Archambault",
    col_nums = c(4, 1, 3),
    custom_pal = NULL) {

  plot_type  <- match.arg(plot_type)
  mode       <- match.arg(mode)
  plot_style <- match.arg(plot_style)
  error_type <- match.arg(error_type)
  group_var <- intersect(c("group", "intCond", "condition"), names(data))[1]
  has_non_completed <- "non_completed" %in% names(data)
  if (length(group_var) == 0 || is.na(group_var)) {
    stop("Data must include one grouping column: group, intCond, or condition")
  }

  to_plt <- data |>
    dplyr::mutate(
      group_raw = as.character(.data[[group_var]]),
      non_completed_plot = if (has_non_completed) as.logical(.data[["non_completed"]]) else FALSE,
      group_plot = dplyr::case_when(
        non_completed_plot ~ "non-completers",
        group_raw == "BA" ~ "behavioural activation",
        group_raw == "CR" ~ "cognitive restructuring",
        TRUE ~ NA_character_
      ),
      group_plot = factor(
        group_plot, levels = c("cognitive restructuring", "behavioural activation", "non-completers")
      ),
      timepoint_plot = factor(
        as.character(.data[["timepoint"]]), levels = c("0", "1"), labels = c("pre-intervention", "post-intervention")
      ),
      valence_plot = tolower(as.character(.data[["valence"]])),
      valence_facet = factor(
        valence_plot,
        levels = c("positive", "negative"),
        labels = c("positive events", "negative events")
      )
    )

  if (plot_type != "both") {
    to_plt <- to_plt |> dplyr::filter(valence_plot == plot_type)
  }

  if (mode == "randomised") {
    to_plt <- to_plt |>
      dplyr::mutate(
        group_plot = dplyr::case_when(
          group_raw == "BA" ~ "behavioural activation",
          group_raw == "CR" ~ "cognitive restructuring",
          TRUE ~ NA_character_
        ),
        group_plot = factor(group_plot, levels = c("cognitive restructuring", "behavioural activation"))
      )
  } else if (mode == "completer_status") {
    to_plt <- to_plt |>
      dplyr::filter(.data[["timepoint"]] == 0) |>
      dplyr::mutate(
        group_plot = dplyr::if_else(
          group_plot == "non-completers", "non-completers", "completers"
        ),
        group_plot = factor(group_plot, levels = c("completers", "non-completers"))
      )
  }

  full_pal <- MetBrewer::met.brewer(brew_col)
  if (!is.null(custom_pal)) {
    pal <- custom_pal
  } else if (is.null(col_nums)) {
    pal <- sample(full_pal, length(unique(to_plt$group_plot)))
  } else {
    pal <- full_pal[col_nums]
  }

  if (is.null(label_colour)) {
    if (plot_type == "positive") label_colour <- full_pal[2]
    else if (plot_type == "negative") label_colour <- full_pal[3]
    else label_colour <- c(full_pal[2], full_pal[3])
  }

  if (plot_type == "both") {
    if (length(label_colour) == 1) label_colour <- rep(label_colour, 2)
    if (length(label_colour) >= 2) {
      to_plt <- to_plt |>
        dplyr::mutate(
          valence_facet = factor(
            valence_plot,
            levels = c("positive", "negative"),
            labels = c(
              paste0("<span style='color:", label_colour[1], ";'>positive events</span>"),
              paste0("<span style='color:", label_colour[2], ";'>negative events</span>")
            )
          )
        )
    }
  }

  plot_df <- to_plt |>
    dplyr::transmute(
      id_plot = as.character(.data[["subID"]]),
      group_plot,
      timepoint_plot,
      valence_facet,
      internalGlobalChosen = as.numeric(.data[["internalGlobalChosen"]]),
      internalSpecificChosen = as.numeric(.data[["internalSpecificChosen"]]),
      externalGlobalChosen = as.numeric(.data[["externalGlobalChosen"]]),
      externalSpecificChosen = as.numeric(.data[["externalSpecificChosen"]])
    ) |>
    tidyr::pivot_longer(
      cols = c(
        internalGlobalChosen,
        internalSpecificChosen,
        externalGlobalChosen,
        externalSpecificChosen
      ),
      names_to = "category",
      values_to = "chosen"
    ) |>
    dplyr::mutate(
      category = factor(
        category,
        levels = c(
          "internalGlobalChosen",
          "internalSpecificChosen",
          "externalGlobalChosen",
          "externalSpecificChosen"
        ),
        labels = c(
          "internal-global",
          "internal-specific",
          "external-global",
          "external-specific"
        )
      )
    ) |>
    dplyr::filter(!is.na(id_plot), !is.na(group_plot), !is.na(chosen), !is.na(timepoint_plot), !is.na(category))

  subj_df <- plot_df |>
    dplyr::group_by(id_plot, group_plot, timepoint_plot, valence_facet, category) |>
    dplyr::summarise(chosen = mean(chosen), .groups = "drop")

  sum_df <- subj_df |>
    dplyr::group_by(group_plot, timepoint_plot, valence_facet, category) |>
    dplyr::summarise(
      chosen_mean = mean(chosen),
      chosen_sd = stats::sd(chosen),
      chosen_se = chosen_sd / sqrt(dplyr::n()),
      .groups = "drop"
    ) |>
    dplyr::mutate(err = if (error_type == "sd") chosen_sd else chosen_se)

  pd <- ggplot2::position_dodge(width = dodge_width)
  legend_key_h <- grid::unit(max(6, 10 * fnt_sz), "pt")
  legend_spacing_y <- grid::unit(max(0, 2 * fnt_sz), "pt")
  title_lab <- if (plot_type == "both") NULL else paste0(plot_type, " events")

  if (plot_style == "bar") {
    p <- ggplot2::ggplot(sum_df, ggplot2::aes(x = category, y = chosen_mean, fill = group_plot, colour = group_plot)) +
      ggplot2::geom_col(
        ggplot2::aes(group = group_plot),
        position = pd,
        width = bar_width,
        alpha = alpha_bar,
        linewidth = 0.3
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = pmax(0, chosen_mean - err),
          ymax = chosen_mean + err,
          group = group_plot
        ),
        linewidth = line_width * 0.6,
        width = err_width,
        alpha = 0.4,
        position = pd
      ) +
      ggplot2::geom_hline(
        yintercept = 0.25, linetype = "dashed", colour = "#cacaca", alpha = 0.8
      ) + {
      if (plot_type == "both") {
        if (mode == "completer_status") {
          ggplot2::facet_wrap(~valence_facet, nrow = 2)
        } else {
          ggplot2::facet_grid(rows = ggplot2::vars(valence_facet), cols = ggplot2::vars(timepoint_plot))
        }
      } else {
        if (mode != "completer_status") ggplot2::facet_wrap(~timepoint_plot, nrow = 1)
      }
    } +
      ggplot2::scale_colour_manual(name = "", values = pal) +
      ggplot2::scale_fill_manual(name = "", values = pal) +
      ggplot2::labs(
        title = title_lab,
        x = NULL,
        y = paste0("proportion of attributions (± ", ifelse(error_type == "sd", "s.d.", "s.e."), ")")
      )
  } else {
    p <- ggplot2::ggplot(
      sum_df,
      ggplot2::aes(
        x = timepoint_plot,
        y = chosen_mean,
        colour = group_plot,
        fill = group_plot,
        shape = group_plot,
        group = group_plot
      )
    ) +
      ggplot2::geom_line(
        linewidth = line_width,
        alpha = 0.7,
        position = pd
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = chosen_mean - err, ymax = chosen_mean + err),
        width = err_width,
        linewidth = line_width * 0.6,
        alpha = 0.6,
        position = pd
      ) +
      ggplot2::geom_point(
        size = point_size,
        stroke = 0.5,
        alpha = 0.95,
        position = pd
      ) +
      ggplot2::geom_hline(
        yintercept = 0.25, linetype = "dashed", colour = "#cacaca", alpha = 0.8
      ) +
      (if (plot_type == "both") {
        ggplot2::facet_grid(rows = ggplot2::vars(valence_facet), cols = ggplot2::vars(category))
      } else {
        ggplot2::facet_wrap(~category, nrow = 1)
      }) +
      ggplot2::scale_x_discrete(
        labels = c("pre-intervention" = "pre", "post-intervention" = "post")
      ) +
      ggplot2::scale_colour_manual(name = "", values = pal) +
      ggplot2::scale_fill_manual(name = "", values = pal) +
      ggplot2::scale_shape_manual(name = "", values = shape_types) +
      ggplot2::labs(
        title = title_lab,
        x = NULL,
        y = paste0("prop. chosen (± ", ifelse(error_type == "sd", "s.d.", "s.e."), ")")
      )
  }

  if (!is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(ylim = ylim)
  }

  p +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = 14 * fnt_sz) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.location = "plot",
      legend.justification = "center",
      legend.key.height = legend_key_h,
      legend.spacing.y = legend_spacing_y,
      legend.justification.bottom = "center",
      legend.background = ggplot2::element_rect(
        fill = "white", colour = "#cacaca", linewidth = 0.5
      ),
      legend.margin = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
      strip.background = ggplot2::element_blank(),
      strip.text.x = if (plot_type == "both" && plot_style == "bar" && mode == "completer_status") {
        ggtext::element_markdown(
          family = "Open Sans SemiBold", colour = "darkslategray4",
          size = 13 * fnt_sz, hjust = 0, margin = ggplot2::margin(b = 0.01, unit = "npc")
        )
      } else {
        ggtext::element_markdown(face = "italic", colour = "darkslategray4", size = 13 * fnt_sz)
      },
      strip.text.y = if (plot_type == "both" && plot_style == "bar" && mode != "completer_status") {
        ggtext::element_markdown(
          family = "Open Sans SemiBold", size = 12 * fnt_sz, hjust = 0,
          margin = ggplot2::margin(b = 0.01, unit = "npc")
        )
      } else {
        ggtext::element_markdown(face = "italic", size = 12 * fnt_sz)
      },
      plot.title = if (!plot_type == "both") {
        ggplot2::element_text(
          size = 14 * fnt_sz,
          face = "plain",
          family = "Open Sans SemiBold",
          colour = label_colour
        )
      } else {
        ggplot2::element_blank()
      },
      axis.title.y = ggplot2::element_text(size = axis_fnt_sc * fnt_sz),
      axis.title.x = ggtext::element_markdown(size = 12 * fnt_sz),
      axis.text.x = if (plot_style == "slope") {
        ggplot2::element_text(size = 14 * fnt_sz, angle = 0, hjust = 0.5, vjust = 0.5)
      } else {
        ggplot2::element_text(size = 10 * fnt_sz, angle = 20, hjust = 0.75, vjust = 0.85)
      },
      axis.text.y = ggplot2::element_text(size = 11 * fnt_sz),
      legend.text = ggplot2::element_text(size = 12 * fnt_sz),
      legend.box = "vertical"
    )
}

extract_indiv_pars <- function(draws, stan_ls, type = "effort", hdi = 0.95) {
  # drop columns of draws starting with a "."
  indiv_pars <- draws |> dplyr::select(-starts_with("."))
  nms <- strsplit(colnames(indiv_pars), "\\[|,|\\]|_p(?!\\w)", perl = TRUE)
  conds <- stan_ls$condition
  compl <- stan_ls$completed
  ret_df <- tibble::tibble(
    "variable" = sapply(nms, function(x) x[1]),
    "id" = sapply(nms, function(x) x[2]),
    "group" = sapply(
      nms,
      function(x) {
        if (x[2] == "") conds[as.integer(x[3])]
        else conds[as.integer(x[2])]
      }
    ),
    "completed" = sapply(
      nms,
      function(x) {
        if (x[2] == "") compl[as.integer(x[3])]
        else compl[as.integer(x[2])]
      }
    ),
    "timepoint" = sapply(nms, function(x) x[3]),
    "mean" = colMeans(indiv_pars),
    "sd" = apply(indiv_pars, 2, sd),
    "sem" = sd / sqrt(nrow(indiv_pars)),
    "hdi_low" = apply(indiv_pars, 2, function(x) bayestestR::hdi(x, ci = hdi)$CI_low),
    "hdi_high" = apply(indiv_pars, 2, function(x) bayestestR::hdi(x, ci = hdi)$CI_high),
  ) |>
    dplyr::mutate(
      id = ifelse(id == "", timepoint, id),
      hdi_width = hdi_high - hdi_low,
      timepoint = dplyr::case_when(
        grepl("delta", variable) ~ "t2-t1",
        timepoint == "1" ~ "t1",
        timepoint == "2" ~ "t2",
        .default = ""
      ),
      t = factor(
        timepoint,
        levels = c("t1", "t2", "t2-t1"),
        labels = c("t1", "t2", "change")
      )
    )

  if (type == "effort") {
    ret_df <- ret_df |>
      dplyr::mutate(
        group = ifelse(group == 1, "BA", "CR"),
        grp = factor(
          group,
          levels = c("CR", "BA"),
          labels = c("bsl", "int")
        )
      )
  } else if (type == "attribution") {
    ret_df <- ret_df |>
      dplyr::mutate(
        group = ifelse(group == 0, "BA", "CR"),
        grp = factor(
          group,
          levels = c("BA", "CR"),
          labels = c("bsl", "int")
        )
      )
  }
  ret_df
}

plot_indiv_pars <- function(pars,
                            var,
                            var_nm,
                            type = "point",
                            tpt = c("t1", "t2"),
                            qnr = NULL,
                            qnr_nm = NULL,
                            compl = FALSE,
                            fnt_sz = 1,
                            sz_scale = c(1, 4),
                            clr = "group",
                            clr_levs = c("BA", "CR"),
                            clr_labs = c("behavioural activation", "cognitive restructuring"),
                            brew_col = "Archambault",
                            col_nums = NULL,
                            custom_pal = NULL,
                            boxp_legend_pos = c(0.3, 0.05),
                            indiv_lines = FALSE) {
  type <- match.arg(type, c("point", "box"))

  if (compl) to_plt <- pars |> dplyr::filter(completed == 1)
  else to_plt <- pars
  if (is.null(col_nums)) pal <- sample(MetBrewer::met.brewer(brew_col), length(unique(to_plt[[clr]])))
  else if (!is.null(custom_pal)) pal <- custom_pal
  else pal <- MetBrewer::met.brewer(brew_col)[col_nums]
  to_plt[["grp"]] <- factor(to_plt[[clr]], levels = clr_levs, labels = clr_labs)
  if (length(var) > 1) {
    to_plt[["var"]] <- factor(to_plt[["variable"]], levels = var, labels = var_nm)
  }
  if (type == "point") {
    if (is.null(qnr)) {
      plt <- to_plt |>
        dplyr::filter(timepoint %in% tpt & variable == var) |>
        dplyr::select(id, grp, timepoint, mean, sd) |>
        tidyr::pivot_wider(names_from = timepoint, values_from = c(mean, sd)) |>
        ggplot2::ggplot(
          ggplot2::aes(x = mean_t1, y = mean_t2, color = grp, fill = grp)
        ) +
        ggplot2::geom_abline(slope = 1, linetype = "dashed", colour = "grey") +
        ggplot2::geom_point(size = 3, alpha = 0.5) +
        ggplot2::geom_errorbar(
          ggplot2::aes(xmin = mean_t1 - sd_t1, xmax = mean_t1 + sd_t1), alpha = .2, orientation = "y"
        ) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_t2 - sd_t2, ymax = mean_t2 + sd_t2), alpha = .2) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
        ggplot2::labs(
          x = paste0("mean (± s.d.) **", var_nm, "** timepoint 1"),
          y = paste0("mean (± s.d.) **", var_nm, "**<br>timepoint 2"),
        ) +
        cowplot::theme_minimal_grid(font_family = "Open Sans", font_size = 16 * fnt_sz)
    } else {
      qn <- rlang::sym(qnr)
      plt <- to_plt |>
        dplyr::filter(timepoint == "t2-t1" & grepl(var, variable)) |>
        ggplot2::ggplot(
          ggplot2::aes(x = !!qn, y = mean, color = grp, fill = grp)
        ) +
        ggplot2::geom_point(ggplot2::aes(size = 1 / sd), alpha = 0.5) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.35, width = 0.6
        ) +
        ggplot2::geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1.5) +
        ggplot2::guides(size = "none") +
        ggplot2::labs(
          x = paste0("change in **", qnr_nm, "**"),
          y = paste0("\u0394 **", var_nm, "**"),
        ) +
        ggplot2::scale_size(range = sz_scale) +
        cowplot::theme_minimal_grid(font_family = "Open Sans", font_size = 20 * fnt_sz)
    }
    plt +
      ggplot2::scale_colour_manual(name = NULL, values = pal) +
      ggplot2::scale_fill_manual(name = NULL, values = pal) +
      ggplot2::theme(
        axis.title.x = ggtext::element_markdown(lineheight = 1.2),
        axis.title.y = ggtext::element_markdown(lineheight = 1.2),
        legend.text = ggplot2::element_text(size = fnt_sz * 18),
        plot.title = ggtext::element_markdown(size = fnt_sz * 20, face = "plain", colour = "slateblue4"),
        plot.subtitle = ggtext::element_markdown(size = fnt_sz * 14, face = "plain", colour = "red4")
      )
  } else if (type == "box") {
    if (length(tpt) == 2) {
      x_aes <- "timepoint"
      to_plt <- to_plt |> dplyr::filter(timepoint %in% tpt & variable == var)
    } else {
      x_aes <- "var"
      to_plt <- to_plt |> dplyr::filter(timepoint %in% tpt)
    }

    plt_base <- to_plt |>
      ggplot2::ggplot(
        ggplot2::aes(x = .data[[x_aes]], y = mean, colour = grp, fill = grp)
      )

    if (length(tpt) == 2) {
      # For 2 timepoints, points are outside, boxplots inside
      set.seed(123) # for reproducible jitter
      to_plt <- to_plt |>
        dplyr::mutate(
          x_base = as.numeric(factor(grp)),
          x_jit = x_base + runif(dplyr::n(), -0.1, 0.1),
          x_nudge = dplyr::case_when(
            timepoint == tpt[1] ~ 0.2,
            timepoint == tpt[2] ~ -0.2
          ),
          x_final = as.numeric(factor(timepoint)) + x_nudge + (x_jit - x_base)
        )

      plt_base <- to_plt |>
        ggplot2::ggplot(
          ggplot2::aes(x = .data[[x_aes]], y = mean, colour = grp, fill = grp)
        ) +
        ggplot2::geom_point(
          ggplot2::aes(x = x_final, group = id, size = 1 / sd),
          alpha = 0.2
        ) +
        ggplot2::geom_boxplot(
          width = 0.15, alpha = 0.5, outlier.shape = NA,
          position = ggpp::position_dodgenudge(
            width = 0.2, x = c(-0.15, -0.15, 0.15, 0.15)
          )
        )
    } else {
      # For 1 timepoint, points are left, boxplots are right
      plt_base <- plt_base +
        ggplot2::geom_boxplot(
          width = 0.2, alpha = 0.5, outlier.shape = NA,
          position = ggpp::position_dodgenudge(width = 0.25, x = 0.15)
        ) +
        ggplot2::geom_point(
          ggplot2::aes(group = id, size = 1 / sd),
          alpha = 0.2,
          position = ggpp::position_jitternudge(
            width = 0.15, x = -0.2, seed = 123,
            nudge.from = "jittered"
          )
        )
    }

    if (length(tpt) == 2 && indiv_lines) {
      lines_df <- to_plt |>
        dplyr::select(id, grp, timepoint, x_final, mean) |>
        tidyr::pivot_wider(
          names_from = timepoint,
          values_from = c(x_final, mean)
        ) |>
        dplyr::rename(
          x1 = paste0("x_final_", tpt[1]),
          x2 = paste0("x_final_", tpt[2]),
          y1 = paste0("mean_", tpt[1]),
          y2 = paste0("mean_", tpt[2])
        ) |>
        tidyr::drop_na()

      plt_base <- plt_base +
        ggplot2::geom_segment(
          data = lines_df,
          ggplot2::aes(x = x1, xend = x2, y = y1, yend = y2, group = id),
          alpha = 0.05
        ) +
        ggplot2::stat_summary(
          fun = mean, geom = "point", ggplot2::aes(group = grp),
          size = 5, alpha = 0.9, shape = 23,
          position = ggpp::position_dodgenudge(
            width = 0.25, x = c(0.2, -0.2, 0.2, -0.2)
          )
        ) +
        ggplot2::stat_summary(
          fun = mean, geom = "line", ggplot2::aes(group = grp),
          linewidth = 1.5, alpha = 0.8, linetype = "dotted",
          position = ggpp::position_dodgenudge(
            width = 0.25, x = c(0.2, -0.2, 0.2, -0.2)
          )
        ) +
        ggplot2::stat_summary(
          fun.data = ggplot2::mean_se, geom = "errorbar",
          width = 0.1, linewidth = 1, alpha = 0.8,
          position = ggpp::position_dodgenudge(
            width = 0.25, x = c(0.2, 0.2, -0.2, -0.2)
          )
        )
    }

    if (length(tpt) == 2) {
      y_lab <- paste0("mean **", var_nm, "**")
      x_scale <- ggplot2::scale_x_discrete(
        labels = c("t1" = "pre-intervention", "t2" = "post-intervention"),
        name = NULL
      )
    } else {
      y_lab <- "mean timepoint 1 estimate"
      x_scale <- ggplot2::scale_x_discrete(name = NULL)
    }

    plt_base <- plt_base +
      ggplot2::ylab(y_lab) +
      x_scale +
      cowplot::theme_half_open(font_family = "Open Sans", font_size = 16 * fnt_sz) +
      ggplot2::scale_colour_manual(name = NULL, labels = clr_labs, values = pal) +
      ggplot2::scale_fill_manual(name = NULL, labels = clr_labs, values = pal) +
      ggplot2::scale_size(range = sz_scale) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(position = "inside", nrow = 1),
        colour = ggplot2::guide_legend(position = "inside", nrow = 1),
        size = "none"
      ) +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = fnt_sz * 18),
        legend.position.inside = boxp_legend_pos,
        legend.key.spacing.x = grid::unit(16 * fnt_sz, "pt"),
        axis.title.x = ggtext::element_markdown(),
        axis.title.y = ggtext::element_markdown()
      )

    if (length(tpt) == 1) {
      plt_base <- plt_base + ggplot2::theme(
        axis.title.x = ggtext::element_markdown(size = fnt_sz * 18),
        axis.title.y = ggtext::element_markdown(size = fnt_sz * 18),
        axis.text.x = ggtext::element_markdown(size = fnt_sz * 16),
        axis.text.y = ggplot2::element_text(size = fnt_sz * 16)
      )
    }
    plt_base
  }
}

# Shared rendering helpers: both plot_reff_change() and plot_sympt_change() call these.
# pop_df must have columns: group (factor), tp ("pre"/"post"), m, lo95, hi95, lo99, hi99, ybase.
plot_group_levels <- function(pop_df,
                              xlim = NULL,
                              x_lab,
                              grp_cols,
                              suppress_labels = FALSE,
                              flip_order = FALSE,
                              fnt_sz = 1,
                              subtitle = "individual-level estimates\n○ pre- ● post-intervention") {
  if (flip_order) {
    nudge_sign   <- pop_df$ybase - as.numeric(pop_df$group)
    pop_df$group <- factor(pop_df$group, levels = rev(levels(pop_df$group)))
    pop_df$ybase <- as.numeric(pop_df$group) + nudge_sign
  }
  grp_levels <- levels(pop_df$group)
  pre  <- pop_df[pop_df$tp == "pre",  ]
  post <- pop_df[pop_df$tp == "post", ]
  plt <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = pop_df,
      ggplot2::aes(x = m, y = ybase, group = group, colour = group),
      linewidth = 0.7, alpha = 0.5, linetype = "dotted"
    ) +
    ggplot2::geom_linerange(
      data = pop_df,
      ggplot2::aes(xmin = lo99, xmax = hi99, y = ybase, colour = group),
      linewidth = 0.9
    ) +
    ggplot2::geom_linerange(
      data = pop_df,
      ggplot2::aes(xmin = lo95, xmax = hi95, y = ybase, colour = group),
      linewidth = 2.0
    ) +
    ggplot2::geom_point(
      data = pre,  ggplot2::aes(x = m, y = ybase, colour = group),
      shape = 21, fill = "white", size = 3.4, stroke = 1.2
    ) +
    ggplot2::geom_point(
      data = post, ggplot2::aes(x = m, y = ybase, colour = group, fill = group),
      shape = 21, size = 3.4, stroke = 1.2
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(grp_levels), labels = gsub(" ", "\n", grp_levels),
      limits = c(0.6, length(grp_levels) + 0.4)
    ) +
    ggplot2::scale_colour_manual(values = grp_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = grp_cols, guide = "none") +
    ggplot2::labs(x = x_lab, y = NULL, subtitle = subtitle) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = 13 * fnt_sz) +
    cowplot::background_grid(major = "x", minor = "none") +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 12 * fnt_sz, colour = "grey30"),
      axis.text.y   = ggplot2::element_text(size = 10 * fnt_sz)
    )
  if (!is.null(xlim)) plt <- plt + ggplot2::coord_cartesian(xlim = xlim)
  if (suppress_labels) {
    plt <- plt + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }
  plt
}

# chg_df must have columns: level (factor, "difference" < "behavioural activation" < "cognitive restructuring"), value.
plot_change_slabs <- function(chg_df,
                              xlim = NULL,
                              x_lab,
                              grp_cols,
                              grp_shapes,
                              diff_comp = "(goal- − restr.)",
                              flip_order = FALSE,
                              suppress_labels = FALSE,
                              fnt_sz = 1,
                              subtitle = "group-level mean pre → post change") {
  if (flip_order) {
    lvls  <- levels(chg_df$level)
    other <- rev(lvls[lvls != "difference"])
    chg_df$level <- factor(chg_df$level, levels = c("difference", other))
  }
  plt <- chg_df |>
    ggplot2::ggplot(
      ggplot2::aes(x = value, y = level, fill = level, colour = level, shape = level)
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1.5, linetype = "dotted", colour = "grey70", linewidth = 2) +
    ggdist::stat_pointinterval(
      .width = c(0.95, 0.99), point_interval = "mean_hdci", interval_size_range = c(1, 3), fatten_point = 2
    ) +
    ggplot2::scale_colour_manual(values = grp_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = grp_cols, guide = "none") +
    ggplot2::scale_shape_manual(values = grp_shapes, guide = "none") +
    ggplot2::scale_y_discrete(labels = c(
      "difference"              = diff_comp,
      "behavioural activation"  = "behavioural\nactivation",
      "cognitive restructuring" = "cognitive\nrestructuring"
    )) +
    ggplot2::labs(x = x_lab, y = NULL, subtitle = subtitle) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = 14 * fnt_sz) +
    cowplot::background_grid(major = "x", minor = "none") +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 12 * fnt_sz, colour = "grey30"),
      axis.text.y   = ggplot2::element_text(size = 10 * fnt_sz),
      # markdown: callers subscript the symptom scale, e.g. "Δ PHQ-9<sub>std</sub>". Plain labels
      # pass through unchanged (no measure or parameter name here contains markdown syntax).
      axis.title.x   = ggtext::element_markdown(size = 12 * fnt_sz)
    )
  if (!is.null(xlim)) plt <- plt + ggplot2::coord_cartesian(xlim = xlim)
  if (suppress_labels) {
    plt <- plt + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }
  plt
}

# Composite reward-effort results plot for one parameter:
#   top    -- group-level population mean per group x timepoint (mu_{param}[t],
#             + {param}_int for BA at t2), pre (open) -> post (filled), dodged,
#             with thick 95% / thin 99% HDI bars. Same raw (pre-Phi_approx)
#             latent scale as the bottom panel and as the reported delta_{param}/
#             {param}_int statistics -- NOT the Phi_approx-transformed, bounded
#             natural scale that individual-level rewSens/effSens draws are on.
#   bottom -- group-level pre->post change: cogn. restructuring (delta), behav.
#             activation (delta + interaction), and the difference (interaction),
#             as slab+interval (95%/99% HDI). Bottom is exactly top(post) minus
#             top(pre) per arm, by construction.
# `draws` are the saved group-level draws (mu_*, delta_*, *_int).
# Returns list(levels_panel, change_panel).
plot_reff_change <- function(draws,
                             param,
                             param_nm,
                             flip_order = TRUE, # BA vs CR default
                             diff_lab = "difference\n(BA − CR)",
                             level_xlim = NULL,
                             change_xlim = NULL,
                             suppress_labels = FALSE,
                             brew_col = "Archambault",
                             col_nums = c(1, 4),
                             diff_col = "#bb9bb4",
                             fnt_sz = 1) {
  pal <- MetBrewer::met.brewer(brew_col)[col_nums]
  grp_cols <- c(
    "behavioural activation" = pal[1], "cognitive restructuring" = pal[2], "difference" = diff_col
  )
  if (flip_order) {
    grp_shapes = c(16, 17, 15)
  } else {
    grp_shapes = c(16, 15, 17)
  }
  nudge <- 0.16

  # --- top: group-level means (mu_{param}[t], + {param}_int for BA @ t2) with 95%/99% HDI ---
  # condition 1 = BA (behavioural activation), 0 = CR (cognitive restructuring, reference)
  mu1 <- draws[[paste0("mu_", param, "[1]")]]  # pre, shared by both arms (no arm-specific t1 param)
  mu2 <- draws[[paste0("mu_", param, "[2]")]]  # post, CR (reference)
  int <- draws[[paste0(param, "_int")]]        # additional BA effect at post

  pop_rows <- list()
  for (g in c(1, 0)) {
    for (t in 1:2) {
      v <- if (t == 1) mu1 else if (g == 1) mu2 + int else mu2
      h95 <- bayestestR::hdi(v, ci = 0.95)
      h99 <- bayestestR::hdi(v, ci = 0.99)
      pop_rows[[length(pop_rows) + 1]] <- data.frame(
        group = ifelse(g == 1, "behavioural activation", "cognitive restructuring"),
        tp = ifelse(t == 1, "pre", "post"),
        m = median(v),
        lo95 = h95$CI_low, hi95 = h95$CI_high,
        lo99 = h99$CI_low, hi99 = h99$CI_high
      )
    }
  }
  pop <- do.call(rbind, pop_rows)
  pop$group <- factor(pop$group, levels = c("behavioural activation", "cognitive restructuring"))
  pop$ybase <- as.numeric(pop$group) + ifelse(pop$tp == "pre", nudge, -nudge)

  # --- bottom: group-level change (slab + 95%/99% HDI), no individuals ---
  # delta_{param} = CR change; param_int = additional BA change vs. CR
  chg <- dplyr::bind_rows(
    data.frame(level = "cognitive restructuring", value = draws[[paste0("delta_", param)]]),
    data.frame(
      level = "behavioural activation",
      value = draws[[paste0("delta_", param)]] + draws[[paste0(param, "_int")]]
    ),
    data.frame(level = "difference", value = draws[[paste0(param, "_int")]])
  )
  chg$level <- factor(
    chg$level, levels = c("difference", "behavioural activation", "cognitive restructuring")
  )

  list(
    levels_panel = plot_group_levels(
      pop, xlim = level_xlim,
      x_lab      = param_nm,
      grp_cols   = grp_cols, fnt_sz = fnt_sz,
      suppress_labels = suppress_labels,
      subtitle   = "group-level estimates\n○ pre- ● post-intervention"
    ),
    change_panel = plot_change_slabs(
      chg, xlim = change_xlim, diff_comp = diff_lab,
      flip_order = flip_order,
      grp_shapes = grp_shapes,
      x_lab      = paste0("change in ", param_nm),
      grp_cols   = grp_cols, fnt_sz = fnt_sz,
      suppress_labels = suppress_labels
    )
  )
}

# Composite PHQ-9 symptom-association results plot for one parameter.
# Mirrors plot_reff_change() structure exactly, using pre-computed posterior
# draws from the 3-way interaction symptom model:
#   top    -- PHQ-9 -> parameter association at t1 and t2 per group (dodged),
#             with thick 95% / thin 99% HDI bars.
#   bottom -- change in association: restructuring, goal-setting, and the
#             group difference (3-way interaction), as slab+interval.
# `draws`  = data frame with columns named "{prefix}_t1_group_ba" etc.
# `prefix` = e.g. "phq_rewSens" or "phq_effSens".
# Returns list(levels_panel, change_panel).
# plot_sympt_change <- function(draws,
#                               prefix,
#                               param_nm,
#                               diff_lab = "difference\n(BA − CR)",
#                               level_xlim = NULL,
#                               change_xlim = NULL,
#                               flip_order = FALSE,
#                               suppress_labels = FALSE,
#                               brew_col = "Archambault",
#                               col_nums = c(1, 4),
#                               diff_col = "#bb9bb4",
#                               fnt_sz = 1) {
#   pal <- MetBrewer::met.brewer(brew_col)[col_nums]
#   grp_cols <- c(
#     "behavioural activation" = pal[1], "cognitive restructuring" = pal[2], "difference" = diff_col
#   )
#   nudge <- 0.16

#   # --- top: PHQ-9 association at t1 and t2 per group ---
#   combos <- list(
#     list(group = "behavioural activation",  tp = "pre",  col = paste0(prefix, "_t1_group_ba")),
#     list(group = "cognitive restructuring", tp = "pre",  col = paste0(prefix, "_t1_group_cr")),
#     list(group = "behavioural activation",  tp = "post", col = paste0(prefix, "_t2_group_ba")),
#     list(group = "cognitive restructuring", tp = "post", col = paste0(prefix, "_t2_group_cr"))
#   )
#   pop_rows <- lapply(combos, function(x) {
#     v   <- draws[[x$col]]
#     h95 <- bayestestR::hdi(v, ci = 0.95)
#     h99 <- bayestestR::hdi(v, ci = 0.99)
#     data.frame(
#       group = x$group, tp = x$tp, m = median(v),
#       lo95 = h95$CI_low, hi95 = h95$CI_high,
#       lo99 = h99$CI_low, hi99 = h99$CI_high
#     )
#   })
#   pop <- do.call(rbind, pop_rows)
#   pop$group <- factor(pop$group, levels = c("behavioural activation", "cognitive restructuring"))
#   pop$ybase <- as.numeric(pop$group) + ifelse(pop$tp == "pre", nudge, -nudge)

#   # --- bottom: change in association per group + group difference ---
#   chg <- dplyr::bind_rows(
#     data.frame(level = "cognitive restructuring", value = draws[[paste0(prefix, "_t2_change_cr")]]),
#     data.frame(level = "behavioural activation",  value = draws[[paste0(prefix, "_t2_change_ba")]]),
#     data.frame(level = "difference",              value = draws[[paste0(prefix, "_t2_change_diff")]])
#   )
#   chg$level <- factor(
#     chg$level, levels = c("difference", "behavioural activation", "cognitive restructuring")
#   )

#   list(
#     levels_panel = plot_group_levels(
#       pop, xlim = level_xlim,
#       x_lab      = paste0("PHQ-9<sub>std</sub> → ", param_nm),
#       grp_cols   = grp_cols, fnt_sz = fnt_sz, flip_order = flip_order,
#       suppress_labels = suppress_labels,
#       subtitle   = paste0(
#         "symptom-parameter association\n○ pre-  ● post-intervention"
#       )
#     ) +
#       # dashed vertical line at 0 for reference
#       ggplot2::geom_vline(
#         xintercept = 0, linetype = "dashed", colour = "grey60"
#       ),
#     change_panel = plot_change_slabs(
#       chg, xlim = change_xlim,
#       x_lab      = paste0("Δ PHQ-9<sub>std</sub> → ", param_nm),
#       grp_cols   = grp_cols, fnt_sz = fnt_sz, diff_comp = diff_lab,
#       flip_order = flip_order, suppress_labels = suppress_labels,
#       subtitle   = paste0("pre → post change in association")
#     )
#   )
# }

# Raw-data companion to plot_sympt_decomp()'s within_panel: one point per
# participant (posterior-mean Δparameter from the decomp model, the same
# quantity averaged into a_cr/a_ba in the mediation figure) against observed
# Δsymptom, coloured by group with per-group lm fits.
# `df` needs columns: delta_par, delta_symptom, group ("behavioural
# activation"/"cognitive restructuring").
plot_indiv_change_scatter <- function(df,
                                      param_nm,
                                      symptom_lab = "Δ PHQ-9<sub>std</sub>",
                                      xlim = NULL,
                                      ylim = NULL,
                                      x_n_breaks = 4,
                                      suppress_labels = FALSE,
                                      flip_order = FALSE,
                                      brew_col = "Archambault",
                                      col_nums = c(1, 4),
                                      fnt_sz = 1) {
  pal <- MetBrewer::met.brewer(brew_col)[col_nums]
  grp_cols <- c("behavioural activation" = pal[1], "cognitive restructuring" = pal[2])
  plt <- df |>
    ggplot2::ggplot(
      ggplot2::aes(x = delta_par, y = delta_symptom, colour = group, shape = group, fill = group)
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(size = 1 / delta_par_sd), alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.15, linewidth = 1.2) +
    ggplot2::scale_colour_manual(values = grp_cols, guide = "none") +
    ggplot2::scale_fill_manual(values = grp_cols, guide = "none") +
    ggplot2::scale_shape_manual(values = c(15, 17), guide = "none") +
    ggplot2::scale_size_continuous(range = c(1, 3)) +
    ggplot2::scale_x_continuous(n.breaks = x_n_breaks) +
    ggplot2::guides(size = "none") +
    ggplot2::labs(x = paste0("Δ ", param_nm), y = symptom_lab) +
    cowplot::theme_minimal_hgrid(font_family = "Open Sans", font_size = 14 * fnt_sz) +
    # markdown so `symptom_lab` can subscript the scaling (PHQ-9<sub>std</sub>). Both titles are
    # merged into the completed theme's existing elements, so size and the y title's 90-degree
    # angle carry over; only the renderer changes.
    ggplot2::theme(
      axis.title.x = ggtext::element_markdown(), axis.title.y = ggtext::element_markdown()
    )
  if (!is.null(xlim) || !is.null(ylim)) {
    plt <- plt + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }
  if (suppress_labels) {
    plt <- plt + ggplot2::theme(
      axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank()
    )
  }
  plt
}

# Companion to plot_sympt_change() for the within/between decomposed symptom
# models (…-symptoms-{rew,eff}-decomp.stan). Those models don't estimate a single
# t2 "level" association (the level is split into baseline + change), so the
# pre/post levels panel isn't well-defined. Returns three panels:
#   scatter_panel -- raw within-person change: one point per participant
#                    (posterior-mean Δparameter) vs. observed Δsymptom
#                    (PHQ-9<sub>std</sub>), by group with per-group lm fits -- the
#                    raw-data view of within_panel below it.
#   between_panel -- between-person (baseline) association = beta_phq_{par}B(+_group)
#   within_panel  -- within-person  (change)   association = beta_phq_{par}W(+_group)
# These are exactly the reported quantities -- what you plot is what you report.
# `draws` must contain the linear-combination columns built by the caller:
#   {prefix}_t1_group_cr / _ba / _t1_diff         -> between (baseline) association
#   {prefix}_t2_change_cr / _ba / _t2_change_diff -> within (change) association
# ...as well as the indexed per-participant `delta_{param}_p[1..nPpts]`
# columns from the same decomp model fit (`param` inferred by stripping
# "phq_" off `prefix`, e.g. prefix = "phq_rewSens" -> looks for
# delta_rewSens_p). `stan_ls` supplies `condition`/`completed` (for
# extract_indiv_pars()) and the z-scored `PHQ9` matrix -- both `stan_ls` and
# the decomp model's `p` index come from the same prep_data() call, so
# `delta_{param}_p[p]` and `PHQ9[p, ]` always refer to the same participant
# and no id/subID join is needed. `type` is passed straight through to
# extract_indiv_pars() -- "effort" (default) for reward-effort's
# condition==1-is-BA encoding, "attribution" for causal-attribution's
# condition==1-is-CR encoding (see prep_data()); get this wrong and BA/CR
# get silently swapped in the scatter panel.
plot_sympt_decomp <- function(draws,
                              stan_ls,
                              prefix,
                              param_nm,
                              type = "effort",
                              diff_lab = "difference\n(BA − CR)",
                              between_xlim = NULL,
                              within_xlim = NULL,
                              scatter_xlim = NULL,
                              scatter_ylim = NULL,
                              flip_order = FALSE,
                              suppress_labels = FALSE,
                              brew_col = "Archambault",
                              col_nums = c(1, 4),
                              diff_col = "#bb9bb4",
                              fnt_sz = 1) {
  pal <- MetBrewer::met.brewer(brew_col)[col_nums]
  grp_cols <- c(
    "behavioural activation" = pal[1], "cognitive restructuring" = pal[2], "difference" = diff_col
  )
  if (flip_order) {
    grp_shapes <- c(16, 17, 15)
  } else {
    grp_shapes <- c(16, 15, 17)
  }

  mediator_par <- sub("^phq_", "", prefix)
  required_cols <- paste0(prefix, c(
    "_t1_group_cr", "_t1_group_ba", "_t1_diff",
    "_t2_change_cr", "_t2_change_ba", "_t2_change_diff"
  ))
  missing_cols <- setdiff(required_cols, names(draws))
  if (length(missing_cols) > 0) {
    stop(
      "plot_sympt_decomp(): `draws` is missing expected columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  delta_cols <- grep(paste0("^delta_", mediator_par, "_p\\["), names(draws), value = TRUE)
  if (length(delta_cols) == 0) {
    stop(
      "plot_sympt_decomp(): `draws` has no columns matching 'delta_", mediator_par, "_p[...]' -- ",
      "make sure 'delta_", mediator_par, "_p' was included in the variables pulled from the fit."
    )
  }

  slab_df <- function(cr_col, ba_col, diff_col_nm) {
    d <- dplyr::bind_rows(
      data.frame(level = "cognitive restructuring", value = draws[[cr_col]]),
      data.frame(level = "behavioural activation",  value = draws[[ba_col]]),
      data.frame(level = "difference",              value = draws[[diff_col_nm]])
    )
    d$level <- factor(
      d$level, levels = c("difference", "behavioural activation", "cognitive restructuring")
    )
    d
  }

  between <- slab_df(
    paste0(prefix, "_t1_group_cr"), paste0(prefix, "_t1_group_ba"), paste0(prefix, "_t1_diff")
  )
  within <- slab_df(
    paste0(prefix, "_t2_change_cr"), paste0(prefix, "_t2_change_ba"), paste0(prefix, "_t2_change_diff")
  )

  # --- raw within-person scatter: Δparameter (posterior mean/participant) vs. ΔPHQ-9_std ---
  indiv_pars <- extract_indiv_pars(
    draws |> dplyr::select(tidyselect::all_of(delta_cols)), stan_ls, type = type
  ) |>
    dplyr::filter(variable == paste0("delta_", mediator_par)) |>
    dplyr::mutate(
      group = ifelse(group == "BA", "behavioural activation", "cognitive restructuring")
    )

  phq_delta <- tibble::tibble(
    id = as.character(seq_len(nrow(stan_ls$PHQ9))),
    delta_symptom = ifelse(
      stan_ls$PHQ9[, 1] == -999 | stan_ls$PHQ9[, 2] == -999,
      NA, stan_ls$PHQ9[, 2] - stan_ls$PHQ9[, 1]
    )
  )

  scatter_df <- indiv_pars |>
    dplyr::select(id, group, delta_par = mean, delta_par_sd = sd) |>
    dplyr::left_join(phq_delta, by = "id") |>
    tidyr::drop_na(delta_symptom)

  list(
    scatter_panel = plot_indiv_change_scatter(
      scatter_df, param_nm = param_nm, xlim = scatter_xlim, ylim = scatter_ylim,
      suppress_labels = suppress_labels, brew_col = brew_col, col_nums = col_nums, fnt_sz = fnt_sz
    ),
    between_panel = plot_change_slabs(
      between, xlim = between_xlim,
      x_lab      = paste0(param_nm, " → PHQ-9<sub>std</sub> (baseline)"),
      grp_cols   = grp_cols, fnt_sz = fnt_sz, diff_comp = diff_lab,
      grp_shapes = grp_shapes,
      flip_order = flip_order, suppress_labels = suppress_labels,
      subtitle   = "between-person (baseline) association"
    ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60"),
    within_panel = plot_change_slabs(
      within, xlim = within_xlim,
      x_lab      = paste0("Δ ", param_nm, " → Δ PHQ-9<sub>std</sub>"),
      grp_cols   = grp_cols, fnt_sz = fnt_sz, diff_comp = diff_lab,
      grp_shapes = grp_shapes,
      flip_order = flip_order, suppress_labels = suppress_labels,
      subtitle = NULL
    )
  )
}

# Small inset interval plot for one mediation path (CR on top, BA below by default --
# flip_order = TRUE puts BA on top instead), used by plot_joint_mediation() -- styled like
# the indirect-effect companion panel so the whole figure reads in one visual language.
# `flip_order` follows the same convention as plot_change_slabs()/plot_sympt_decomp()
# elsewhere: unflipped = CR on top (the causal-attribution "arm of interest" default),
# flip_order = TRUE = BA on top (reward-effort's "arm of interest" convention).
# `cr = NULL` switches to SINGLE-SERIES mode: one row, labelled `single_lab` and drawn in the
# neutral difference colour. Needed for the ANCOVA symptom submodels, whose direct (c') effect is a
# single baseline-adjusted arm contrast rather than a per-arm pair -- an ANCOVA's outcome is the t2
# LEVEL, so there is no within-arm change term to report per arm.
mini_interval_plot <- function(ba, cr = NULL, title, grp_cols, fnt_sz = 1, flip_order = FALSE,
                               single_lab = "BA \u2212 CR", single_col = "#bb9bb4") {
  if (is.null(cr)) {
    d <- data.frame(level = factor(single_lab, levels = single_lab), value = ba)
    pal <- stats::setNames(single_col, single_lab)
    shp <- stats::setNames(18, single_lab)
    y_labs <- stats::setNames(single_lab, single_lab)
  } else {
    d <- dplyr::bind_rows(
      data.frame(level = "behavioural activation",  value = ba),
      data.frame(level = "cognitive restructuring", value = cr)
    )
    lvls <- if (flip_order) {
      c("cognitive restructuring", "behavioural activation")
    } else {
      c("behavioural activation", "cognitive restructuring")
    }
    d$level <- factor(d$level, levels = lvls)
    pal <- grp_cols
    shp <- stats::setNames(if (flip_order) c(17, 15) else c(15, 17), lvls)
    y_labs <- c("behavioural activation" = "BA", "cognitive restructuring" = "CR")
  }
  ggplot2::ggplot(d, ggplot2::aes(x = value, y = level, colour = level, shape = level)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    ggdist::stat_pointinterval(
      .width = c(0.95, 0.99), point_interval = "mean_hdci", interval_size_range = c(0.75, 2)
    ) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_shape_manual(values = shp, guide = "none") +
    ggplot2::scale_x_continuous(n.breaks = 3) +
    ggplot2::scale_y_discrete(labels = y_labs) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = 12 * fnt_sz) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9.5 * fnt_sz, colour = "grey30"),
      plot.title = ggplot2::element_text(
        size = 11 * fnt_sz, face = "italic", colour = "slategray", hjust = 0.5
      )
    )
}

# Combined indirect-effect panel for a two-mediator joint model: two dodged
# pointranges per row (mediator 1 above, mediator 2 below -- matching the path
# diagram, where mediator 1's path arcs above the direct line and mediator 2's
# mirrors below). Colour (`mediator_cols`) also tints the M1/M2 boxes in
# plot_joint_mediation()'s diagram, so the diagram doubles as a visual key --
# but this panel additionally draws its own bottom, centre-justified colour
# legend (mediator1_nm/mediator2_nm as labels, falling back to the raw suffix
# if a display name isn't supplied) so the panel is self-explanatory on its
# own, styled like the behavioural plots' legends (white background, thin grey
# border, centred -- see reff_plot()/cattr_behav_plot()/plot_qnr_measure()).
# x-axis is deliberately generic ("Δ parameters", not naming one mediator)
# since colour + legend + diagram already say which is which.
#   mediator1_suffix/mediator2_suffix: the `_{suffix}_cr`/`_{suffix}_ba` column
#     suffixes in `draws_med` (e.g. "rew"/"eff", or "int_neg"/"glob_neg").
#   mediator1_nm/mediator2_nm: display names for the legend (e.g. "reward
#     sensitivity"); NULL (default) falls back to the raw suffix.
#   mediator_cols: named vector, names == c(mediator1_suffix, mediator2_suffix).
plot_joint_indirect <- function(draws_med,
                                mediator1_suffix, mediator2_suffix,
                                mediator_cols,
                                mediator1_nm = NULL, mediator2_nm = NULL,
                                diff_lab = "difference\n(BA − CR)",
                                flip_order = FALSE,
                                fnt_sz = 1) {
  if (!all(c(mediator1_suffix, mediator2_suffix) %in% names(mediator_cols))) {
    stop(
      "plot_joint_indirect(): `mediator_cols` must be named with `mediator1_suffix`/",
      "`mediator2_suffix` (\"", mediator1_suffix, "\"/\"", mediator2_suffix, "\") -- got names: ",
      paste(names(mediator_cols), collapse = ", ")
    )
  }
  mediator1_lab <- if (is.null(mediator1_nm)) mediator1_suffix else mediator1_nm
  mediator2_lab <- if (is.null(mediator2_nm)) mediator2_suffix else mediator2_nm
  mediator_labs <- stats::setNames(
    c(mediator1_lab, mediator2_lab), c(mediator1_suffix, mediator2_suffix)
  )
  row_vals <- function(mediator, level, value) {
    data.frame(mediator = mediator, level = level, value = value)
  }
  ind_df <- dplyr::bind_rows(
    row_vals(mediator1_suffix, "behavioural activation",  draws_med[[paste0("indirect_", mediator1_suffix, "_ba")]]),
    row_vals(mediator1_suffix, "cognitive restructuring", draws_med[[paste0("indirect_", mediator1_suffix, "_cr")]]),
    row_vals(mediator1_suffix, "difference",              draws_med[[paste0("index_mod_med_", mediator1_suffix)]]),
    row_vals(mediator2_suffix, "behavioural activation",  draws_med[[paste0("indirect_", mediator2_suffix, "_ba")]]),
    row_vals(mediator2_suffix, "cognitive restructuring", draws_med[[paste0("indirect_", mediator2_suffix, "_cr")]]),
    row_vals(mediator2_suffix, "difference",              draws_med[[paste0("index_mod_med_", mediator2_suffix)]])
  )
  # unflipped = CR on top (causal-attribution default), flip_order = TRUE = BA on top
  # (reward-effort default) -- same convention as mini_interval_plot()/plot_change_slabs()
  other <- if (flip_order) {
    c("cognitive restructuring", "behavioural activation")
  } else {
    c("behavioural activation", "cognitive restructuring")
  }

  ind_df$level <- factor(ind_df$level, levels = c("difference", other))
  nudge <- 0.16
  # mediator 1 sits above mediator 2 within each row, mirroring the diagram
  # (mediator 1's path arcs above the direct line, mediator 2's mirrors below)
  ind_df$ybase <- as.numeric(ind_df$level) + ifelse(ind_df$mediator == mediator1_suffix, nudge, -nudge)

  ind_df |>
    ggplot2::ggplot(
      ggplot2::aes(x = value, y = ybase, colour = mediator, group = interaction(mediator, level))
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1.5, linetype = "dotted", colour = "grey70", linewidth = 2) +
    ggdist::stat_pointinterval(
      .width = c(0.95, 0.99), point_interval = "mean_hdci", interval_size_range = c(1, 3), fatten_point = 2
    ) +
    ggplot2::scale_colour_manual(
      values = mediator_cols, labels = mediator_labs,
      breaks = c(mediator1_suffix, mediator2_suffix), name = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = 1:3,
      labels = c(diff_lab, gsub(" ", "\n", other)),
      limits = c(0.6, 3.4)
    ) +
    ggplot2::labs(
      x = "indirect path (a × b): Δ symptoms via Δ parameters",
      y = NULL
    ) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = 16 * fnt_sz) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12 * fnt_sz),
      axis.text.y = ggplot2::element_text(size = 10 * fnt_sz),
      legend.position = "bottom",
      legend.justification.bottom = "center",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11 * fnt_sz),
      legend.background = ggplot2::element_rect(
        fill = "white", colour = "#cacaca", linewidth = 0.5
      ),
      legend.box.spacing = grid::unit(16, "pt")
    )
}

# Two-mediator (parallel) path diagram for a joint mediation model (…-symptoms-
# joint-decomp.stan). Both mediators enter the same PHQ9 regression there, so
# each b-path (and indirect effect) is already adjusted for the other -- this
# just gives that joint model the two-mediator path diagram it needs: mediator
# 1 drawn as an arc *above* the direct (c') path, mediator 2 mirrored *below*
# it, both converging on the same symptom-change box. The direct (c') path is
# nudged below the X/Y boxes' centreline to make room for its own inset (titled
# "direct effect"; note what that means here -- c' is what is LEFT of the module ->
# symptom association after adjusting for both mediators drawn, and for a model with
# more than 2 mediators overall, potentially for more than this particular diagram
# shows). The M1/M2 boxes are filled with their own colour
# (`mediator_cols`, at `box_fill_alpha`) so they double as the key for the combined
# indirect-effect panel below (plot_joint_indirect()), which also draws its own
# bottom, centre-justified colour legend (mediator1_nm/mediator2_nm as labels) so
# that panel reads on its own rather than relying solely on the diagram's box tint;
# box-label text colour is set via `med_col_text` (one colour for both M1/M2 boxes).
# Generic over any 2-mediator pair sharing the
# `{quantity}_{suffix}_{cr,ba}` / `index_mod_med_{suffix}` naming convention
# (reward-effort's `_rew_`/`_eff_`, or e.g. causal attribution's
# `_int_neg_`/`_glob_neg_`) -- no reward/effort-specific naming is hardcoded.
#   mediator1_suffix/mediator2_suffix: the `_{suffix}_cr`/`_{suffix}_ba` column
#     suffixes in `draws_med` (e.g. "rew"/"eff").
#   mediator1_nm/mediator2_nm: display names (e.g. "reward sensitivity").
#   mediator_cols: named vector, names == c(mediator1_suffix, mediator2_suffix).
#   box_fill_alpha: opacity of the M1/M2 box fills (1 = solid mediator colour;
#     lower = a lighter wash of it). `med_col_text` sets the box-label text colour
#     (default "white"; pick a dark colour instead if using a light `box_fill_alpha`).
#   flip_order: FALSE (default) draws CR on top/BA below in every path inset and
#     the indirect panel (causal-attribution's "arm of interest" convention);
#     TRUE flips to BA on top (reward-effort's convention) -- same meaning as
#     flip_order elsewhere (plot_change_slabs()/plot_sympt_decomp()).
#   title/title_col: optional bold title drawn in the diagram's own whitespace
#     above the M1 box (e.g. "positive events"/"negative events" for the
#     causal-attribution valence split); NULL (default) draws nothing.
#   draws_med: df with a_{suffix}_cr/ba, b_{suffix}_cr/ba, direct_cr/ba,
#              indirect_{suffix}_cr/ba, index_mod_med_{suffix} for both
#              mediators (exactly the columns pulled from the joint model fit).
# Returns list(diagram, a1_path, b1_path, a2_path, b2_path, direct_path, indirect) of
# ggplots. The main script insets the five path plots onto `diagram` (by their arrows)
# and stacks `indirect` alongside it.
plot_joint_mediation <- function(draws_med,
                                 mediator1_suffix, mediator2_suffix,
                                 mediator1_nm, mediator2_nm,
                                 mediator_cols,
                                 diff_lab = "difference\n(BA − CR)",
                                 brew_col = "Archambault",
                                 col_nums = c(1, 4),
                                 box_fill_alpha = 1,
                                 flip_order = FALSE,
                                 title = NULL,
                                 title_col = "grey15",
                                 med_col_text = "white",
                                 direct_lab = NULL,
                                 fnt_sz = 1) {
  if (!all(c(mediator1_suffix, mediator2_suffix) %in% names(mediator_cols))) {
    stop(
      "plot_joint_mediation(): `mediator_cols` must be named with `mediator1_suffix`/",
      "`mediator2_suffix` (\"", mediator1_suffix, "\"/\"", mediator2_suffix, "\") -- got names: ",
      paste(names(mediator_cols), collapse = ", ")
    )
  }
  pal <- MetBrewer::met.brewer(brew_col)[col_nums]
  grp_cols <- c("behavioural activation" = pal[1], "cognitive restructuring" = pal[2])

  # --- path diagram skeleton: mediator 1 arcs above the c' path, mediator 2 mirrors below ---
  hw <- 0.95
  hh <- 0.52
  mediator1_lab <- sub(" ", "\n", mediator1_nm)
  mediator2_lab <- sub(" ", "\n", mediator2_nm)
  # M1/M2 boxes are filled with their own colour at box_fill_alpha so they double as the
  # key for plot_joint_indirect()'s combined panel; X/Y boxes stay neutral grey. Box-label
  # text colour is picked per-box from how the fill will actually look blended with the
  # (white) background, so it stays legible at any hue/alpha combination.
  med_cols <- unname(mediator_cols[c(mediator1_suffix, mediator2_suffix)])
  boxes <- data.frame(
    x    = c(0.0, 4.0, 8.0, 4.0),
    y    = c(0.0, 1.8, 0.0, -1.8),
    fill = c("grey97", med_cols[1], "grey97", med_cols[2]),
    alpha = c(1, box_fill_alpha, 1, box_fill_alpha),
    text_col = c("grey15", med_col_text, "grey15", med_col_text),
    lab  = c(
      "iCBT modules", paste0("Δ ", mediator1_lab),
      "Δ symptoms<br>(PHQ-9<sub>std</sub>)", paste0("Δ ", mediator2_lab)
    )
  )
  # The outcome box carries a "<sub>std</sub>" scaling suffix, so the two PLAIN boxes (X and Y) are
  # drawn as richtext. The MEDIATOR boxes are not: their labels are the caller's display names with
  # sub(" ", "\n", .) applied, i.e. real newlines, which an HTML renderer would collapse to a space.
  # Splitting the layer keeps the caller's interface unchanged.
  plain_boxes <- boxes[c(1, 3), ]
  med_boxes <- boxes[c(2, 4), ]
  # arrow segments: a1/b1 (upper path, mediator 1) and a2/b2 (lower path, mediator 2)
  # solid; c' (X->Y) faint dashed, nudged below the X/Y boxes' centreline (rather than
  # straight through it) so its own inset has a clear, unobstructed slot above it
  seg_ab <- data.frame(
    x    = c(0.95, 4.95, 0.95,  4.95),
    y    = c(0.52,  1.28, -0.52, -1.28),
    xend = c(3.05, 7.03, 3.05,  7.03),
    yend = c(1.28,  0.52, -1.28, -0.52)
  )
  seg_c <- data.frame(x = 0.95, y = -0.35, xend = 7.05, yend = -0.35)
  # path letters near each arrow (anchor the matching inset). a/b labels mirror the
  # paper's notation (a_g^(k)/b_g^(k), group g and parameter k) via a parsed plotmath
  # expression -- deliberately generic (no actual group/parameter substituted in) since
  # the nearby box labels already say which mediator/parameter this path belongs to.
  # Rendered as its own geom_text layer (parse = TRUE) since plotmath expressions can't
  # mix with the plain-text "c′" label in one layer, and ggplot2's `fontface` aesthetic
  # is ignored under parse = TRUE -- italics have to be requested inside the expression
  # itself (`italic(...)`) to match the plain c′ label's fontface = "italic" styling.
  plett_ab <- data.frame(
    x = c(2.1, 5.85, 2.1, 5.85), y = c(0.64, 0.64, -0.64, -0.64),
    lab = c("italic(a[g]^{(k)})", "italic(b[g]^{(k)})", "italic(a[g]^{(k)})", "italic(b[g]^{(k)})")
  )
  plett_c <- data.frame(x = 4, y = -0.55, lab = "c′")

  diagram <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = boxes,
      ggplot2::aes(xmin = x - hw, xmax = x + hw, ymin = y - hh, ymax = y + hh, fill = fill, alpha = alpha),
      colour = "grey45", linewidth = 0.5
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::geom_segment(
      data = seg_c, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linetype = "dashed", colour = "grey60", linewidth = 0.5,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.018, "npc"), type = "closed")
    ) +
    ggplot2::geom_segment(
      data = seg_ab, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = "grey25", linewidth = 0.7,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.022, "npc"), type = "closed")
    ) +
    ggtext::geom_richtext(
      data = plain_boxes, ggplot2::aes(x = x, y = y, label = lab, colour = text_col),
      family = "Open Sans", size = 4.6 * fnt_sz, lineheight = 1.1,
      fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt")
    ) +
    ggplot2::geom_text(
      data = med_boxes, ggplot2::aes(x = x, y = y, label = lab, colour = text_col),
      family = "Open Sans", size = 4.6 * fnt_sz, lineheight = 0.9
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(
      data = plett_ab, ggplot2::aes(x = x, y = y, label = lab),
      family = "Open Sans", colour = "grey35", size = 4 * fnt_sz, parse = TRUE
    ) +
    ggplot2::geom_text(
      data = plett_c, ggplot2::aes(x = x, y = y, label = lab),
      family = "Open Sans", fontface = "italic", colour = "grey35", size = 4.8 * fnt_sz
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-1.3, 9.4), ylim = c(-2.55, 2.55), clip = "off", expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

  # optional bold title in the diagram's own whitespace above the M1 box (e.g. valence
  # labels for causal attribution); drawn as a fixed-colour annotation so it's unaffected
  # by the box-label scale_colour_identity() above
  if (!is.null(title)) {
    diagram <- diagram + ggplot2::annotate(
      "text", x = 0, y = 2.7, label = title, colour = title_col,
      family = "Open Sans", fontface = "bold", size = 5.4 * fnt_sz, hjust = 0.5
    )
  }

  # --- inset interval plots for each path (CR on top / BA below by default; see flip_order) ---
  a1_path <- mini_interval_plot(
    draws_med[[paste0("a_", mediator1_suffix, "_ba")]], draws_med[[paste0("a_", mediator1_suffix, "_cr")]],
    paste0("module → ", mediator1_nm), grp_cols, fnt_sz, flip_order = flip_order
  )
  b1_path <- mini_interval_plot(
    draws_med[[paste0("b_", mediator1_suffix, "_ba")]], draws_med[[paste0("b_", mediator1_suffix, "_cr")]],
    paste0(mediator1_nm, " → symptoms"), grp_cols, fnt_sz, flip_order = flip_order
  )
  a2_path <- mini_interval_plot(
    draws_med[[paste0("a_", mediator2_suffix, "_ba")]], draws_med[[paste0("a_", mediator2_suffix, "_cr")]],
    paste0("module → ", mediator2_nm), grp_cols, fnt_sz, flip_order = flip_order
  )
  b2_path <- mini_interval_plot(
    draws_med[[paste0("b_", mediator2_suffix, "_ba")]], draws_med[[paste0("b_", mediator2_suffix, "_cr")]],
    paste0(mediator2_nm, " → symptoms"), grp_cols, fnt_sz, flip_order = flip_order
  )
  # --- companion: both mediators' indirect effects on one panel (dodged, mediator 1 above mediator 2) ---
  indirect <- plot_joint_indirect(
    draws_med, mediator1_suffix, mediator2_suffix, mediator_cols,
    mediator1_nm = mediator1_nm, mediator2_nm = mediator2_nm,
    diff_lab = diff_lab, flip_order = flip_order, fnt_sz = fnt_sz
  )

  # residual (c') path -- what's left of the module -> symptoms association once both
  # mediators are adjusted for; coloured by arm (grp_cols) like the a/b insets, not by
  # mediator, since it doesn't belong to either one.
  #
  # TWO SHAPES, depending on the symptom submodel the draws came from:
  #   * the compound-symmetry models emit a per-arm pair (direct_ba = beta_phq_time,
  #     direct_cr = that plus the arm x time interaction), drawn as two rows;
  #   * the ANCOVA models (...-joint-decomp-ancova.stan) emit a single `direct_arm`, because their
  #     outcome is the t2 LEVEL and there is no within-arm change term to split per arm. that is
  #     drawn as ONE row in the neutral difference colour.
  # `direct_lab` names that single row; left NULL it follows flip_order, which is how both call
  # sites already encode arm direction (flip_order = TRUE for reward-effort, whose contrast is
  # BA - CR; FALSE for causal attribution, whose contrast is CR - BA).
  has_arm_pair <- all(c("direct_ba", "direct_cr") %in% names(draws_med))
  if (!has_arm_pair && !"direct_arm" %in% names(draws_med)) {
    stop(
      "plot_joint_mediation(): `draws_med` has neither the per-arm pair (direct_ba/direct_cr) nor ",
      "the ANCOVA single direct effect (direct_arm) -- was `direct_*` included in the variables ",
      "pulled from the fit?"
    )
  }
  if (is.null(direct_lab)) direct_lab <- if (flip_order) "BA \u2212 CR" else "CR \u2212 BA"
  direct_path <- if (has_arm_pair) {
    mini_interval_plot(
      draws_med$direct_ba, draws_med$direct_cr, "direct effect", grp_cols, fnt_sz,
      flip_order = flip_order
    )
  } else {
    mini_interval_plot(
      draws_med$direct_arm, NULL, "direct effect", grp_cols, fnt_sz,
      flip_order = flip_order, single_lab = direct_lab
    )
  }

  list(
    diagram = diagram,
    a1_path = a1_path, b1_path = b1_path,
    a2_path = a2_path, b2_path = b2_path,
    direct_path = direct_path,
    indirect = indirect
  )
}

# N-mediator generalisation of plot_joint_indirect(): one dodged pointrange per
# mediator within each group row (mediators ordered top->bottom by `suffixes`,
# matching the diagram). Unlike the 2-mediator version, the mediator colours are
# shown as a small legend, because the grouped-box diagram it accompanies can't
# tint a per-mediator box to serve as the key.
#   suffixes:      mediator column suffixes, top->bottom (e.g. c("rew","eff","alpha")).
#   mediator_nms:  display names, same length/order as `suffixes` (legend labels).
#   mediator_cols: named vector, names == `suffixes`.
plot_joint_indirect_n <- function(draws_med,
                                  suffixes,
                                  mediator_nms,
                                  mediator_cols,
                                  diff_lab = "difference\n(BA − CR)",
                                  flip_order = FALSE,
                                  legend = TRUE,
                                  fnt_sz = 1) {
  if (!all(suffixes %in% names(mediator_cols))) {
    stop(
      "plot_joint_indirect_n(): `mediator_cols` must be named with every entry of ",
      "`suffixes` (", paste(suffixes, collapse = ", "), ") -- got names: ",
      paste(names(mediator_cols), collapse = ", ")
    )
  }
  n_med <- length(suffixes)
  row_vals <- function(suffix, level, value) {
    data.frame(suffix = suffix, level = level, value = value)
  }
  ind_df <- dplyr::bind_rows(lapply(suffixes, function(sfx) {
    dplyr::bind_rows(
      row_vals(sfx, "behavioural activation",  draws_med[[paste0("indirect_", sfx, "_ba")]]),
      row_vals(sfx, "cognitive restructuring", draws_med[[paste0("indirect_", sfx, "_cr")]]),
      row_vals(sfx, "difference",              draws_med[[paste0("index_mod_med_", sfx)]])
    )
  }))

  # unflipped = CR on top; flip_order = TRUE = BA on top (as elsewhere)
  other <- if (flip_order) {
    c("cognitive restructuring", "behavioural activation")
  } else {
    c("behavioural activation", "cognitive restructuring")
  }
  ind_df$level <- factor(ind_df$level, levels = c("difference", other))
  # dodge mediators within each row: first suffix highest, last lowest
  step <- min(0.22, 0.66 / n_med)
  offs <- ((n_med + 1) / 2 - seq_len(n_med)) * step
  names(offs) <- suffixes
  ind_df$ybase <- as.numeric(ind_df$level) + offs[ind_df$suffix]
  ind_df$mediator <- factor(ind_df$suffix, levels = suffixes, labels = mediator_nms)
  leg_cols <- stats::setNames(unname(mediator_cols[suffixes]), mediator_nms)

  ind_df |>
    ggplot2::ggplot(
      ggplot2::aes(x = value, y = ybase, colour = mediator, group = interaction(mediator, level))
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1.5, linetype = "dotted", colour = "grey70", linewidth = 2) +
    ggdist::stat_pointinterval(
      .width = c(0.95, 0.99), point_interval = "mean_hdci", interval_size_range = c(1, 3), fatten_point = 2
    ) +
    ggplot2::scale_colour_manual(
      values = leg_cols, name = NULL, guide = if (legend) "legend" else "none"
    ) +
    ggplot2::scale_y_continuous(
      breaks = 1:3, labels = c(diff_lab, gsub(" ", "\n", other)), limits = c(0.5, 3.5)
    ) +
    ggplot2::labs(x = "indirect path (a × b): Δ symptoms via Δ parameters", y = NULL) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = 16 * fnt_sz) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12 * fnt_sz),
      axis.text.y = ggplot2::element_text(size = 10 * fnt_sz),
      legend.position = if (legend) "bottom" else "none",
      legend.text = ggplot2::element_text(size = 11 * fnt_sz),
      legend.justification = "center"
    )
}

# Grouped three-mediator path diagram for the extended joint-decomp model
# (…-symptoms-joint-decomp-acceptbias.stan). Keeps the clean two-box look: the
# TOP box holds the two value-sensitivity mediators (reward + effort) with a
# *pair* of a-arrows in and b-arrows out, each coloured by its mediator; the
# BOTTOM box is the acceptance-bias mediator with a single a/b arc. The direct
# (c') path runs straight along the centreline between the boxes. The mediator
# colours live on the arrows and on the companion indirect panel's legend (a
# single grouped box can't double as a per-mediator key like the 2-mediator
# diagram does). `mediator_cols` must be named with all three suffixes.
# Returns list(diagram, indirect) -- compact (no per-path insets).
plot_joint_mediation_grouped <- function(draws_med,
                                         top_suffixes = c("rew", "eff"),
                                         top_nms = c("reward sensitivity", "effort sensitivity"),
                                         top_box_lab = "Δ reward /\neffort sensitivity",
                                         bottom_suffix = "alpha",
                                         bottom_nm = "acceptance bias",
                                         bottom_box_lab = "Δ acceptance\nbias",
                                         mediator_cols,
                                         diff_lab = "difference\n(BA − CR)",
                                         flip_order = FALSE,
                                         fnt_sz = 1) {
  all_suffixes <- c(top_suffixes, bottom_suffix)
  if (!all(all_suffixes %in% names(mediator_cols))) {
    stop(
      "plot_joint_mediation_grouped(): `mediator_cols` must be named with every mediator ",
      "suffix (", paste(all_suffixes, collapse = ", "), ") -- got names: ",
      paste(names(mediator_cols), collapse = ", ")
    )
  }
  c_rew   <- mediator_cols[[top_suffixes[1]]]
  c_eff   <- mediator_cols[[top_suffixes[2]]]
  c_alpha <- mediator_cols[[bottom_suffix]]

  hw <- 0.95
  hh <- 0.52
  boxes <- data.frame(
    x        = c(0.0, 4.0, 8.0, 4.0),
    y        = c(0.0, 1.8, 0.0, -1.8),
    fill     = c("grey97", "grey94", "grey97", "grey94"),
    text_col = c("grey15", "grey15", "grey15", "grey15"),
    lab      = c("iCBT modules", top_box_lab, "Δ symptoms<br>(PHQ-9<sub>std</sub>)", bottom_box_lab)
  )
  # as in plot_joint_mediation(): richtext for the plain X/Y boxes (the outcome carries a
  # "<sub>std</sub>" suffix), geom_text for the caller-supplied mediator labels, whose "\n" line
  # breaks HTML would collapse
  plain_boxes <- boxes[c(1, 3), ]
  med_boxes <- boxes[c(2, 4), ]
  # top box: paired a-arrows (X->M) and b-arrows (M->Y), one per top mediator,
  # coloured by mediator; bottom box: single a/b arrows for acceptance bias
  seg_ab <- data.frame(
    x    = c(0.98, 0.98, 4.98, 4.98,  0.98,  4.98),
    y    = c(0.45, 0.15, 1.95, 1.60, -0.40, -1.70),
    xend = c(3.02, 3.02, 7.02, 7.02,  3.02,  7.02),
    yend = c(1.95, 1.60, 0.45, 0.15, -1.70, -0.40),
    col  = c(c_rew, c_eff, c_rew, c_eff, c_alpha, c_alpha)
  )
  seg_c <- data.frame(x = 0.95, y = 0, xend = 7.05, yend = 0)
  plett <- data.frame(
    x   = c(1.95, 6.05, 1.95, 6.05, 4.0),
    y   = c(1.35, 1.35, -1.18, -1.18, 0.18),
    lab = c("a", "b", "a", "b", "c′")
  )

  diagram <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = boxes,
      ggplot2::aes(xmin = x - hw, xmax = x + hw, ymin = y - hh, ymax = y + hh, fill = fill),
      colour = "grey45", linewidth = 0.5
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_segment(
      data = seg_c, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linetype = "dashed", colour = "grey60", linewidth = 0.5,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.018, "npc"), type = "closed")
    ) +
    ggplot2::geom_segment(
      data = seg_ab, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = col),
      linewidth = 0.8,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"), type = "closed")
    ) +
    ggtext::geom_richtext(
      data = plain_boxes, ggplot2::aes(x = x, y = y, label = lab, colour = text_col),
      family = "Open Sans", size = 4.6 * fnt_sz, lineheight = 0.9,
      fill = NA, label.color = NA, label.padding = grid::unit(rep(0, 4), "pt")
    ) +
    ggplot2::geom_text(
      data = med_boxes, ggplot2::aes(x = x, y = y, label = lab, colour = text_col),
      family = "Open Sans", size = 4.6 * fnt_sz, lineheight = 0.9
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_text(
      data = plett, ggplot2::aes(x = x, y = y, label = lab),
      family = "Open Sans", fontface = "italic", colour = "grey35", size = 4.8 * fnt_sz
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-1.3, 9.4), ylim = c(-2.55, 2.55), clip = "off", expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

  indirect <- plot_joint_indirect_n(
    draws_med,
    suffixes = all_suffixes,
    mediator_nms = c(top_nms, bottom_nm),
    mediator_cols = mediator_cols,
    diff_lab = diff_lab, flip_order = flip_order, fnt_sz = fnt_sz
  )

  list(diagram = diagram, indirect = indirect)
}

# Session-1 completer vs. non-completer parameter comparison.
# Three slab rows: completers mean, non-completers mean, and their difference.
# `draws` must contain columns "mu_{param}[1]" and "{param}_ncpl".
# Returns a single ggplot (call plot_change_slabs directly).
plot_ncpl_diff <- function(draws,
                           param,
                           param_nm,
                           mu_col          = NULL,
                           ncpl_col        = NULL,
                           diff_lab        = "difference in\nnon-completers",
                           xlim            = NULL,
                           flip_order      = FALSE,
                           suppress_labels = FALSE,
                           brew_col        = "Archambault",
                           col_nums        = c(2, 5),
                           diff_col        = "#bb9bb4",
                           fnt_sz          = 1) {
  if (is.null(mu_col))   mu_col   <- paste0("mu_", param, "[1]")
  if (is.null(ncpl_col)) ncpl_col <- paste0(param, "_ncpl")

  if (flip_order) {
    grp_shapes <- c(21, 25, 23)
  } else {
    grp_shapes <- c(21, 23, 25)
  }

  pal      <- MetBrewer::met.brewer(brew_col)[col_nums]
  grp_cols <- c(
    "completers"     = pal[1],
    "non-completers" = pal[2],
    "difference"     = diff_col
  )
  title <- if (flip_order) {
    paste0("session 1: non-completers vs. completers\n(", param_nm, ")")
  } else {
    paste0("session 1: completers vs. non-completers\n(", param_nm, ")")
  }

  comp_vec  <- draws[[mu_col]]
  ncpl_diff <- draws[[ncpl_col]]

  chg <- dplyr::bind_rows(
    data.frame(level = "completers",     value = comp_vec),
    data.frame(level = "non-completers", value = comp_vec + ncpl_diff),
    data.frame(level = "difference",     value = ncpl_diff)
  )
  chg$level <- factor(
    chg$level, levels = c("difference", "non-completers", "completers")
  )

  plot_change_slabs(
    chg, xlim = xlim,
    x_lab = param_nm,
    grp_shapes = grp_shapes,
    grp_cols = grp_cols, fnt_sz = fnt_sz, flip_order = flip_order,
    suppress_labels = suppress_labels,
    subtitle = title
  ) +
    ggplot2::scale_y_discrete(labels = c(
      "completers"     = "completers",
      "non-completers" = "non-completers",
      "difference"     = diff_lab
    ))
}

# Session-1 baseline symptom/questionnaire scores: completers vs. non-completers,
# one dodged bar pair per measure. Designed to replace a "no baseline differences"
# table -- `value_col` should already be z-scored *within each measure* (so
# differently-scaled questionnaires, e.g. PHQ-9 /27 vs. AMI-BA /4, sit on one common
# axis); z-scoring is left to the caller (e.g. `dplyr::group_by(measure) |>
# dplyr::mutate(value_z = as.numeric(scale(value)))`) rather than done here, since the
# "pool across whom" question (all enrolled vs. completers only) is a data decision,
# not a plotting one. Same "completers"/"non-completers" naming + default Hokusai3
# palette (col_nums = c(2, 5)) as plot_ncpl_diff(), so this reads as the same session-1
# non-completer comparison already established elsewhere in fig_5.
#   df: long-format, one row per participant x measure.
#   measure_col: column with the measure label (x-axis categories, e.g. "PHQ-9").
#   value_col: the (already z-scored) score column.
#   completer_col: logical/0-1 column, TRUE/1 = non-completer.
plot_ncpl_symptom_bars <- function(df,
                                   measure_col = "measure",
                                   value_col = "value_z",
                                   completer_col = "non_completer",
                                   error_type = c("sd", "se"),
                                   bar_width = 0.65,
                                   err_width = 0.18,
                                   dodge_width = 0.75,
                                   alpha_bar = 0.8,
                                   fnt_sz = 1,
                                   legend_pos = "bottom",
                                   brew_col = "Archambault",
                                   col_nums = c(2, 5)) {
  error_type <- match.arg(error_type)
  pal <- MetBrewer::met.brewer(brew_col)[col_nums]
  grp_cols <- c("completers" = pal[1], "non-completers" = pal[2])

  plt_df <- df |>
    dplyr::transmute(
      measure = .data[[measure_col]],
      value = .data[[value_col]],
      group_plot = factor(
        ifelse(as.logical(.data[[completer_col]]), "non-completers", "completers"),
        levels = c("completers", "non-completers")
      )
    ) |>
    dplyr::filter(!is.na(value), !is.na(group_plot))

  sum_df <- plt_df |>
    dplyr::group_by(measure, group_plot) |>
    dplyr::summarise(
      mean = mean(value), sd = stats::sd(value), se = sd / sqrt(dplyr::n()), .groups = "drop"
    ) |>
    dplyr::mutate(err = if (error_type == "sd") sd else se)

  pd <- ggplot2::position_dodge(width = dodge_width)
  ggplot2::ggplot(sum_df, ggplot2::aes(x = measure, y = mean, fill = group_plot, colour = group_plot)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "#cacaca", alpha = 0.8) +
    ggplot2::geom_col(position = pd, width = bar_width, alpha = alpha_bar, linewidth = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - err, ymax = mean + err),
      position = pd, width = err_width, linewidth = 0.5, alpha = 0.5
    ) +
    ggplot2::scale_colour_manual(name = "", values = grp_cols) +
    ggplot2::scale_fill_manual(name = "", values = grp_cols) +
    ggplot2::labs(
      x = NULL,
      y = paste0("z-score (± ", ifelse(error_type == "sd", "s.d.", "s.e."), ")")
    ) +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = 14 * fnt_sz) +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.justification = "center"
    )
}

# Group-level population mean + pre→post change for a causal attribution parameter.
# Mirrors plot_reff_change() but uses the CA model's naming convention:
#   delta_{param}      = behavioural activation (BA) group change
#   theta_int_{param}  = additional change in cognitive restructuring (CR) vs BA
# `param` = e.g. "internal_neg". Top and bottom panels are both on the raw
# (untransformed) logit scale -- same scale as the reported delta_{param}/
# theta_int_{param} statistics -- since top is now built from the same
# mu_{mu_prefix}[t]/theta_int_{param} group-level draws as the bottom panel,
# rather than averaging individual theta_{param}[id,t] draws. `mu_prefix`
# (e.g. "mu_internal_theta_neg") is needed because the Stan mu_* names put
# "internal"/"global" *before* "theta", unlike `param`'s "internal_neg"/
# "global_neg" ordering, so it can't be derived by pasting onto `param` --
# same reason plot_ncpl_diff() takes an explicit `mu_col` override.
# Returns list(levels_panel, change_panel).
plot_ca_change <- function(draws,
                           param,
                           param_nm,
                           flip_order = FALSE, # CR on top default
                           mu_prefix = NULL,
                           diff_lab = "difference\n(CR − BA)",
                           level_xlim = NULL,
                           change_xlim = NULL,
                           suppress_labels = FALSE,
                           brew_col = "Archambault",
                           col_nums = c(1, 4),
                           diff_col = "#bb9bb4",
                           fnt_sz = 1) {
  if (is.null(mu_prefix)) mu_prefix <- paste0("mu_", param)
  pal <- MetBrewer::met.brewer(brew_col)[col_nums]
  grp_cols <- c(
    "behavioural activation" = pal[1], "cognitive restructuring" = pal[2], "difference" = diff_col
  )
  if (flip_order) {
    grp_shapes = c(16, 17, 15)
  } else {
    grp_shapes = c(16, 15, 17)
  }
  nudge <- 0.16

  # --- top: group-level means (mu_prefix[t], + theta_int_{param} for CR @ t2) with 95%/99% HDI ---
  # For CA: condition 1 = CR (cognitive restructuring), 0 = BA (behavioural activation, reference)
  mu1 <- draws[[paste0(mu_prefix, "[1]")]]  # pre, shared by both arms
  mu2 <- draws[[paste0(mu_prefix, "[2]")]]  # post, BA (reference)
  int_col <- paste0("theta_int_", param)
  int <- draws[[int_col]]                   # additional CR effect at post

  pop_rows <- list()
  for (g in c(1, 0)) {
    for (t in 1:2) {
      v <- if (t == 1) mu1 else if (g == 1) mu2 + int else mu2
      h95 <- bayestestR::hdi(v, ci = 0.95)
      h99 <- bayestestR::hdi(v, ci = 0.99)
      pop_rows[[length(pop_rows) + 1]] <- data.frame(
        group = ifelse(g == 1, "cognitive restructuring", "behavioural activation"),
        tp    = ifelse(t == 1, "pre", "post"),
        m     = median(v),
        lo95  = h95$CI_low, hi95 = h95$CI_high,
        lo99  = h99$CI_low, hi99 = h99$CI_high
      )
    }
  }
  pop <- do.call(rbind, pop_rows)
  pop$group <- factor(pop$group, levels = c("behavioural activation", "cognitive restructuring"))
  pop$ybase <- as.numeric(pop$group) + ifelse(pop$tp == "pre", nudge, -nudge)

  # --- bottom: group-level change (slab + 95%/99% HDI) ---
  # For CA: BA = reference (condition 0), so delta_{param} = BA change;
  # theta_int_{param} = additional change in CR (cognitive restructuring) vs BA.
  chg <- dplyr::bind_rows(
    data.frame(level = "behavioural activation",  value = draws[[paste0("delta_", param)]]),
    data.frame(
      level = "cognitive restructuring",
      value = draws[[paste0("delta_", param)]] + draws[[int_col]]
    ),
    data.frame(level = "difference", value = draws[[int_col]])
  )
  chg$level <- factor(
    chg$level, levels = c("difference", "behavioural activation", "cognitive restructuring")
  )

  list(
    levels_panel = plot_group_levels(
      pop, xlim = level_xlim,
      x_lab    = param_nm,
      grp_cols = grp_cols, fnt_sz = fnt_sz,
      flip_order = flip_order,
      suppress_labels = suppress_labels,
      subtitle = "group-level estimates\n○ pre- ● post-intervention"
    ),
    change_panel = plot_change_slabs(
      chg, xlim = change_xlim, diff_comp = diff_lab,
      x_lab    = paste0("change in ", param_nm),
      grp_shapes = grp_shapes, flip_order = flip_order,
      grp_cols = grp_cols, fnt_sz = fnt_sz,
      suppress_labels = suppress_labels
    )
  )
}

# Posterior correlations between model parameters and symptom scores, computed
# for three timepoints (`t1`, `t2` = per-timepoint levels; `change` = pre-post
# deltas) and four groupings of the `group` dimension:
#   "BA" / "CR"  -- within each intervention arm
#   "both"       -- pooled across arms
#   "difference" -- between-arm difference in correlation (BA - CR), per draw
# Returns a list with:
#   $cor_df     -- per-draw correlations (long: group, par, questionnaire,
#                  timepoint, .draw, correlation)
#   $summary_df -- posterior mean + 95% HDI per (group, par, questionnaire,
#                  timepoint); the canonical input to plot_correlation_heatmap()
#   $panels     -- $summary_df split by group ($BA/$CR/$both/$difference), each a
#                  ready-to-plot slice (pass with int_group = NULL, or ignore and
#                  filter $summary_df yourself via int_group)
get_symptom_corrs <- function(data_long,
                              model_type,
                              draws,
                              stan_ls,
                              qqnrs,
                              delta_par_pattern,
                              par_labels,
                              indiv_par_map = NULL,
                              qns_to_invert = c()) {
  # Prepare questionnaire data: pre/post levels ("_t1"/"_t2") as well as
  # pre-post change ("_delta"), all inverted consistently for measures where
  # a higher score reflects lower symptom burden.
  invert_bases <- unique(gsub("_(delta|t1|t2)$", "", qns_to_invert))
  invert_pattern <- paste0("^(", paste(invert_bases, collapse = "|"), ")_(t1|t2|delta)$")

  qns_wide <- qqnrs |>
    dplyr::select(subID, group, tidyselect::matches("_total_(t1|t2|delta)$")) |>
    dplyr::rename_with(~ gsub("_total_", "_", .)) |>
    dplyr::mutate(
      dplyr::across(tidyselect::matches(invert_pattern), ~ . * -1)
    )

  # Prepare individual parameter draws
  ids_dt <- prep_data(data_long, type = model_type, ret_ids = TRUE) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::select(id, subID) |>
    data.table::as.data.table()

  # long parameter-draws table (id/par/par_value/.draw) x long questionnaire
  # table (subID/group/questionnaire/score) -> posterior correlations, within
  # group x par x questionnaire x draw.
  #
  # Correlation is computed from closed-form sums (n, sum(x), sum(y), sum(xy),
  # sum(x^2), sum(y^2)) via data.table group-by aggregation instead of
  # tidyr::nest() + purrr::map_dbl(cor); the latter materialises one tibble
  # per (group, par, questionnaire, .draw) combination -- millions of them at
  # full draw counts -- and calls cor() on each, which dominates runtime.
  # Those sums are additive across disjoint subsets, so the pooled "both"
  # group is derived by summing the already-aggregated BA/CR rows rather than
  # re-scanning the full joined table a second time.
  corr_from_sums <- function(agg) {
    agg[, den := sqrt(pmax(n * sxx - sx^2, 0) * pmax(n * syy - sy^2, 0))]
    agg[, correlation := data.table::fifelse(
      n < 2 | den == 0, NA_real_, (n * sxy - sx * sy) / den
    )]
    agg[, .(group, par, questionnaire, .draw, correlation)]
  }

  correlate_draws <- function(par_draws, qn_long, timepoint_label) {
    pd <- data.table::as.data.table(par_draws)
    pd[, id := as.character(id)]
    pd_ids <- merge(pd, ids_dt, by = "id", all.x = TRUE, allow.cartesian = TRUE)
    joined <- merge(
      data.table::as.data.table(qn_long), pd_ids,
      by = "subID", all.x = TRUE, allow.cartesian = TRUE
    )
    joined[, c("id", "subID") := NULL]
    joined <- joined[!is.na(par_value) & !is.na(score)]

    agg_grp <- joined[
      , .(
        n   = .N,
        sx  = sum(par_value),
        sy  = sum(score),
        sxy = sum(par_value * score),
        sxx = sum(par_value * par_value),
        syy = sum(score * score)
      ),
      by = .(group, par, questionnaire, .draw)
    ]
    # per-group correlations, plus a "both" pool across groups -- useful e.g.
    # at t1 where BA/CR haven't yet diverged and splitting is uninformative
    agg_both <- agg_grp[
      , .(n = sum(n), sx = sum(sx), sy = sum(sy), sxy = sum(sxy), sxx = sum(sxx), syy = sum(syy)),
      by = .(par, questionnaire, .draw)
    ][, group := "both"]

    out <- corr_from_sums(rbind(agg_grp, agg_both, use.names = TRUE, fill = TRUE))
    out[, timepoint := timepoint_label]
    tibble::as_tibble(out)
  }

  qn_long_for <- function(suffix) {
    qns_wide |>
      dplyr::select(subID, group, tidyselect::ends_with(suffix)) |>
      tidyr::drop_na(tidyselect::ends_with(suffix)) |>
      tidyr::pivot_longer(
        cols = tidyselect::ends_with(suffix),
        names_to = "questionnaire",
        values_to = "score"
      )
  }

  ## ---- change: delta parameter vs. delta symptom score --------------------
  delta_draws <- draws |>
    dplyr::select(tidyselect::starts_with("delta_"), .draw)

  if (length(par_labels) > 1) {
    delta_draws <- delta_draws |>
      tidyr::pivot_longer(
        cols = -".draw",
        names_to = c("par", "id"),
        names_pattern = delta_par_pattern,
        values_to = "par_value"
      )
  } else {
    delta_draws <- delta_draws |>
      tidyr::pivot_longer(
        cols = -".draw",
        names_to = "id",
        names_pattern = delta_par_pattern,
        values_to = "par_value"
      ) |>
      dplyr::mutate(par = names(par_labels))
  }
  delta_draws <- delta_draws |>
    dplyr::mutate(id = as.integer(id))

  change_cor <- correlate_draws(delta_draws, qn_long_for("_delta"), "change")

  ## ---- levels: parameter vs. symptom score, within each timepoint --------
  level_draws <- draws |>
    dplyr::select(-tidyselect::starts_with("delta_"), -tidyselect::any_of(c(".chain", ".iteration"))) |>
    tidyr::pivot_longer(
      cols = -".draw",
      names_to = c("raw_par", "id", "tp"),
      names_pattern = "^(.*)\\[(\\d+),\\s*(\\d+)\\]$",
      values_to = "par_value"
    ) |>
    dplyr::mutate(
      id = as.integer(id),
      par = if (is.null(indiv_par_map)) raw_par else dplyr::recode(raw_par, !!!indiv_par_map),
      timepoint = dplyr::case_when(tp == "1" ~ "t1", tp == "2" ~ "t2", TRUE ~ NA_character_)
    ) |>
    dplyr::filter(par %in% names(par_labels)) |>
    dplyr::select(-raw_par, -tp)

  level_cor <- lapply(c("t1", "t2"), function(tpt) {
    correlate_draws(
      level_draws |> dplyr::filter(timepoint == tpt) |> dplyr::select(-timepoint),
      qn_long_for(paste0("_", tpt)),
      tpt
    )
  }) |>
    dplyr::bind_rows()

  # per-draw correlations for the per-arm ("BA"/"CR") and pooled ("both") groups
  per_draw <- dplyr::bind_rows(change_cor, level_cor)

  # between-group difference in correlation (BA - CR), computed per draw so its
  # posterior/HDI properly reflects the joint uncertainty -- another value of the
  # `group` dimension, so it flows through the plot function like any other.
  diff_draw <- per_draw |>
    dplyr::filter(group %in% c("BA", "CR")) |>
    tidyr::pivot_wider(names_from = group, values_from = correlation) |>
    dplyr::mutate(correlation = BA - CR, group = "difference") |>
    dplyr::select(-BA, -CR)

  # combine everything, then apply the shared factor levels/labels
  posterior_cor_df <- dplyr::bind_rows(per_draw, diff_draw) |>
    dplyr::mutate(
      par = factor(par, levels = names(par_labels), labels = par_labels),
      group = factor(group, levels = c("BA", "CR", "both", "difference")),
      timepoint = factor(timepoint, levels = c("t1", "t2", "change"))
    )

  # Summarise posterior correlations
  cor_summary_df <- posterior_cor_df |>
    dplyr::summarise(
      mean_corr = mean(correlation),
      .lower = bayestestR::hdi(correlation, ci = 0.95)$CI_low,
      .upper = bayestestR::hdi(correlation, ci = 0.95)$CI_high,
      # p_direction drives the tile shading in plot_correlation_heatmap(): weak-evidence cells fade
      # toward white rather than reading as strongly coloured, with no hard cutoff anywhere
      p_direct = as.numeric(bayestestR::p_direction(correlation, method = "direct")),
      .by = c(group, par, questionnaire, timepoint)
    ) |>
    dplyr::rename(correlation = mean_corr) |>
    dplyr::ungroup()

  # `summary_df`/`cor_df` are the canonical long tables (every quantity is a
  # (group, timepoint) slice); `panels` is a convenience view splitting the
  # summary by the quantity of interest, each element ready to hand straight to
  # plot_correlation_heatmap() (with the matching int_group, or int_group = NULL).
  panels <- split(cor_summary_df, cor_summary_df$group)

  list(cor_df = posterior_cor_df, summary_df = cor_summary_df, panels = panels)
}

plot_correlation_heatmap <- function(cor_summary_df,
                                     int_group = NULL,
                                     par_labels,
                                     qn_labels,
                                     timepoint = c("change", "t1", "t2"),
                                     cor_lims = NULL,
                                     htmp_title = NULL,
                                     title_clr = "black",
                                     txt_sz = 4,
                                     fnt_sz = 1,
                                     label_type = c("full", "r_only", "sig_only", "none"),
                                     alpha_range = c(0.15, 1)) {
  label_type <- match.arg(label_type)
  timepoint <- match.arg(timepoint)
  # p_direct is added by get_symptom_corrs(); a frame cached before that change has no such column,
  # and without this guard the missing aes surfaces as an opaque ggplot error several frames down
  if (!"p_direct" %in% names(cor_summary_df)) {
    stop(
      "plot_correlation_heatmap(): `cor_summary_df` has no `p_direct` column, which the tile ",
      "shading maps. It predates the p_direction addition to get_symptom_corrs() -- re-run that ",
      "(or delete the cached summary) rather than plotting unshaded tiles.",
      call. = FALSE
    )
  }
  par_prefix <- if (timepoint == "change") "\u0394 " else ""

  # Reshape the correlation data for plotting. int_group filters the group
  # dimension ("BA"/"CR"/"both"/"difference"); pass int_group = NULL when the
  # supplied frame already holds a single group (e.g. one of `panels`).
  cor_plot_df <- cor_summary_df |>
    dplyr::filter(
      timepoint == !!timepoint,
      if (is.null(int_group)) TRUE else group %in% !!int_group
    ) |>
    dplyr::mutate(
      par = factor(par, labels = paste0(par_prefix, par_labels)),
      questionnaire = gsub("_(delta|t1|t2)$", "", questionnaire),
      questionnaire = factor(
        questionnaire,
        levels = rev(names(qn_labels)),
        labels = rev(qn_labels)
      )
    )

  if (is.null(cor_lims)) {
    cor_lims <- c(
      min(cor_plot_df$correlation) - 0.1,
      max(cor_plot_df$correlation) + 0.1
    )
  }

  cor_plot_df <- cor_plot_df |>
    dplyr::mutate(
      cell_label = if (label_type == "full") {
        sprintf("%.2f\n(%.2f, %.2f)", correlation, .lower, .upper)
      } else if (label_type == "r_only") {
        sprintf("%.2f", correlation)
      } else if (label_type == "sig_only") {
        ifelse(.lower > 0 | .upper < 0, sprintf("%.2f", correlation), "")
      } else {
        ""
      }
    )

  # alpha (not just fill) carries p_direction, so weak-evidence cells fade toward white rather than
  # reading as strongly coloured -- an uncertainty cue with no hard cutoff, which also makes these
  # panels legible at a distance (the fully-opaque version read as uniformly strong from afar)
  base_plt <- cor_plot_df |>
    ggplot2::ggplot(ggplot2::aes(x = par, y = questionnaire, fill = correlation, alpha = p_direct)) +
    ggplot2::geom_tile(color = "white", linewidth = 1.5)

  if (label_type != "none") {
    # show.legend = FALSE: this layer also inherits the alpha aes (so the printed r fades along with
    # its tile), but its own key glyph ("a") must not win the pd legend over geom_tile's rect
    base_plt <- base_plt + ggplot2::geom_text(
      ggplot2::aes(label = cell_label),
      color = "black", size = txt_sz, family = "Open Sans", lineheight = 0.8, show.legend = FALSE
    )
  }

  base_plt +
    ggplot2::scale_fill_gradientn(
      name = "<em>r</em>",
      limits = cor_lims,
      colors = c("#85D4E3", "white", "#E39485"),
      values = scales::rescale(c(cor_lims[1], 0, cor_lims[2]))
    ) +
    # pd is bounded in [0.5, 1] (bayestestR::p_direction never dips below 0.5, by definition), so the
    # breaks span that full range rather than implying a wider one
    ggplot2::scale_alpha_continuous(
      name = "<em>P</em><sub>direction</sub>", range = alpha_range, limits = c(0.5, 1),
      breaks = c(0.5, 0.75, 1),
      guide = ggplot2::guide_legend(override.aes = list(fill = "grey20"))
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(title = htmp_title, x = NULL, y = NULL) +
    cowplot::theme_minimal_grid(font_size = fnt_sz * 18, font_family = "Open Sans") +
    ggplot2::theme(
      axis.text.x = ggtext::element_markdown(size = fnt_sz * 14),
      axis.text.x.top = ggtext::element_markdown(size = fnt_sz * 14),
      axis.text.y = ggtext::element_markdown(size = fnt_sz * 16),
      legend.title = ggtext::element_markdown(size = fnt_sz * 12, face = "plain"),
      legend.text = ggplot2::element_text(size = fnt_sz * 11),
      legend.key.height = grid::unit(1.25, "lines"),
      plot.title = ggtext::element_markdown(
        family = "Open Sans SemiBold", face = "plain", hjust = 0.5,
        size = fnt_sz * 18, color = title_clr
      ),
      panel.grid.major = ggplot2::element_blank(),
      legend.position = "right"
    )
}

plot_reff_recovery <- function(
  sim_params,
  sim_data,
  sim_fit,
  stan_ls_sim,
  vars,
  pal = "Archambault",
  recode_effort = FALSE
) {
  sim_draws <- sim_fit$draws(format = "df", variables = vars)
  sim_recov <- extract_indiv_pars(sim_draws, stan_ls_sim)

  ids_sim <- prep_data(
    sim_data,
    type = "effort",
    ret_ids = TRUE,
    excl_eff_task = "none",
    excl_quest = "none",
    recode_effort = recode_effort
  )

  sim_true <- sim_params |>
    dplyr::mutate(timepoint = dplyr::if_else(timepoint == 0, "t1", "t2")) |>
    dplyr::select(subID, timepoint, dplyr::all_of(vars))

  sim_recov_wide <- sim_recov |>
    dplyr::filter(variable %in% vars, timepoint %in% c("t1", "t2")) |>
    dplyr::mutate(id = as.integer(id)) |>
    dplyr::left_join(ids_sim, by = c("id" = "id")) |>
    dplyr::select(subID, timepoint, variable, mean) |>
    tidyr::pivot_wider(names_from = variable, values_from = mean)

  sim_compare <- sim_true |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(vars),
      names_to = "param",
      values_to = "true"
    ) |>
    dplyr::left_join(
      sim_recov_wide |>
        tidyr::pivot_longer(
          cols = dplyr::all_of(vars),
          names_to = "param",
          values_to = "recovered"
        ),
      by = c("subID", "timepoint", "param"), suffix = c("_true", "_rec")
    ) |>
    dplyr::mutate(param = factor(param, levels = vars))

  sim_cor <- sim_compare |>
    dplyr::group_by(param, timepoint) |>
    dplyr::summarise(r = stats::cor(true, recovered, use = "pairwise.complete.obs"), .groups = "drop") |>
    dplyr::mutate(label = paste0("italic(r)~'='~", sprintf("%.2f", r)))

  scatter_plot <- sim_compare |>
    ggplot2::ggplot(ggplot2::aes(x = true, y = recovered, colour = param, fill = param)) +
    ggplot2::geom_point(alpha = 0.3, size = 2) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    ggplot2::geom_text(
      data = sim_cor,
      ggplot2::aes(x = -Inf, y = Inf, label = label),
      hjust = -0.05, vjust = 1.1, inherit.aes = FALSE,
      size = 6, family = "Open Sans", parse = TRUE
    ) +
    ggplot2::scale_colour_manual(
      name = "", values = MetBrewer::met.brewer(pal, 3)[seq_along(vars)]
    ) +
    ggplot2::scale_fill_manual(
      name = "", values = MetBrewer::met.brewer(pal, 3)[seq_along(vars)]
    ) +
    ggplot2::guides(colour = "none", fill = "none") +
    ggplot2::facet_grid(
      timepoint ~ param, scales = "free",
      labeller = ggplot2::labeller(
        timepoint = ggplot2::as_labeller(
          c(t1 = "italic('timepoint 1')", t2 = "italic('timepoint 2')"),
          ggplot2::label_parsed
        ),
        param = ggplot2::as_labeller(
          c(
            rewSens = "italic('reward sensitivity')",
            effSens = "italic('effort sensitivity')",
            alpha = "italic('acceptance bias')"
          ),
          ggplot2::label_parsed
        )
      )
    ) +
    ggplot2::labs(
      x = "true parameter value",
      y = "recovered posterior mean"
    ) +
    cowplot::theme_minimal_grid(font_family = "Open Sans", font_size = 18) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = 20)
    )

  label_map <- c(
    rewSens = "reward sensitivity",
    effSens = "effort sensitivity",
    alpha = "acceptance bias"
  )
  cor_labels <- ifelse(vars %in% names(label_map), label_map[vars], vars)

  build_cor_df <- function(df, label) {
    df |>
      dplyr::filter(timepoint %in% c("t1", "t2")) |>
      dplyr::group_by(timepoint) |>
      dplyr::group_modify(~ {
        mat <- stats::cor(.x[, vars, drop = FALSE], use = "pairwise.complete.obs")
        as.data.frame(as.table(mat)) |>
          dplyr::rename(param_x = Var1, param_y = Var2, correlation = Freq)
      }) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        source = label,
        param_x = factor(param_x, levels = vars, labels = cor_labels),
        param_y = factor(param_y, levels = vars, labels = cor_labels)
      )
  }

  cor_true_df <- build_cor_df(sim_true, "true")
  cor_recov_df <- build_cor_df(sim_recov_wide, "recovered")

  tri_df <- dplyr::bind_rows(cor_true_df, cor_recov_df) |>
    dplyr::mutate(
      x_idx = as.integer(param_x),
      y_idx = as.integer(param_y),
      keep = dplyr::case_when(
        source == "true" & y_idx > x_idx ~ TRUE,
        source == "recovered" & y_idx < x_idx ~ TRUE,
        TRUE ~ FALSE
      )
    ) |>
    dplyr::filter(keep & x_idx != y_idx)

  true_col <- "#1b9e77"
  rec_col <- "#d95f02"
  n_params <- length(cor_labels)
  diag_df <- tibble::tibble(
    timepoint = c("t1", "t2"),
    x = 0.5,
    y = 0.5,
    xend = n_params + 0.5,
    yend = n_params + 0.5
  )
  label_df <- tibble::tibble(
    timepoint = rep(c("t1", "t2"), each = 2),
    source = rep(c("true", "recovered"), times = 2),
    label = rep(c("true", "recovered"), times = 2),
    x = rep(c(cor_labels[1], cor_labels[length(cor_labels)]), times = 2),
    y = rep(c(cor_labels[length(cor_labels)], cor_labels[1]), times = 2),
    nudge_x = rep(c(-0.3, 0.15), times = 2),
    nudge_y = rep(c(0.7, -0.7), times = 2)
  )

  heatmap_plot <- tri_df |>
    ggplot2::ggplot(ggplot2::aes(x = param_x, y = param_y, fill = correlation)) +
    ggplot2::geom_tile(
      ggplot2::aes(color = source),
      linewidth = 1
    ) +
    ggplot2::geom_segment(
      data = diag_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      color = "white",
      linewidth = 3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", correlation)),
      size = 4.5, family = "Open Sans", color = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      low = "#2e294e",
      mid = "white",
      high = "#e71d36",
      midpoint = 0,
      name = "r"
    ) +
    ggplot2::scale_color_manual(values = c(true = true_col, recovered = rec_col), guide = "none") +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = label, color = source),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold",
      family = "Open Sans",
      nudge_x = label_df$nudge_x,
      nudge_y = label_df$nudge_y
    ) +
    ggplot2::facet_wrap(~ timepoint, labeller = ggplot2::as_labeller(
      c(t1 = "timepoint 1", t2 = "timepoint 2")
    )) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +
    cowplot::theme_minimal_grid(font_family = "Open Sans", font_size = 18) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )

  cor_plot_df <- dplyr::bind_rows(cor_true_df, cor_recov_df)

  list(scatter = scatter_plot, heatmap = heatmap_plot, cor_df = cor_plot_df)
}

# extract fixed effects
get_hdi_pd_brms <- function(mod, rel = "BA", var = NULL) {
  if (!is.null(var) && is.null(rel)) {
    var <- brms::fixef(mod, summary = FALSE)[, paste0(var, "TRUE")] |> exp()
    out <- c(
      mean = signif(mean(var), digits = 3),
      hdi_low = signif(bayestestR::hdi(var, ci = 0.95)[[2]], digits = 3),
      hdi_high = signif(bayestestR::hdi(var, ci = 0.95)[[3]], digits = 3),
      pd = signif(bayestestR::p_direction(var, method = "direct", null = 1)[[2]], digits = 3)
    )
  } else {
    vms <- c(paste0("intCond", rel, ":timepoint"), "timepoint", paste0("intCond", rel))
    vars <- brms::fixef(mod, summary = FALSE)[, vms] |> exp()
    out <- sapply(
      seq_len(ncol(vars)),
      function(i) {
        c(
          mean = signif(mean(vars[, i]), digits = 3),
          hdi_low = signif(bayestestR::hdi(vars[, i], ci = 0.95)[[2]], digits = 3),
          hdi_high = signif(bayestestR::hdi(vars[, i], ci = 0.95)[[3]], digits = 3),
          pd = signif(bayestestR::p_direction(vars[, i], method = "direct", null = 1)[[2]], digits = 3)
        )
      }
    )
    colnames(out) <- vms
  }
  out
}

get_hdi_pd <- function(draws, var) {
  idx <- draws[[var]]
  mean_idx <- mean(idx)
  hdi_idx <- paste0(
    "[",
    signif(bayestestR::hdi(idx, ci = 0.95)[[2]], 4),
    ", ",
    signif(bayestestR::hdi(idx, ci = 0.95)[[3]], 4),
    "]"
  )
  pd_idx <- bayestestR::p_direction(idx, method = "direct")[[2]]
  print(
    paste0(
      var, ": mean = ", signif(mean_idx, 3),
      ", 95\\% HDI = ", hdi_idx,
      ", p_direction = ", signif(pd_idx, 3)
    )
  )
}

## Function to analyze post-task ratings with zero-one-inflated beta regression
# analyze_post_ratings <- function(data, question_type_filter) {
#   # Filter data
#   data_filtered <- data |> dplyr::filter(question_type == question_type_filter)

#   # Fit model
#   model <- brms::brm(
#     brms::bf(
#       rating ~ intCond * timepoint + block_no + (1 + block_no | subID),
#       phi ~ (1 | subID),
#       zoi ~ intCond * timepoint + block_no + (1 | subID),
#       coi ~ intCond * timepoint + block_no + (1 | subID)
#     ),
#     family = brms::zero_one_inflated_beta(link = "logit", link_phi = "log"),
#     data = data_filtered,
#     warmup = 2000,
#     iter = 12000,
#     cores = 4,
#     backend = "cmdstanr"
#   )
#   # Get fixed effects
#   fixed_effects <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)

#   # Get interaction effect draws
#   interaction_effect <- brms::as_draws_df(model) |>
#     dplyr::select(dplyr::matches("b_.*intCond.*timepoint|b_.*timepoint.*intCond"))

#   # Get conditional effects plot
#   cond_effects_plot <- brms::conditional_effects(
#     model,
#     effects = "timepoint:intCond",
#     re_formula = NA
#   )

#   # Create newdata for both conditions at timepoint 1
#   newdata_compare <- data.frame(
#     intCond = factor(c("CR", "BA"), levels = c("CR", "BA")),
#     timepoint = 1,
#     block_no = mean(data_filtered$block_no)
#   )

#   # Extract zoi for both conditions
#   zoi_both <- brms::posterior_epred(
#     model,
#     newdata = newdata_compare,
#     dpar = "zoi",
#     re_formula = NA
#   )

#   # Extract coi for both conditions
#   coi_both <- brms::posterior_epred(
#     model,
#     newdata = newdata_compare,
#     dpar = "coi",
#     re_formula = NA
#   )

#   # Calculate one-inflation and zero-inflation for each condition
#   # Column 1 = CR, Column 2 = BA
#   one_inflation_cr <- zoi_both[, 1] * coi_both[, 1]
#   one_inflation_ba <- zoi_both[, 2] * coi_both[, 2]
#   zero_inflation_cr <- zoi_both[, 1] * (1 - coi_both[, 1])
#   zero_inflation_ba <- zoi_both[, 2] * (1 - coi_both[, 2])

#   # Create tibbles with differences
#   one_inflation_diff <- tibble::tibble(
#     difference = one_inflation_ba - one_inflation_cr
#   )
#   zero_inflation_diff <- tibble::tibble(
#     difference = zero_inflation_ba - zero_inflation_cr
#   )

#   # Return results
#   list(
#     model = model,
#     fixed_effects = fixed_effects,
#     interaction_effect = interaction_effect,
#     conditional_effects_plot = cond_effects_plot,
#     one_inflation_diff = one_inflation_diff,
#     zero_inflation_diff = zero_inflation_diff
#   )
# }

## Questionnaire ANCOVA: baseline-adjusted arm effects =========================
## Ported from the sibling paper repo (causality_training_mh/analyses/model_fns.R) and adapted to
## this repo's BA/CR arms. Replaces the earlier change-score analysis
## (paired t-test + lm(delta ~ group)): a baseline-adjusted model is more efficient whenever
## test-retest r < 1, and on this data PHQ-9's r is 0.41, so change-score residual variance is
## 2(1-r) = 1.18 against ANCOVA's (1-r^2) = 0.83.
##
## Deliberate deviation from the reference implementation: these take `ref`/`comp` arguments
## rather than repo-wide ARM_REF/ARM_INT constants, because the two task models use OPPOSITE
## condition coding (prep_data(): effort BA = 1, attribution CR = 1) and a global constant would
## invite a silent sign error when these are reused alongside them. The questionnaire convention
## is comp - ref = BA - CR, matching the earlier change-score analysis's compute_interaction().
##
## Fit on RAW scores: the data are never transformed, and cross-scale comparability is a reporting
## step (divide the contrast draws by sd(t2)) rather than something baked into the likelihood.

# shared interval widths, matching plot_change_slabs()/mini_interval_plot()'s mean_hdci convention
ANCOVA_WIDTHS <- c(0.95, 0.99)

# per-measure constants, ALWAYS computed on the primary (ITT) frame and reused unchanged for every
# other frame (PP, cLDA) so all frames sit on one reference grid and their estimates are directly
# comparable. `sd_t2` doubles as the standardising divisor. NOT sd(t1): baseline is range-restricted
# by the study's own screening criterion (PHQ-9 sd 3.29 at t1 vs 3.86 at t2), so sd(t1) is not a
# population sd and dividing by it inflates the PHQ-9 effect.
ancova_constants <- function(dat, prefix) {
  t1 <- dat[[paste0(prefix, "_total_t1")]]
  t2 <- dat[[paste0(prefix, "_total_t2")]]
  keep <- !is.na(t1) & !is.na(t2)
  stopifnot(sum(keep) > 1)
  sd_t1 <- stats::sd(t1[keep])
  sd_t2 <- stats::sd(t2[keep])
  list(
    t1_mean = mean(t1[keep]), t2_mean = mean(t2[keep]),
    sd_t1 = sd_t1, sd_t2 = sd_t2,
    # d_av denominator, computed ONCE per measure over both arms rather than within each -- a
    # per-arm denominator would put the arms on different scales, and the panel that uses this
    # exists to be compared across arms by eye. descriptive only; sd_t2 stays the model divisor.
    sd_av = (sd_t1 + sd_t2) / 2,
    r_t1t2 = stats::cor(t1[keep], t2[keep]),
    n = sum(keep)
  )
}

# modelling frame for one measure x analysis frame. complete cases only (an ANCOVA needs both
# timepoints); `t1_mean` comes from ancova_constants() on the primary frame, never recomputed here.
make_ancova_df <- function(dat, prefix, t1_mean, ref = "CR", comp = "BA", extra_cols = NULL) {
  d <- data.frame(
    subID = dat$subID,
    group = factor(as.character(dat$group), levels = c(ref, comp)),
    t1    = dat[[paste0(prefix, "_total_t1")]],
    t2    = dat[[paste0(prefix, "_total_t2")]]
  )
  for (nm in extra_cols) d[[nm]] <- dat[[nm]]
  d <- d[stats::complete.cases(d), ]
  d$t1_c <- d$t1 - t1_mean
  stopifnot(nrow(d) > 0, levels(d$group)[1] == ref, comp %in% levels(d$group))
  d
}

# Long (one row per observed score) frame for the cLDA missing-data sensitivity analysis. The
# critical difference from make_ancova_df(): it does NOT drop participants missing t2 -- that is the
# entire point, since those are the people a complete-case ANCOVA cannot see. Everyone has a t1
# score, so all randomised participants contribute at least one row.
#
# `post` is deliberately NUMERIC 0/1 rather than a factor. With `~ post + post:group`, R generates
# exactly one group column (post:group{COMP}) which is identically 0 at baseline, imposing the cLDA
# constraint that the arms share a baseline mean (justified by randomisation). The seemingly
# equivalent factor form `~ time + time:group` does NOT: it expands to time{t1}:group{COMP} AND
# time{t2}:group{COMP}, restoring a free baseline arm difference and silently destroying the
# constraint. The model still fits and still reports a plausible number, so this is asserted at the
# call site, not trusted.
make_clda_df <- function(dat, prefix, ref = "CR", comp = "BA") {
  base <- data.frame(
    subID = dat$subID,
    group = factor(as.character(dat$group), levels = c(ref, comp))
  )
  d <- rbind(
    cbind(base, post = 0L, score = dat[[paste0(prefix, "_total_t1")]]),
    cbind(base, post = 1L, score = dat[[paste0(prefix, "_total_t2")]])
  )
  d <- d[!is.na(d$score), ]
  # factor twin of `post`. REQUIRED, not decorative: the cLDA's unstructured covariance
  # (`unstr(time = tf, gr = subID)`) and its per-timepoint residual sd (`sigma ~ 0 + tf`) both key
  # on it, as does the nlme::gls cross-check's varIdent(~ 1 | tf).
  d$tf <- factor(d$post, levels = c(0L, 1L), labels = c("t1", "t2"))
  stopifnot(nrow(d) > 0, levels(d$group)[1] == ref, all(d$post %in% 0:1))
  d
}

# Posterior draws for a SUM of named model coefficients. Every quantity these analyses report is a
# linear combination of one or two named coefficients, because each model has at most a two-level
# factor and its interaction. Taken from the draws directly rather than via an emmeans reference
# grid: fewer moving parts, no draw-order shuffling, and every reported number traces to a
# coefficient that can be named in the write-up. Correctness is checked against lm()/gls() -- an
# independent implementation -- rather than against another layer of the same fit.
coef_draws <- function(fit, ...) {
  nms <- c(...)
  m <- posterior::as_draws_matrix(fit$fit)
  stopifnot(all(nms %in% colnames(m)))
  as.vector(rowSums(m[, nms, drop = FALSE]))
}

# the arm contrast: `group` is a two-level factor with `ref` as reference, so the comp - ref effect
# IS its coefficient. list-with-$contrast shape retained so the summarise/forest call sites read the
# same as the reference implementation's.
ancova_contrast_draws <- function(fit, comp = "BA") {
  v <- coef_draws(fit, paste0("b_group", comp))
  list(contrast = tibble::tibble(.draw = seq_along(v), .value = v))
}

# one tidy row per fit: raw-point effect (primary reporting), plus the two standardised variants so
# the denominator choice is visible rather than hidden. d_sd_t2 divides by a CONSTANT (pure
# relabelling of the posterior); d_sigma divides draw-wise by residual sigma, which the baseline
# adjustment shrinks, so it runs larger and is NOT what is reported. sd(t2) is the reported one
# (Cochrane's ANCOVA guidance). Point estimate is the MEAN and intervals are HDIs, matching the
# mean_hdci convention the rest of this file's plotting layer uses -- so the numbers here are the
# same summary the figures draw, not a different one.
summarise_ancova_draws <- function(fit, contrast_draws, const, measure, frame, n_obs) {
  v <- contrast_draws$.value
  sg <- as.vector(posterior::as_draws_matrix(fit$fit)[, "sigma"])
  h95 <- bayestestR::hdi(v, ci = 0.95)
  h90 <- bayestestR::hdi(v, ci = 0.90)
  d95 <- bayestestR::hdi(v / const$sd_t2, ci = 0.95)

  np <- brms::nuts_params(fit)
  dg <- posterior::summarise_draws(
    posterior::subset_draws(posterior::as_draws_df(fit$fit), variable = c("b_", "sigma"), regex = TRUE),
    "rhat", "ess_bulk", "ess_tail"
  )

  tibble::tibble(
    measure = measure, frame = frame, n = n_obs, sd_t2 = const$sd_t2, r_t1t2 = const$r_t1t2,
    est = mean(v), lo95 = h95$CI_low, hi95 = h95$CI_high,
    lo90 = h90$CI_low, hi90 = h90$CI_high,
    d_sd_t2 = mean(v) / const$sd_t2, d_lo95 = d95$CI_low, d_hi95 = d95$CI_high,
    d_sigma = mean(v / sg),
    pd = as.numeric(bayestestR::p_direction(v, method = "direct")),
    rhat_max = max(dg$rhat), ess_bulk_min = min(dg$ess_bulk), ess_tail_min = min(dg$ess_tail),
    divergent = sum(np$Value[np$Parameter == "divergent__"])
  )
}

# lighter summariser for when two DIFFERENT models feed one table (ANCOVA vs cLDA, primary vs a
# sensitivity refit): takes draws rather than a fit, so it reports no diagnostics of its own.
summarise_effect_draws <- function(v, sd_t2) {
  h95 <- bayestestR::hdi(v, ci = 0.95)
  tibble::tibble(
    est = mean(v), lo95 = h95$CI_low, hi95 = h95$CI_high,
    d_sd_t2 = mean(v) / sd_t2, sd_t2 = sd_t2,
    pd = as.numeric(bayestestR::p_direction(v, method = "direct"))
  )
}

# worst-case sampler diagnostics across a brms fit's population-level and scale parameters
brms_diag <- function(fit, vars = c("b_", "sigma", "sd_")) {
  np <- brms::nuts_params(fit)
  dg <- posterior::summarise_draws(
    posterior::subset_draws(posterior::as_draws_df(fit$fit), variable = vars, regex = TRUE),
    "rhat", "ess_bulk", "ess_tail"
  )
  tibble::tibble(
    rhat_max = max(dg$rhat), ess_bulk_min = min(dg$ess_bulk), ess_tail_min = min(dg$ess_tail),
    divergent = sum(np$Value[np$Parameter == "divergent__"])
  )
}

# Forest of standardised arm-effect posteriors, one row per measure. Bayesian replacement for
# the earlier change-score analysis's plot_forest() (which drew a frequentist change-score Cohen's d
# with significance stars). Follows this file's plotting conventions: mean_hdci intervals at
# ANCOVA_WIDTHS, a slab+interval for a single series, and a plain pointinterval where rows are dodged
# (slabs from two dodged series overlap illegibly) -- the same split plot_change_slabs() already makes
# between its single- and multi-series callers.
#   draws_df: long -- measure (factor; level order is BOTTOM-to-top, so pass rev(main_measures)),
#             .value (already standardised by the caller), plus `dodge_by`'s column if used.
#   dodge_by: optional column dodged within each measure row (e.g. ANCOVA vs cLDA).
#   colour_by: defaults to the measure, which is right when every dodged series is the same quantity
#             under a different frame. Pass the dodge column instead when the series ARE the thing
#             being compared, and the colour scale then becomes a shown legend.
#   secondary_measures: optional measure names to separate below a horizontal rule (an outcome-tier
#             split). NULL -- the default -- draws no rule; this repo has no pre-registered outcome
#             hierarchy in code, and inventing one would be a substantive claim, not a style choice.
plot_ancova_forest <- function(draws_df, labels, measure_cols, dodge_by = NULL, tag = "B",
                              xlab = NULL, base_size = 18, dodge_w = 0.55,
                              colour_by = NULL, colour_values = NULL,
                              dodge_shapes = NULL, dodge_alphas = NULL, dodge_labels = NULL,
                              secondary_measures = NULL, tier_labels = c("primary", "secondary"),
                              legend_pos = "bottom") {
  clr <- colour_by %||% "measure"
  clr_vals <- colour_values %||% measure_cols
  # when the dodge IS the colour, the colour scale has to be a shown guide so the keys carry the
  # series colour, and must share the shape/alpha scales' name and labels so ggplot merges the three
  # into one legend rather than stacking three copies. when colour tracks the measure instead it
  # stays hidden -- the measure is already the y axis.
  clr_in_legend <- !is.null(dodge_by) && identical(clr, dodge_by)

  p <- ggplot2::ggplot(draws_df, ggplot2::aes(x = .value, y = measure, colour = .data[[clr]])) +
    ggplot2::geom_vline(xintercept = 0, linetype = "32", colour = "grey50")

  if (is.null(dodge_by)) {
    p <- p +
      ggdist::stat_slabinterval(
        ggplot2::aes(fill = .data[[clr]]),
        .width = ANCOVA_WIDTHS, point_interval = "mean_hdci", slab_alpha = 0.45, scale = 0.7,
        point_size = 3, interval_size_range = c(1, 3), fatten_point = 2
      ) +
      ggplot2::scale_fill_manual(values = clr_vals, guide = "none")
  } else {
    n_lvl <- nlevels(draws_df[[dodge_by]])
    stopifnot(n_lvl >= 2)
    shapes <- dodge_shapes %||% c(16, 21, 15, 17, 18)[seq_len(n_lvl)]
    alphas <- dodge_alphas %||% if (n_lvl == 2) c(1, 0.6) else seq(1, 0.65, length.out = n_lvl)
    labs_d <- dodge_labels %||% levels(draws_df[[dodge_by]])
    # reverse = TRUE so the FIRST level sits on top within each measure row, matching the
    # primary-first ordering used everywhere else
    p <- p +
      ggdist::stat_pointinterval(
        ggplot2::aes(shape = .data[[dodge_by]], alpha = .data[[dodge_by]]),
        .width = ANCOVA_WIDTHS, point_interval = "mean_hdci",
        point_size = 3, interval_size_range = c(1, 3), fatten_point = 2,
        position = ggplot2::position_dodge(width = dodge_w, reverse = TRUE)
      ) +
      ggplot2::scale_shape_manual(values = shapes, name = NULL, labels = labs_d) +
      ggplot2::scale_alpha_manual(values = alphas, name = NULL, labels = labs_d)
  }

  # optional outcome-tier rule, between the secondary block (bottom) and the primary block (top).
  # y levels run bottom-to-top, so the rule sits above however many secondary rows there are.
  if (!is.null(secondary_measures)) {
    stopifnot(length(tier_labels) == 2, all(secondary_measures %in% levels(draws_df$measure)))
    n_sec <- length(secondary_measures)
    tier_ann <- function(lbl, y, vjust) {
      ggplot2::annotate(
        "text", x = Inf, y = y, label = lbl, hjust = 1.05, vjust = vjust,
        family = "Open Sans", size = base_size / 4.6, colour = "grey45", fontface = "italic"
      )
    }
    p <- p +
      ggplot2::geom_hline(yintercept = n_sec + 0.5, colour = "grey75", linewidth = 0.5) +
      tier_ann(tier_labels[[1]], n_sec + 0.62, 0) +
      tier_ann(tier_labels[[2]], n_sec + 0.38, 1)
  }

  p <- p +
    ggplot2::scale_y_discrete(labels = labels, expand = ggplot2::expansion(add = 0.7)) +
    ggplot2::labs(x = xlab, y = NULL, tag = tag) +
    ggplot2::coord_cartesian(clip = "off")

  p <- p + if (clr_in_legend) {
    ggplot2::scale_colour_manual(values = clr_vals, name = NULL, labels = dodge_labels %||% ggplot2::waiver())
  } else {
    ggplot2::scale_colour_manual(values = clr_vals, guide = "none")
  }

  p +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = base_size) +
    ggplot2::theme(
      plot.tag     = ggplot2::element_text(face = "bold", size = base_size + 6, colour = "black"),
      axis.title.x = ggtext::element_markdown(size = base_size - 2, lineheight = 1.15),
      axis.text.y  = ggtext::element_markdown(),
      legend.position = if (is.null(dodge_by)) "none" else legend_pos,
      legend.justification.bottom = "center",
      legend.background = ggplot2::element_rect(fill = "white", colour = "#cacaca", linewidth = 0.5),
      legend.margin = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
      plot.margin  = ggplot2::margin(2, 30, 2, 10)
    )
}

# Per-participant standardised change, all measures on one effect-size axis, dodged by arm.
# Condenses the six per-measure pre->post slope panels (plot_qnr_measure()) into one panel: same
# descriptive result, same visual grammar as those panels (quasirandom individual points + mean and
# 95% CI) so it still reads as OBSERVED data rather than as a model posterior -- the forest beside it
# is the model-estimated one. Coloured by ARM, and supplies the figure's arm legend.
#   df: long -- measure (factor, plotting order left->right), group (factor, `ref` first),
#       .value (already standardised BY THE CALLER; the panel is agnostic about the divisor and the
#       axis label is where the units are declared).
#   raw_change: optional tibble(measure, label) for a single pooled (both arms) headline annotation
#       per measure -- raw-score change plus d_av, so a reader can judge clinical magnitude on the
#       measure's own scale. Deliberately NOT per arm: that would double the value labels already on
#       the panel. Positioned above each measure's OWN labels rather than at a shared panel-wide row.
plot_change_distributions <- function(df, labels, tag = "A", ylab = NULL, base_size = 18,
                                      grp_cols, grp_labs, grp_shapes, dodge_w = 0.6, note = NULL,
                                      note_pos = c("tr", "br"),
                                      show_values = FALSE, value_digits = 2, value_nudge = 0.28,
                                      value_size = base_size / 3.6,
                                      raw_change = NULL, raw_change_nudge = 0.55,
                                      raw_change_size = base_size / 4.7,
                                      top_expand = 0.14) {
  note_pos <- match.arg(note_pos)
  dsumm <- df |>
    dplyr::group_by(measure, group) |>
    dplyr::summarise(
      mean_chg = mean(.value), n = dplyr::n(), sd_chg = stats::sd(.value),
      se = sd_chg / sqrt(n), ci = stats::qt(0.975, n - 1) * se, .groups = "drop"
    )

  pd <- ggplot2::position_dodge(width = dodge_w)

  p <- ggplot2::ggplot(dsumm, ggplot2::aes(x = measure, y = mean_chg, colour = group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "32", colour = "grey50") +
    ggplot2::geom_point(
      data = df, ggplot2::aes(x = measure, y = .value, colour = group, shape = group),
      position = ggbeeswarm::position_quasirandom(width = 0.13, dodge.width = dodge_w),
      alpha = 0.34, size = 1.3, show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_chg - ci, ymax = mean_chg + ci),
      position = pd, width = 0.14, linewidth = 0.9
    ) +
    ggplot2::geom_point(ggplot2::aes(shape = group), position = pd, size = 3) +
    ggplot2::scale_colour_manual(values = grp_cols, labels = grp_labs, name = NULL) +
    ggplot2::scale_shape_manual(values = grp_shapes, labels = grp_labs, name = NULL) +
    ggplot2::scale_x_discrete(labels = labels) +
    # headroom at the top for the value / raw-change labels, plus extra on whichever side `note`
    # sits; the corner it is given is empty in the data either way
    # `top_expand` must keep pace with value_nudge/raw_change_nudge: the labels are placed at a
    # data value, so raising them without raising the panel ceiling clips them
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(if (!is.null(note) && note_pos == "br") 0.12 else 0.05, top_expand)
      )
    ) +
    ggplot2::labs(x = NULL, y = ylab, tag = tag)

  # numeric effect size per measure x arm, taken from `dsumm` -- i.e. from the SAME summary that
  # positions the markers, so the printed number cannot drift from the plotted point. value only:
  # the 95% CI is already drawn as the errorbar. semi-transparent fill so the label reads over the
  # beeswarm without punching a hole in it.
  if (show_values) {
    p <- p + ggtext::geom_richtext(
      data = dsumm,
      ggplot2::aes(
        x = measure, y = mean_chg + ci + value_nudge, colour = group,
        label = formatC(mean_chg, format = "f", digits = value_digits, flag = "+")
      ),
      position = pd, inherit.aes = FALSE, show.legend = FALSE,
      size = value_size, family = "Open Sans", lineheight = 1.2,
      fill = grDevices::adjustcolor("white", alpha.f = 0.72), label.color = NA,
      label.padding = grid::unit(c(0.06, 0.12, 0.06, 0.12), "lines"),
      label.r = grid::unit(0.1, "lines")
    )
  }

  if (!is.null(raw_change)) {
    # sits above the taller of the two group value labels (if shown) or the taller CI otherwise, so
    # it never collides with either -- computed from `dsumm`, not hardcoded, so it tracks whatever
    # headroom show_values/value_nudge already carved out
    rc_y <- dsumm |>
      dplyr::group_by(measure) |>
      dplyr::summarise(y = max(mean_chg + ci), .groups = "drop") |>
      dplyr::mutate(y = y + (if (show_values) value_nudge else 0) + raw_change_nudge)
    rc <- dplyr::left_join(raw_change, rc_y, by = "measure")
    stopifnot(!anyNA(rc$y))

    p <- p + ggtext::geom_richtext(
      data = rc, ggplot2::aes(x = measure, y = y, label = label), inherit.aes = FALSE,
      size = raw_change_size, family = "Open Sans", fontface = "italic", colour = "grey30",
      lineheight = 1.15,
      fill = grDevices::adjustcolor("white", alpha.f = 0.72), label.color = NA,
      label.padding = grid::unit(c(0.06, 0.12, 0.06, 0.12), "lines"),
      label.r = grid::unit(0.1, "lines")
    )
  }

  if (!is.null(note)) {
    p <- p + ggtext::geom_richtext(
      data = data.frame(x = Inf, y = if (note_pos == "br") -Inf else Inf, label = note),
      ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
      hjust = 1, vjust = if (note_pos == "br") -0.5 else 1.5, fill = NA, label.color = "slategrey",
      label.padding = grid::unit(rep(0.2, 4), "lines"), label.r = grid::unit(0.4, "lines"),
      colour = "grey20", lineheight = 1.2, size = base_size / 4, family = "Open Sans"
    )
  }

  p +
    cowplot::theme_half_open(font_family = "Open Sans", font_size = base_size) +
    cowplot::background_grid(major = "y", minor = "none") +
    ggplot2::theme(
      plot.tag        = ggplot2::element_text(face = "bold", size = base_size + 6, colour = "black"),
      axis.title.y    = ggtext::element_markdown(size = base_size - 2),
      axis.text.x     = ggtext::element_markdown(),
      legend.position = "bottom",
      legend.justification.bottom = "center",
      legend.background = ggplot2::element_rect(fill = "white", colour = "#cacaca", linewidth = 0.5),
      legend.margin   = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm"),
      plot.margin     = ggplot2::margin(2, 10, 2, 10)
    )
}

## Test-retest, covariance and mediator collinearity ===========================
## Ported from the sibling repo. Every model in this repo builds its individual-level parameters as
## diag_pre_multiply(sigma, R_chol) * raw and emits only the CORRELATION matrix in generated
## quantities, so a covariance has to be reconstructed in R from that plus the sd vector -- hence
## get_cov_matrix_tbl() alongside get_corr_matrix_tbl().

# tidy (mean, 95% HDI, pd) summary of every unique off-diagonal entry of a k x k correlation-matrix
# parameter's draws (e.g. var = "R_theta_neg", with columns "R_theta_neg[i,j]" in `draws`). `labels`
# names the k row/col indices. One row per unique pair (upper triangle, i < j) -- the diagonal
# (always 1) is dropped.
get_corr_matrix_tbl <- function(draws, var, labels) {
  pairs <- utils::combn(length(labels), 2)
  rows <- lapply(seq_len(ncol(pairs)), function(p) {
    i <- pairs[1, p]
    j <- pairs[2, p]
    nm <- sprintf("%s[%d,%d]", var, i, j)
    if (!nm %in% names(draws)) stop("get_corr_matrix_tbl(): draws are missing column '", nm, "'")
    v <- draws[[nm]]
    hdi <- bayestestR::hdi(v, ci = 0.95)
    tibble::tibble(
      var = var, label_i = labels[i], label_j = labels[j],
      mean = mean(v), hdi_95_lo = hdi$CI_low, hdi_95_hi = hdi$CI_high,
      pd = as.numeric(bayestestR::p_direction(v, method = "direct"))
    )
  })
  dplyr::bind_rows(rows)
}

# covariance counterpart: cov[i,j] = R[i,j] * sigma[i] * sigma[j], formed PER DRAW so the summary
# carries the joint uncertainty of the correlation and both sds. `sd_var` is the vector-valued sd
# parameter indexed on the SAME ordering as the correlation matrix (e.g. corr_var = "R_theta_neg"
# with sd_var = "pars_sigma_neg", both ordered internal_t1, global_t1, internal_t2, global_t2).
# Unlike the correlation table the DIAGONAL is kept: those are the variances, which carry information
# where a correlation's diagonal of 1 does not.
get_cov_matrix_tbl <- function(draws, corr_var, sd_var, labels) {
  k <- length(labels)
  sd_cols <- sprintf("%s[%d]", sd_var, seq_len(k))
  miss <- setdiff(sd_cols, names(draws))
  if (length(miss)) stop("get_cov_matrix_tbl(): draws are missing column(s) ", paste(miss, collapse = ", "))
  idx <- which(upper.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
  rows <- lapply(seq_len(nrow(idx)), function(p) {
    i <- idx[p, "row"]
    j <- idx[p, "col"]
    r <- if (i == j) 1 else {
      nm <- sprintf("%s[%d,%d]", corr_var, i, j)
      if (!nm %in% names(draws)) stop("get_cov_matrix_tbl(): draws are missing column '", nm, "'")
      draws[[nm]]
    }
    v <- r * draws[[sd_cols[i]]] * draws[[sd_cols[j]]]
    hdi <- bayestestR::hdi(v, ci = 0.95)
    tibble::tibble(
      var = corr_var, label_i = labels[i], label_j = labels[j],
      quantity = if (i == j) "variance" else "covariance",
      mean = mean(v), hdi_95_lo = hdi$CI_low, hdi_95_hi = hdi$CI_high
    )
  })
  dplyr::bind_rows(rows)
}

# raw per-draw correlation (length = ndraws) between two per-participant vector-valued generated
# quantities (e.g. var_x = "delta_rewSens_p", columns "delta_rewSens_p[id]" in `draws`), computed PER
# DRAW across participants -- so shrinkage/measurement-error uncertainty in both quantities
# propagates into the resulting distribution, unlike naively correlating posterior means, which
# understates it. Sum-of-squares identity (same trick as get_symptom_corrs()'s corr_from_sums) keeps
# this vectorised over draws. `ids` (Stan participant indices, NULL = everyone) restricts to a
# subset, e.g. one arm -- pooling arms would let between-arm mean differences leak into what is meant
# to be a within-person association.
get_draws_corr_raw <- function(draws, var_x, var_y, ids = NULL) {
  extract_sorted <- function(v) {
    cols <- grep(paste0("^", v, "\\[\\d+\\]$"), names(draws), value = TRUE)
    if (length(cols) == 0) stop("get_draws_corr_raw(): no columns matching '", v, "[...]'")
    col_ids <- as.integer(sub(".*\\[(\\d+)\\]$", "\\1", cols))
    ord <- order(col_ids)
    cols <- cols[ord]
    col_ids <- col_ids[ord]
    if (!is.null(ids)) cols <- cols[col_ids %in% ids]
    as.matrix(draws[, cols])
  }
  x <- extract_sorted(var_x)
  y <- extract_sorted(var_y)
  stopifnot(ncol(x) == ncol(y), ncol(x) > 1)
  n <- ncol(x)
  sx <- rowSums(x)
  sy <- rowSums(y)
  sxx <- rowSums(x^2)
  syy <- rowSums(y^2)
  sxy <- rowSums(x * y)
  (n * sxy - sx * sy) / sqrt((n * sxx - sx^2) * (n * syy - sy^2))
}

# mean/95% HDI/pd summary of a raw per-draw statistic (e.g. from get_draws_corr_raw(), or the
# elementwise difference of two such vectors, to test whether an association differs credibly
# between two subsets).
summarise_draws_r <- function(r) {
  hdi <- bayestestR::hdi(r, ci = 0.95)
  tibble::tibble(
    mean = mean(r), hdi_95_lo = hdi$CI_low, hdi_95_hi = hdi$CI_high,
    pd = as.numeric(bayestestR::p_direction(r, method = "direct"))
  )
}

get_draws_corr <- function(draws, var_x, var_y, ids = NULL) {
  summarise_draws_r(get_draws_corr_raw(draws, var_x, var_y, ids = ids))
}

# Per-draw correlations between every pair of a set of per-participant delta quantities, within each
# arm and pooled. This is the collinearity check that decides whether a JOINT (mutually-adjusted)
# mediation model is the right call: the sibling repo abandoned its joint 4-mediator attribution
# model for one-fit-per-mediator precisely because its four delta-theta domains ran r ~ 0.8 on
# posterior means (0.49-0.70 per draw), which flipped b-path signs. `vars` is a NAMED vector of
# draws-column stems (names = display labels); `group_ids` a named list of Stan id vectors per arm.
mediator_collinearity <- function(draws, vars, group_ids = NULL, label = NA_character_) {
  stopifnot(length(vars) >= 2, !is.null(names(vars)))
  pairs <- utils::combn(length(vars), 2)
  sets <- c(list(pooled = NULL), if (is.null(group_ids)) list() else group_ids)
  dplyr::bind_rows(lapply(names(sets), function(g) {
    dplyr::bind_rows(lapply(seq_len(ncol(pairs)), function(p) {
      i <- pairs[1, p]
      j <- pairs[2, p]
      r <- get_draws_corr_raw(draws, vars[[i]], vars[[j]], ids = sets[[g]])
      n_used <- if (is.null(sets[[g]])) {
        # pooled: every participant the quantity is defined for
        length(grep(paste0("^", vars[[i]], "\\[\\d+\\]$"), names(draws)))
      } else {
        length(sets[[g]])
      }
      dplyr::bind_cols(
        tibble::tibble(
          comparison = label, group = g,
          par_x = names(vars)[i], par_y = names(vars)[j], n_ppts = n_used
        ),
        summarise_draws_r(r)
      )
    }))
  }))
}

# reliability (lambda) of an individual-level measurement from its per-participant posterior
# summaries: the share of between-participant variance that is signal rather than shrinkage/
# measurement noise. lambda = tau^2 / (tau^2 + sigma^2). Classical-test-theory reliability; identical
# in form to an ICC and to the empirical-Bayes shrinkage weight (Kelley's weight).
#
# estimating it from posterior summaries: for a normal-normal hierarchy the posterior mean is
# m_i = lambda*y_i + (1 - lambda)*mu with posterior variance v_i = lambda*sigma^2, so
#   var_i(m_i)  = lambda^2 * (tau^2 + sigma^2) = lambda * tau^2
#   mean_i(v_i) = lambda * sigma^2
# hence var(m)/mean(v) = tau^2/sigma^2 and lambda = ratio/(1 + ratio). the lambda factors cancel,
# which is why this works on shrunken summaries at all. it assumes lambda is not wildly heterogeneous
# across participants, so treat it as approximate.
reliability <- function(hat, se) {
  stopifnot(length(hat) == length(se), length(hat) > 1)
  ratio <- stats::var(hat) / mean(se^2)
  ratio / (1 + ratio)
}

# Tidy direct / indirect / total decomposition of a joint ANCOVA mediation fit, with a cross-check of
# the model's own total_arm against the mediator-free questionnaire ANCOVA.
#
# One row per mediator x component: the a-paths (ITT and observed-t2-only), the b-paths per arm AND
# their between-arm difference (`b_diff`), the per-arm indirect effects and their difference
# (`index_mod_med`), plus the model-level `direct_arm` and `total_arm`. Every contrast row is
# oriented to `arm_dir`.
#
# total_arm is the model's EXACT reconstruction of the marginal arm effect on the baseline-adjusted
# outcome (the baseline-mediator-imbalance terms included, not an approximation -- see the .stan
# header). direct_arm + sum(index_mod_med_*) is NOT guaranteed to equal it when baseline mediator
# level is arm-imbalanced, which it is here, so the two are reported separately.
#
# SIGN TRAP, asserted rather than assumed: `arm_dir` says which contrast this fit's condition coding
# produces. prep_data() codes condition == 1 as BA for type = "effort" but as CR for
# type = "attribution", so direct_arm/total_arm mean BA - CR in the reward-effort models and CR - BA
# in the causal-attribution ones, while the questionnaire ANCOVA forest is always BA - CR. The
# reference value is flipped here to match, rather than the comparison being made on trust.
mediation_direct_total <- function(draws, suffixes, mediator_nms,
                                   arm_dir = c("BA - CR", "CR - BA"),
                                   ancova_path = "outputs/ancova_primary.csv",
                                   measure = "PHQ-9", tol = 0.05, out_path = NULL) {
  arm_dir <- match.arg(arm_dir)
  stopifnot(!is.null(names(mediator_nms)), all(suffixes %in% names(mediator_nms)))
  need <- function(col) {
    if (!col %in% names(draws)) stop("mediation_direct_total(): draws are missing '", col, "'")
    draws[[col]]
  }
  row_of <- function(mediator, component, v) {
    h <- bayestestR::hdi(v, ci = 0.95)
    tibble::tibble(
      mediator = mediator, component = component, mean = mean(v),
      hdi_95_lo = h$CI_low, hdi_95_hi = h$CI_high,
      pd = as.numeric(bayestestR::p_direction(v, method = "direct"))
    )
  }
  # per-mediator paths. the a_*_obs rows are the observed-t2-only companions to the ITT a-paths:
  # where they diverge, the ITT quantity is leaning on model-imputed t2 values for non-completers.
  per_med <- dplyr::bind_rows(lapply(suffixes, function(sfx) {
    nm <- unname(mediator_nms[[sfx]])
    obs <- paste0("a_", sfx, c("_ba_obs", "_cr_obs"))
    b_ba <- need(paste0("b_", sfx, "_ba"))
    b_cr <- need(paste0("b_", sfx, "_cr"))
    # between-arm difference in the b path: how much the mediator -> symptom slope itself differs by
    # arm, which is the model's beta_*W_group up to sign. Taken DRAW-WISE, so it carries the full
    # posterior of the difference -- it is not recoverable from the two per-arm rows, whose intervals
    # overlap far more than the difference's does whenever the two share posterior uncertainty.
    #
    # Oriented to `arm_dir`, so every contrast row in this table (b_diff, index_mod_med, direct_arm,
    # total_arm) reads in the SAME direction. That is not a convention imposed here: the .stan files
    # already define index_mod_med as indirect_ba - indirect_cr in the reward-effort model and
    # indirect_cr - indirect_ba in the attribution one, matching each fit's condition coding. The
    # a_*/b_* rows are labelled by arm and so are unaffected either way.
    b_diff <- if (arm_dir == "BA - CR") b_ba - b_cr else b_cr - b_ba
    dplyr::bind_rows(
      row_of(nm, "a_ba", need(paste0("a_", sfx, "_ba"))),
      row_of(nm, "a_cr", need(paste0("a_", sfx, "_cr"))),
      if (all(obs %in% names(draws))) row_of(nm, "a_ba_obs", draws[[obs[1]]]),
      if (all(obs %in% names(draws))) row_of(nm, "a_cr_obs", draws[[obs[2]]]),
      row_of(nm, "b_ba", b_ba),
      row_of(nm, "b_cr", b_cr),
      row_of(nm, "b_diff", b_diff),
      row_of(nm, "indirect_ba", need(paste0("indirect_", sfx, "_ba"))),
      row_of(nm, "indirect_cr", need(paste0("indirect_", sfx, "_cr"))),
      row_of(nm, "index_mod_med", need(paste0("index_mod_med_", sfx)))
    )
  }))
  # direct and total belong to the whole model, not to a mediator
  shared <- dplyr::bind_rows(
    row_of("(all mediators)", "direct_arm", need("direct_arm")),
    row_of("(all mediators)", "total_arm", need("total_arm"))
  )
  tbl <- dplyr::bind_rows(per_med, shared) |>
    dplyr::mutate(arm_dir = arm_dir, .after = component)

  # cross-check total_arm against the mediator-free ANCOVA, if that table has been written
  gap <- NA_real_
  if (file.exists(ancova_path)) {
    ref <- readr::read_csv(ancova_path, show_col_types = FALSE) |>
      dplyr::filter(measure == !!measure)
    if (nrow(ref) == 1) {
      ref_val <- if (arm_dir == "BA - CR") ref$d_sd_t2 else -ref$d_sd_t2
      gap <- mean(need("total_arm")) - ref_val
      cat(sprintf(
        "\n  total_arm (%s) = %+.3f vs questionnaire ANCOVA %s b_std = %+.3f  -> gap %+.3f SD\n",
        arm_dir, mean(need("total_arm")), measure, ref_val, gap
      ))
      cat(sprintf(
        "  (ANCOVA n = %d on sd_t2 = %.3f; the mediation frame may differ slightly after the task\n",
        ref$n, ref$sd_t2
      ))
      cat("   catch exclusions, so a small gap is expected -- a large one is not)\n")
      if (abs(gap) > tol) {
        warning(
          "total_arm differs from the questionnaire ANCOVA by ", sprintf("%+.3f", gap),
          " SD (> ", tol, "): check the symptom standardisation (z_sd_ref) and the arm direction",
          call. = FALSE
        )
      }
    }
  } else {
    message("  ", ancova_path, " not found -- skipping the total_arm cross-check")
  }
  tbl$gap_vs_ancova <- ifelse(tbl$component == "total_arm", gap, NA_real_)

  if (!is.null(out_path)) readr::write_csv(tbl, out_path)
  tbl
}

# Per-draw correlation between two per-participant quantities that live in DIFFERENT fits (e.g. the
# reward-effort model's delta_rewSens_p against the causal-attribution model's
# delta_internal_neg_p). Needed because prep_data() numbers participants from each task's own long
# file, so the two `id` vectors are NOT the same people -- the alignment here is on subID, and doing
# it on id would silently correlate one participant's reward change with another's theta change.
#
# `ids_x` / `ids_y` are the prep_data(..., ret_ids = TRUE) frames for the two fits (columns id,
# subID). Draws are paired position-wise across the two fits: they are independent posteriors, so the
# joint factorises and pairing draw j with draw j gives a valid posterior for the correlation. If the
# two have different numbers of draws the longer is truncated, which is reported.
get_cross_task_corr_raw <- function(draws_x, var_x, ids_x, draws_y, var_y, ids_y,
                                    sub_ids = NULL, quiet = FALSE) {
  pull_mat <- function(draws, var, ids) {
    cols <- grep(paste0("^", var, "\\[\\d+\\]$"), names(draws), value = TRUE)
    if (length(cols) == 0) stop("get_cross_task_corr_raw(): no columns matching '", var, "[...]'")
    col_ids <- as.integer(sub(".*\\[(\\d+)\\]$", "\\1", cols))
    sub <- ids$subID[match(col_ids, ids$id)]
    if (anyNA(sub)) stop("get_cross_task_corr_raw(): '", var, "' has ids absent from its ids frame")
    m <- as.matrix(draws[, cols])
    colnames(m) <- sub
    m
  }
  mx <- pull_mat(draws_x, var_x, ids_x)
  my <- pull_mat(draws_y, var_y, ids_y)
  common <- intersect(colnames(mx), colnames(my))
  if (!is.null(sub_ids)) common <- intersect(common, sub_ids)
  if (length(common) < 3) stop("get_cross_task_corr_raw(): fewer than 3 shared participants")
  mx <- mx[, common, drop = FALSE]
  my <- my[, common, drop = FALSE]
  nd <- min(nrow(mx), nrow(my))
  if (nrow(mx) != nrow(my) && !quiet) {
    message(sprintf(
      "  get_cross_task_corr_raw(): %d vs %d draws -- truncating both to %d", nrow(mx), nrow(my), nd
    ))
  }
  mx <- mx[seq_len(nd), , drop = FALSE]
  my <- my[seq_len(nd), , drop = FALSE]
  n <- length(common)
  sx <- rowSums(mx)
  sy <- rowSums(my)
  r <- (n * rowSums(mx * my) - sx * sy) /
    sqrt((n * rowSums(mx^2) - sx^2) * (n * rowSums(my^2) - sy^2))
  attr(r, "n_ppts") <- n
  r
}

# tidy wrapper: one row per (x, y) pair, with the number of shared participants carried through
cross_task_collinearity <- function(draws_x, vars_x, ids_x, draws_y, vars_y, ids_y,
                                    sub_ids = NULL, group = "pooled") {
  stopifnot(!is.null(names(vars_x)), !is.null(names(vars_y)))
  dplyr::bind_rows(lapply(names(vars_x), function(nx) {
    dplyr::bind_rows(lapply(names(vars_y), function(ny) {
      r <- get_cross_task_corr_raw(
        draws_x, vars_x[[nx]], ids_x, draws_y, vars_y[[ny]], ids_y, sub_ids = sub_ids, quiet = TRUE
      )
      dplyr::bind_cols(
        tibble::tibble(group = group, par_x = nx, par_y = ny, n_ppts = attr(r, "n_ppts")),
        summarise_draws_r(r)
      )
    }))
  }))
}
