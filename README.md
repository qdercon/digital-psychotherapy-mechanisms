# digital-psychotherapy-mechanisms

## Data and analysis code for a randomised online study of two digital psychotherapy modules

Participants were randomised to a behavioural activation (BA) or a cognitive restructuring (CR) module, and completed a reward-effort decision-making task, a causal attribution task and a set of questionnaires before and after the intervention.

### This study builds on:

> Norbury, A., Hauser, T. U., Fleming, S. M., Dolan, R. J., & Huys, Q. J. M. (2024). Different components of cognitive-behavioral therapy affect specific cognitive mechanisms. *Science Advances*, **10**(13), eadk3222. https://doi.org/10.1126/sciadv.adk3222

The tasks are adapted from those in the accompanying repository: [agnesnorbury/cognitive-mechanisms-psychotherapy](https://github.com/agnesnorbury/cognitive-mechanisms-psychotherapy).

## Repository structure

| Path | Contents |
| --- | --- |
| `STAND-prescreen/` | prescreening questionnaire (web code) |
| `STAND-tasks-1/` | session 1 (pre-intervention) tasks and questionnaires |
| `STAND-tasks-2/` | session 2 (post-intervention) tasks and questionnaires, the combined long-format task data, and the Stan models (`analysis/stan-models/`) |
| `model_fns.R` | all data preparation, modelling and plotting functions (sourced by every script) |
| `qnr_analysis_ancova.R` | questionnaire outcomes: baseline-adjusted (ANCOVA) arm effects (Figure 2) |
| `modelling_reward-effort.R` | reward-effort task: model-free and model-based analyses, mediation, sensitivity analyses and parameter recovery (Figures 3 and 5, Supplementary Figures 1–3) |
| `modelling_causal-attr.R` | causal attribution task: model-free and model-based analyses and mediation (Figures 4 and 5) |
| `task_covariance_collinearity.R` | test-retest reliability of task parameters, reliability of change scores, and mediator collinearity |

Each `public/` folder contains the web code for that session's study (built with [jsPsych](https://www.jspsych.org/) and [Phaser](https://phaser.io/)).

## Data

All data are provided as `.csv` files ending `_scrambled.csv`. Participant IDs have been replaced with random IDs (`sub_XXXXX`), consistent across all files, and platform identifiers have been removed. Within the analysis scripts, files are loaded with `read_stand_data()` (defined in `model_fns.R`).

| File | Contents |
| --- | --- |
| `STAND-tasks-{1,2}/analysis/STAND-tasks-{1,2}-self-report-data_scrambled.csv` | questionnaire responses at each session |
| `STAND-tasks-2/analysis/qqnrs_t2_scrambled.csv` | questionnaire totals at both sessions, one row per participant |
| `STAND-tasks-2/analysis/stand-tasks-both-rew-eff-task-data-long_scrambled.csv` | reward-effort task, trial-level, both sessions |
| `STAND-tasks-2/analysis/stand-tasks-both-causal-attr-task-data-long_scrambled.csv` | causal attribution task, trial-level, both sessions |
| `STAND-tasks-{1,2}/analysis/STAND-tasks-{1,2}-rew-eff-task-post-task-data-long_scrambled.csv` | post-block ratings from the reward-effort task |
| `STAND-tasks-2/analysis/worksheet_qc_scrambled.csv` | coded worksheet-quality ratings for intervention completers |

## Running the analyses

### Requirements

- [R](https://www.r-project.org/) (tested with 4.5.2)
- [CmdStan](https://mc-stan.org/cmdstanr/) via `cmdstanr` (tested with CmdStan 2.39.0), and [`brms`](https://paul-buerkner.github.io/brms/) for the questionnaire models
- the R packages used by the scripts:

```R
install.packages(c(
  "bayestestR", "brms", "broom.mixed", "cowplot", "data.table", "dplyr", "effectsize", "ggbeeswarm",
  "ggdist", "ggplot2", "ggpp", "ggtext", "lme4", "loo", "MetBrewer", "nlme", "patchwork", "posterior",
  "psych", "purrr", "readr", "rlang", "scales", "stringr", "tibble", "tidybayes", "tidyr", "tidyselect"
))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

Figures use the Open Sans font family, which needs to be installed for them to render as in the paper.

### Order

Run each script from the repository root, which it treats as the working directory:

1. `qnr_analysis_ancova.R`. Run this first: the mediation summaries in the two modelling scripts check their results against its output (`outputs/ancova_primary.csv`).
2. `modelling_reward-effort.R`, then `modelling_causal-attr.R`, **in the same R session**. Figure 5 combines panels from both scripts.
3. `task_covariance_collinearity.R`, which reads the model fits saved by step 2.

No fitted models are included: every model is fitted from scratch by MCMC, so a full run takes some time. Fits and draws are saved to `STAND-tasks-2/analysis/stan-fits/`, and tables to `outputs/`. Because the IDs are scrambled, participants are ordered differently from the original analysis. Group-level results will match the paper up to MCMC error, but individual-level quantities will be ordered differently.
