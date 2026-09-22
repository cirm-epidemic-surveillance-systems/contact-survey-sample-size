# Handoff: bootstrap-vs-sample-size analysis (2026-09-17 → resume after restart)

Written because the machine is about to restart mid-task. Everything below survives
the restart (it's all on disk); only the live Claude Code conversation context is lost.
Read this file, then read the **Provenance** section at the bottom of
`R/analysis_4_bootstrap_sample_size.qmd` for the full prompt-by-prompt history — this
file is just a shorter orientation on top of that.

## Status: done, rendered, verified

The task (bootstrap uncertainty vs. sample size, using this repo's own paper's
activity-heterogeneity formalism rather than the simpler bristol-vaccine-centre
reference-repo approach) is **complete and has been successfully rendered end-to-end**
at full scale (`numreps = 200`, both datasets, all 3 assortativity values). Nothing is
mid-run. If you're resuming, the most likely next steps are refinements/extensions
(see "Possible next steps" below), not finishing something broken.

- New file (untracked, not committed): `R/analysis_4_bootstrap_sample_size.qmd`
- Rendered output: `R/analysis_4_bootstrap_sample_size.html` (gitignored, regenerate with
  the render command below)
- Figures: `figures/bootstrap_*.png` (gitignored, 9 files: per-dataset σ_ind fig2 +
  σ_ind uncertainty + eigenvalue fig2 + structural-vs-statistical, plus one
  cross-dataset comparison)
- No existing repo files were modified. `data/french_connection/` and R packages (below)
  were added as side effects.

## What the analysis does (one paragraph)

For Hong Kong and French Connection (the two local datasets with genuine repeated
contact measurements), at each of ~30 sample sizes, run 200 cluster-bootstrap
replicates: resample participants with replacement (relabelling repeats as distinct
virtual participants, keeping each one's full repeat-day/wave record together),
rebuild the age-structured contact matrix, **re-estimate σ_ind** (activity
heterogeneity) via `lmer(log1p(total_contacts) ~ (1|age_band) + (1|virtual_participant))`,
build activity classes via Gaussian quadrature, build the activity matrix via **Tom's**
iterative-normalisation kernel method at 3 fixed α values (5, 10, 20 — the paper's own
sensitivity grid), Kronecker-combine with the age matrix, and take the dominant
eigenvalue. This separates two kinds of uncertainty: **statistical** (in σ_ind and in
the eigenvalue at fixed α — shrinks with sample size) vs. **structural** (the gap
between α values — does not shrink with sample size, since α isn't identifiable from
the data at all). Confirmed in the actual results: all four model series (age-only,
α=5/10/20) stay flat and separated in mean-eigenvalue-vs-n, while all four SD-vs-n
curves shrink together.

## Environment setup already done (don't redo)

R packages installed into the standard user library (`~/Library/R/arm64/4.6/library`,
auto-discovered by plain `Rscript`/`quarto`, no `.libPaths()` hacks needed):
`tidyverse`, `furrr`, `patchwork`, `lme4`, `gaussquad` (+ their dependencies). These
persist across a reboot. Note: `R/packages.R` also wants `conmat`, `usedist`, `cowplot`,
`RCurl` — these were **not** installed, and the new qmd deliberately avoids
`source("R/packages.R")` for exactly this reason (it would fail on `library(conmat)`,
a GitHub-only package). It only does `source("R/functions.R")` (safe — just function
definitions) plus its own explicit `library()` calls.

Data: `data/french_connection/` was downloaded from Zenodo (record 3886590) into the
layout `R/functions.R`'s existing `download_fc_data()` expects, and is gitignored.
`data/hongkong/HongKongData.csv` was already present before this task started.

Quarto CLI is not on `PATH` in this environment — use RStudio's bundled copy:
```
/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto render R/analysis_4_bootstrap_sample_size.qmd
```
Run from the repo root. Full render (both datasets, 200 reps, 3 α values) takes under
10 minutes. Per-chunk caching (`#| cache: true`) means a second render only re-runs
chunks whose code actually changed.

## Possible next steps (not started, just ideas)

- **Fill the paper's own TODO**: `Paper/main.tex:224` has a literal placeholder —
  *"We estimate the value of the activity heterogeneity parameter to be σ=?? in the
  French dataset and σ=?? in the Hong Kong dataset."* — our bootstrap's σ_ind mean at
  the largest sample size for each dataset is a real candidate answer to that. Worth
  pointing out to Leon; not yet extracted/reported anywhere as a headline number.
- Could add a without-replacement comparison (declined earlier — "with-replacement
  only" was the explicit choice).
- Could push `numreps` higher than 200 if smoother curves are wanted; timing budget is
  well characterised (~0.07s/replicate Hong Kong, ~0.083s/replicate French Connection,
  linear scaling, zero hard `lmer` failures at any sample size tested).
- `README.md` is still just a placeholder ("Team 5: Contact surveys — Watch this space")
  and was deliberately left untouched; might eventually want a pointer to this analysis.
- Nothing has been committed to git. `R/analysis_4_bootstrap_sample_size.qmd` and its
  `_cache/` dir are untracked (cache dir is arguably not something you want to commit —
  Leon hasn't been asked about this yet).

## Delete this file

This file (`HANDOFF.md`) is scratch/orientation, not a deliverable — safe to delete once
you've picked the thread back up and don't need it anymore.
