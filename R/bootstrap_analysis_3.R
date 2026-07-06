# Participant-level (cluster) bootstrap for the analysis-3 two-stage pipeline.
#
# Each replicate re-imputes the contact ages (a fresh call to
# process_fc_data_conmat()), resamples the participants with replacement giving
# each draw a fresh part_id, and re-runs the whole downstream pipeline
# (age-structured contact model -> residual heterogeneity sigma -> activity
# matrix -> R0s -> final sizes for the three scenarios). We store, per replicate,
# the resulting quantities plus the participant membership (which participants
# were drawn), so that a BCa acceleration term can be estimated later via
# jackknife-after-bootstrap without re-running anything.
#
# This is expensive (~15 s / replicate; ~40 min for B = 1000 across 7 workers),
# so it is kept out of the Quarto document: run this script once and the qmd
# reads the cached results from data/boot_analysis_3.rds.
#
# Usage:  Rscript R/bootstrap_analysis_3.R [B]     (B defaults to 1000)

suppressMessages({
  source("R/packages.R")
  source("R/functions.R")
})

# number of bootstrap replicates (override from the command line for testing)
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 1000L

# fixed configuration, shared with the point-estimate analysis in the qmd
beta <- 0.1
alpha <- 0.5
epsilon <- 0.5
n_activity_bins <- 10
max_age <- 90

# population / demography is external census data, held fixed across bootstraps
fc_population <- age_population(
  data = socialmixr::wpp_age(),
  location_col = country,
  location = c("France"),
  age_col = lower.age.limit,
  year_col = year,
  year = 2010
)

# canonical participant ordering, so membership indices are comparable across
# replicates. The set of participants is fixed (imputation only changes contact
# ages), so this can be established once.
proc0 <- process_fc_data_conmat()
canonical_ids <- as.character(unique(proc0$part_id))
rm(proc0)

# helper to load, re-impute and rename the participant-level data
load_fc_participant_data <- function() {
  process_fc_data_conmat() |>
    rename(age_from = part_age, age_to = cont_age) |>
    mutate(setting = "all", part_id = as.character(part_id))
}

# one bootstrap replicate: re-impute, resample participants, run the pipeline
one_boot <- function(i) {
  contact_data <- load_fc_participant_data()
  resample <- resample_participants(contact_data, canonical_ids)
  result <- tryCatch(
    fc_final_size_pipeline(
      resample$data,
      population = fc_population,
      beta = beta,
      alpha = alpha,
      epsilon = epsilon,
      n_activity_bins = n_activity_bins,
      max_age = max_age
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  list(rep = i, membership = resample$membership, result = result)
}

# also compute the point estimate on the full (unresampled) data, so the doc's
# central line and the BCa bias-correction use exactly the same estimator as the
# replicates
point_data <- load_fc_participant_data()
point_estimate <- fc_final_size_pipeline(
  point_data,
  population = fc_population,
  beta = beta,
  alpha = alpha,
  epsilon = epsilon,
  n_activity_bins = n_activity_bins,
  max_age = max_age
)
rm(point_data)

# run the replicates in parallel, with parallel-safe RNG streams
n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)

# date-mnemonic seed (see repo convention); with furrr_options(seed = TRUE) this
# seeds reproducible independent L'Ecuyer streams for the workers
set.seed(2026-07-06)

message(sprintf("running %d bootstrap replicates on %d workers...", B, n_workers))
boot_start <- Sys.time()

replicates <- future_map(
  seq_len(B),
  one_boot,
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)

plan(sequential)
boot_secs <- as.numeric(difftime(Sys.time(), boot_start, units = "secs"))

# report how many replicates succeeded
failed <- vapply(replicates,
                 function(x) !is.null(x$result$error),
                 logical(1))
message(sprintf("done in %.1f min; %d / %d replicates succeeded (%d failed)",
                boot_secs / 60, sum(!failed), B, sum(failed)))

boot_results <- list(
  point_estimate = point_estimate,
  replicates = replicates,
  canonical_ids = canonical_ids,
  config = list(beta = beta, alpha = alpha, epsilon = epsilon,
                n_activity_bins = n_activity_bins, max_age = max_age, B = B),
  n_failed = sum(failed)
)

if (!dir.exists("data")) dir.create("data")
saveRDS(boot_results, "data/boot_analysis_3.rds")
message("saved data/boot_analysis_3.rds")
