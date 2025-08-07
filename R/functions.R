# Download all French Connection study data:

# Béraud G, Kazmercziak S, Beutels P, Levy-Bruhl D, Lenne X, Mielcarek N, et al.
# (2015) The French Connection: The First Large Population-Based Contact Survey
# in France Relevant for the Spread of Infectious Diseases. PLoS ONE 10(7):
# e0133203. https://doi.org/10.1371/journal.pone.0133203

# data on zenodo here: https://zenodo.org/records/3886590

fc_data_dir <- function() {
  "data/french_connection"
}

fc_data_filepath <- function(who = c("participant", "contact"),
                             details = c("common", "extra")) {
  who <- match.arg(who)
  details <- match.arg(details)
  
  filename <- paste0(
    "2015_Beraud_France_",
    who,
    "_",
    details,
    ".csv"
  )
  file.path(
    fc_data_dir(),
    filename
  )
}

fc_data_filepaths <- function() {
  c(
    fc_data_filepath("participant", "common"),
    fc_data_filepath("participant", "extra"),
    fc_data_filepath("contact", "common"),
    fc_data_filepath("contact", "extra")
  )
}

# create a directory, download the file, and unzip it
download_fc_data <- function() {
  dir.create(fc_data_dir(),
             recursive = TRUE,
             showWarnings = FALSE)
  fc_filepath <- file.path(fc_data_dir(),
                           "everything.zip")
  download.file("https://zenodo.org/api/records/3886590/files-archive",
                fc_filepath)
  unzip(fc_filepath,
        exdir = fc_data_dir())
}

# check whether the French Connection data has been downloaded, and if not,
# download it
maybe_download_fc_data <- function() {
  files_exist <- vapply(fc_data_filepaths(),
                        file.exists,
                        FUN.VALUE = TRUE)
  if (!all(files_exist)) {
    message("Downloading French Connection data")
    download_fc_data()
  }
}

# Download all Hong Kong Survey data

# Kwok Kin On, Cowling Ben, Wei Vivian, Riley Steven and Read Jonathan M.
# 2018 Temporal variation of human encounters and the number of locations in
# which they occur: a longitudinal study of Hong Kong residentsJ. R. Soc.
# Interface. 15:20170838 https://doi.org/10.1098/rsif.2017.0838

# direct csv download here:
# https://royalsocietypublishing.org/action/downloadSupplement?doi=10.1098%2Frsif.2017.0838&file=rsif20170838supp2.csv

hk_data_dir <- function() {
  "data/hongkong"
}

hk_data_filepath <- function() {
  file.path(
    hk_data_dir(),
    "HongKongData.csv"
  )
}

# create a directory and download the CSV there
download_hk_data <- function() {
  dir.create(hk_data_dir(),
             recursive = TRUE,
             showWarnings = FALSE)
  download.file("https://royalsocietypublishing.org/action/downloadSupplement?doi=10.1098%2Frsif.2017.0838&file=rsif20170838supp2.csv",
                hk_data_filepath())
}

# check whether the Hong Kong data has been downloaded, and if not, download it
maybe_download_hk_data <- function() {
  files_exist <- vapply(hk_data_filepaths(),
                        file.exists,
                        FUN.VALUE = TRUE)
  if (!all(files_exist)) {
    message("Downloading Hong Kong data")
    download_hk_data()
  }
}

# load and process the French Connection data for the variance partitioning
# analysis
process_fc_data_varpart <- function() {
  
  # download the French Connection data, if needed
  maybe_download_fc_data()
  
  # participant info
  fc_participants <- read_csv(
    fc_data_filepath("participant",
                     "common"))
  fc_participants_extra <- read_csv(
    fc_data_filepath("participant",
                     "extra"))
  
  # get all observed combinations of participant and wave, to pad contacts with
  # 0s
  fc_participants_observed <- fc_participants_extra |>
    select(part_id,
           wave,
           participant_occupation) |>
    rename(
      part_occupation = participant_occupation 
    ) |>
    expand_grid(
      studyDay = 1:2
    ) |>
    left_join(
      fc_participants,
      by = "part_id"
    ) |>
    select(
      -hh_id
    )
  
  # contact events
  fc_contacts <- read_csv(
    fc_data_filepath("contact",
                     "common"))
  fc_contacts_extra <- read_csv(
    fc_data_filepath("contact",
                     "extra"))
  
  # collapse contact events to count contacts per participant, per wave, per
  # studyDay
  fc_contacts_sry <- fc_contacts |>
    # add on wave and day columns
    left_join(fc_contacts_extra,
              by = "cont_id") |>
    # count contacts per participant/wave/day
    group_by(
      part_id,
      wave,
      studyDay
    ) |>
    summarise(
      contacts = n(),
      .groups = "drop"
    )
  
  # add observed contacts, fill in others with 0s
  fc_participants_observed |>
    left_join(fc_contacts_sry,
              by = c("part_id", "wave", "studyDay")) |>
    mutate(
      contacts = replace_na(contacts, 0)
    ) |>
    # set factors
    mutate(
      part_id = factor(part_id),
      part_age = factor(part_age),
      part_gender = factor(part_gender),
      part_occupation = factor(part_occupation),
      wave = factor(wave),
      studyDay = factor(studyDay)
    ) |>
    relocate(
      part_age,
      part_gender,
      part_occupation,
      .after = part_id
    )
}

# load and process the Hong Kong data (pre-prepared by Leon) for the variance
# partitioning analysis
process_hk_data_varpart <- function() {
  
  # download the Hong Kong data if needed
  maybe_download_hk_data()
  
  hk_data_filepath() |>
    read_csv() |>
    select(
    part_id = pid,
    part_age = age,
    part_gender = sex,
    contacts = n.contact.total
  ) |> 
    mutate(
      part_id = as_factor(part_id),
      part_age = as_factor(part_age),
      part_gender = as_factor(part_gender)
    )
  
}

# a hack because this errors with the native pipe:
#  x |> `[[`(y)
# so we can do:
#  x |> extracting(y)
extracting <- function(x, y) {
  x[[y]]
}

# given a fitted random effects model, perform the variance partitioning
partition_variance_lmer <- function (model) {
  model |>
    summary() |>
    extracting("varcor") |>
    as_tibble() |>
    mutate(
      var = sdcor ^ 2,
      proportion = var / sum(var)
    ) |>
    select(
      grp,
      var,
      proportion
    ) |>
    mutate(
      grp = case_when(
        grp == "part_age" ~ "between ages",
        grp == "part_gender" ~ "between genders",
        grp == "part_occupation" ~ "between occupations",
        grp == "Residual" ~ "within individuals",
        .default = "between individuals (unexplained)",
      )
    ) |>
    rename(
      partition = grp
    )
  
}

make_barplot <- function(partitioning) {
  
  # set up colour palette for bar charts
  colour_palette <- c(
    "within individuals" = "pink",
    "between individuals (unexplained)" = scales::alpha("blue", 0.1),
    "between genders" = scales::alpha("blue", 0.4),
    "between ages" = scales::alpha("blue", 0.5),
    "between occupations" =  scales::alpha("blue", 0.6)
  )
  
  partitioning |>
    # set the factor order for partition to how we want to plot it
    mutate(
      partition = factor(partition,
                         levels = names(colour_palette))
    ) |>
    # sort by partition (to the factor order), so we can place the labels in the
    # right places
    arrange(
      partition
    ) |>
    group_by(
      Study
    ) |>
    # for labels
    mutate(
      text_percent = scales::label_percent()(proportion),
      text_position = rev(cumsum(rev(proportion))) - proportion / 2
    ) |>
    ungroup() |>
    # make some better names
    rename(
      `Variance explained` = proportion,
      Component = partition
    ) |>
    ggplot(
      aes(
        x = Study,
        y = `Variance explained`,
        fill = Component
      )
    ) + 
    geom_bar(
      stat = "identity"
    ) +
    geom_text(
      aes(
        label = text_percent,
        y = text_position
      )
    ) +
    scale_y_continuous(
      labels = scales::label_percent()
    ) +
    scale_fill_manual(values = colour_palette) +
    theme_minimal()
}

# load and process the French Connection data for learning the age-structured
# contct matrix, augmented with additional heterogeneity classes
process_fc_data_conmat <- function() {
  
  # download the French Connection data, if needed
  maybe_download_fc_data()
  
  # participant info
  fc_participants <- read_csv(
    fc_data_filepath("participant",
                     "common"))
  fc_participants_extra <- read_csv(
    fc_data_filepath("participant",
                     "extra"))
  
  # get all observed combinations of participant and wave, to pad contacts with
  # 0s
  fc_participants_observed <- fc_participants_extra |>
    select(part_id,
           wave) |>
    expand_grid(
      studyDay = 1:2
    ) |>
    left_join(
      fc_participants,
      by = "part_id"
    ) |>
    select(
      -hh_id,
      -part_gender
    )
  
  # contact events
  fc_contacts <- read_csv(
    fc_data_filepath("contact",
                     "common"))
  fc_contacts_extra <- read_csv(
    fc_data_filepath("contact",
                     "extra"))  
  
  # collapse contact events to count contacts per participant, per wave, per
  # studyDay, per contact age
  fc_contacts_sry <- fc_contacts |>
    # deal with missing and interval contact ages but imputing from a uniform
    # distribution (in a vectorised way, for marginal computational efficiency
    # gain)
    mutate(
      # if missing the age range (127 of ~39K), take the outer limits of ages
      # from the survey
      cnt_age_est_min = replace_na(cnt_age_est_min,
                                   min(cnt_age_est_min, na.rm = TRUE)),
      cnt_age_est_max = replace_na(cnt_age_est_max,
                                   max(cnt_age_est_max, na.rm = TRUE)),
      # impute the contact age for all of these
      cont_age_range = cnt_age_est_max - cnt_age_est_min,
      cont_age_imputed = cnt_age_est_min + cont_age_range * runif(n()),
      # use the exact age preferentially, or else the imputed age
      cont_age = case_when(
        !is.na(cnt_age_exact) ~ cnt_age_exact,
        .default = round(cont_age_imputed)
      ),
      .after = everything()
    ) |>
    # add on wave and day columns
    left_join(fc_contacts_extra,
              by = "cont_id") |>
    # count contacts per participant/wave/day
    group_by(
      part_id,
      wave,
      studyDay,
      cont_age
    ) |>
    summarise(
      contacts = n(),
      .groups = "drop"
    )
  
  # get the range of contact ages, and pad the observations with explicit 0s
  possible_ages <- seq(
    from = min(fc_contacts_sry$cont_age),
    to = max(fc_contacts_sry$cont_age),
    by = 1)
  
  fc_contacts_complete <- fc_contacts_sry |>  
    complete(
      # only the participant/wave/studyDay combinations in the data
      nesting(part_id, wave, studyDay),
      # but for all possible ages
      cont_age = possible_ages,
      # if the value is not int he contacts file, there must be 0 contacts
      fill = list(contacts = 0)
    )
    
  # add observed contacts, fill in others with 0s
  fc_participants_observed |>
    left_join(fc_contacts_complete,
              by = c("part_id", "wave", "studyDay")) |>
    mutate(
      contacts = replace_na(contacts, 0)
    ) |>
    # set factors
    mutate(
      part_id = factor(part_id),
      wave = factor(wave),
      studyDay = factor(studyDay)
    ) |>
    relocate(
      part_age,
      .before = cont_age
    )
}

# supporting function for quantizing
measure_mean_by_contact_quantile <- function(contacts, n_quantiles) {
  contacts |> 
    mutate(
      means_by_quantile = ntile(X, n_quantiles)
    ) |> 
    group_by(
      means_by_quantile
    ) |> 
    summarise(
      means_by_quantile = mean(X)
    ) |> 
    mutate(
      quantile = row_number()
    )
}

generate_proportionate_mixing_matrix <- function(activity_classes) {

  # average contacts, accounting for uneven population fractions
  mean_degree <- sum(activity_classes$activity * activity_classes$fraction)
  
  cross_join(
    activity_classes,
    activity_classes
  ) |> 
    mutate(
      # for each contact, the probability the contactee is in the y compartment
      prob_contactee = activity.y * fraction.y / mean_degree,
      weight = activity.x * prob_contactee,
    ) |> 
    select(
      from = class.x,
      to = class.y,
      weight
    )
}

# Full assortativity - if the class is the same, they get that activity level.
# No effect of class population fractions on mixing
generate_fully_assortative_mixing_matrix <- function(activity_classes) {
  
  cross_join(
    activity_classes,
    activity_classes
  ) |>
    rename(
      from = class.x,
      to = class.y
    ) |> 
    mutate(
      weight = ifelse(from == to,
                      activity.x,
                      0)
    ) |> 
    select(
      from,
      to,
      weight
    )
}

# use numerical integration to estimate the mean activity level in a bin bounded
# by lower and upper activity classes, for a given value of sigma
bin_mean_activity <- function(lower, upper, sigma) {
  f <- function(x) {
    x * dlnorm(x, meanlog = 0, sdlog = sigma)
  }
  # compute the expectation
  int <- integrate(f,
                   lower = lower,
                   upper = upper)
  
  # correct for truncation (normalises the densities to the truncated range)
  cdf_upper <- plnorm(upper, meanlog = 0, sdlog = sigma)
  cdf_lower <- plnorm(lower, meanlog = 0, sdlog = sigma)
  interval <- cdf_upper - cdf_lower
  
  int$value / interval
}

# quantile method, but with finite support and numerical integration to get mean
# activity level per bin
quantile2_classes_mean <- function(sigma,
                                   n_classes,
                                   alpha = 1e-5) {
  
  # upper limit for integration
  max_quantile <- 1 - alpha
  # max_activity <- qlnorm(max_quantile, meanlog = 0, sdlog = sigma)
  
  # define quantiles and get lower and upper bounds of the bins
  probs <- seq(0, max_quantile, length.out = n_classes + 1)
  breaks <- qlnorm(probs, meanlog = 0, sdlog = sigma)
  lower_bounds <- breaks[-length(breaks)]
  upper_bounds <- breaks[-1]
  
  # compute the mean activity level in each bin
  activity <- mapply(bin_mean_activity,
                     lower_bounds,
                     upper_bounds,
                     MoreArgs = list(sigma = sigma))
  
  # get the population fraction in each bin, accounting for upper truncation
  cdf_upper <- plnorm(upper_bounds,
                      meanlog = 0,
                      sdlog = sigma)
  cdf_lower <- plnorm(lower_bounds,
                      meanlog = 0,
                      sdlog = sigma)
  
  fraction <- (cdf_upper - cdf_lower) / max_quantile
  
  # return the tibble
  tibble(
    class = seq_len(n_classes),
    activity = activity,
    fraction = fraction
  )
  
}


# given an activity level standard deviation sigma and target number of classes,
# use the Monte Carlo quantile method to the mean activity level in each class.
# Control the precision of the MC approximation with n_per_class (total sims =
# n_classes * n_per_class)
quantile_classes_mean <- function(sigma,
                                  n_classes,
                                  n_per_class = 10000) {
  
  n_sims <- n_classes * n_per_class
  normal_dist_test.df <- as_tibble(
    x = exp(rnorm(n = n_sims,
                  mean = 0,
                  sd = sigma))
  ) |> 
    rename(
      X = value
    )
  
  s <- measure_mean_by_contact_quantile(normal_dist_test.df,
                                        n_classes)
  
  s |>
    mutate(
      fraction = 1 / n_classes
    ) %>%
    select(
      class = quantile,
      activity = means_by_quantile,
      fraction
    )
  
}

# given an activity level standard deviation sigma and target number of classes,
# use Gaussian Legendre quadrature to define activity classes and compute the
# population fraction (integration weight) and average activity level
# (integration point location) in each class. alpha controls the bounds of ther
# integration (from alpha/2 to 1-alpha/2), and the integral is corrected for the
# missing tails
gaussquad_classes <- function(sigma,
                              n_classes,
                              alpha = 0.01) {
  quantile_bounds <- c(alpha / 2, 1 - alpha / 2)
  log_activity_bounds <- qnorm(quantile_bounds,
                               mean = 0,
                               sd = sigma)
  log_activity_range <- log_activity_bounds[2] - log_activity_bounds[1]
  
  # compute the integration points
  if (n_classes > 1) {
    quads <- gaussquad::hermite.he.quadrature.rules(n_classes)[[n_classes]]
  } else {
    # handle single-class situation
    quads <- data.frame(
      x = 0,
      w = 2.506628
    )
  }
  
  quads |>
    as_tibble() |>
    # normalise the weights to population fractions
    mutate(
      fraction = w / sum(w)
      # fraction = w
    ) |>
    mutate(
      # convert form standard normal to lognormal with sigma
      activity = exp(rev(x * sigma))
    ) |>
    # add classes
    mutate(
      class = row_number()
    ) |>
    select(
      class,
      activity,
      fraction
    )

  # missing tail density needs to be added on, by dividing activities by 1 -
  # alpha
  
}

# given sigma, a number of classes, and a method, return a binning scheme with
# activity classes
get_activity_classes <- function(sigma,
                                 n_classes,
                                 method = c("gaussquad",
                                            "quantile_mean",
                                            "quantile2_mean")) {
  switch(
    method,
    quantile_mean = quantile_classes_mean(sigma,
                                          n_classes),
    quantile2_mean = quantile2_classes_mean(sigma,
                                            n_classes),
    gaussquad = gaussquad_classes(sigma,
                                  n_classes)
  )
}

# estimate the mean of a lognormal distribution with parameter `sigma`,
# using `n_classes` bins defined by method `method`.
bin_estimate_mean <- function(sigma,
                              n_classes,
                              method = method) {
  
  activity_classes <- get_activity_classes(sigma = sigma,
                                           n_classes = n_classes,
                                           method = method)
  sum(activity_classes$activity * activity_classes$fraction)
  
}
# estimate the sd of a lognormal distribution with parameter `sigma`,
# using `n_classes` bins defined by method `method`.
bin_estimate_sd <- function(sigma,
                            n_classes,
                            method = method) {
  
  activity_classes <- get_activity_classes(sigma = sigma,
                                           n_classes = n_classes,
                                           method = method)
  
  mean <- sum(activity_classes$activity * activity_classes$fraction)
  
  sqrt(sum(activity_classes$activity ^ 2 * activity_classes$fraction) - mean ^ 2)
  
}

# wrapper function
# sigma is the normal distribution standard deviation
# assort is the amount of assortativity
# n_classes is the number of classes to use to approximate the activity
#   level variance
# method is the method to use for integration: 'quantile_mean' to use quantiles
#   of equal size and quantile mean activity levels estimated by Monte Carlo
#   simulation, 'quantile2_mean' does the same by numerical integration,
#   'gaussquad' uses Gaussian Legendre quadrature points, with different weights
#   (population fractions of classes) and fixed quadrature points.
generate_matrix <- function(sigma = 1,
                            assort = 1,
                            n_classes = 2,
                            method = c("gaussquad",
                                       "quantile_mean",
                                       "quantile2_mean")) {
  
  activity_classes <- get_activity_classes(sigma = sigma,
                                           n_classes = n_classes,
                                           method = method)
  
  assortative_matrix <- generate_fully_assortative_mixing_matrix(
    activity_classes)
  proportionate_matrix <- generate_proportionate_mixing_matrix(
    activity_classes)
  
  assortative_matrix |> 
    rename(
      weight_assortative = weight
    ) |> 
    mutate(
      weight_proportionate = proportionate_matrix$weight
    ) |> 
    mutate(
      weight_combined = assort * weight_assortative +
        (1 - assort) * weight_proportionate
    ) |> 
    select(
      from,
      to,
      weight_combined
    )
  
}

# M is the matrix in long format, eigen_value is which eigenvalue to extract 
# defaults to eigen_value = 1 for the dominant 
matrix_to_eigenvalue <- function(M, eigen_value = 1) {
  decomp <- M |>
    usedist::pivot_to_numeric_matrix(
      from,
      to,
      weight_combined
    ) |> 
    eigen() 
  Re(decomp$values[eigen_value])
}

map_to_eigen <- function(sigma = 1,
                         assort = 1,
                         n_classes = 3,
                         eigen_value = 1,
                         beta = 1,
                         method = c("gaussquad",
                                    "quantile_mean",
                                    "quantile2_mean")) {
  method <- match.arg(method)
  M <- generate_matrix(sigma = sigma,
                       assort = assort,
                       n_classes = n_classes,
                       method = method)
  M |>
    mutate(
      weight_combined = weight_combined * beta
    )
  
  # NOTE: the line above has no effect - presumably unintended, but currently
  # not used
  
  matrix_to_eigenvalue(M, eigen_value = eigen_value)
}

# make a proportionate mixing contact matrix from class contact rates
make_contact_matrix <- function(class_means) {
  class_means %*% t(class_means) / sum(class_means)
}

# get the dominant eigenvalue of a matrix
get_eigenval <- function(matrix) {
  Re(eigen(matrix)$values[1])
}

# get it faster (but mor approximately) using the power method
fast_eigenval <- function(matrix, tol = 0.001, maxiter = 1000) {
  eig <- fastmatrix::power.method(matrix,
                                  maxiter = maxiter,
                                  tol = tol)  
  eig$value
}
