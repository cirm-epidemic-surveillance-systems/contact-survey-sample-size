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
  files_exist <- vapply(hk_data_filepath(),
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
  
  # get the French Connection datasets
  fc_datasets <- load_fc_datasets()
  
  # get all observed combinations of participant and wave, to pad contacts with
  # 0s
  fc_participants_observed <- fc_datasets$fc_participants_extra |>
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
      fc_datasets$fc_participants,
      by = "part_id"
    ) |>
    select(
      -hh_id
    )
  
  # collapse contact events to count contacts per participant, per wave, per
  # studyDay
  fc_contacts_sry <- fc_datasets$fc_contacts |>
    # add on wave and day columns
    left_join(fc_datasets$fc_contacts_extra,
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
  
  hk_data <- hk_data_filepath() |>
    read_csv(
      col_types = cols(
        pid = col_double(),
        hid = col_double(),
        age = col_double(),
        sex = col_character(),
        n.samples = col_double(),
        wave = col_character(),
        reporting.day = col_character(),
        n.contact.total = col_double(),
        n.loc.group = col_character(),
        `n.contact.0-4` = col_double(),
        `n.contact.5-19` = col_double(),
        `n.contact.20-39` = col_double(),
        `n.contact.40-64` = col_double(),
        `n.contact.65+` = col_double()
      )
    )
  
  hk_data |>
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

# given a fitted Poisson-log GLMM that includes an observation-level random
# effect named "obs", perform the variance partitioning on the latent (log)
# scale. Unlike a Gaussian LMM, a Poisson-log model has no residual variance
# parameter: the within-individual, day-to-day variation is split into (i) an
# observation-level random effect capturing extra-Poisson (overdispersion)
# variation, and (ii) a "distribution-specific" term for the Poisson sampling
# variance, expressed on the log scale. Both are pooled into the "within
# individuals" partition so the result is directly comparable to
# partition_variance_lmer(). See Nakagawa, Johnson & Schielzeth (2017).
partition_variance_glmer <- function (model) {

  # random-effect variances, on the log (latent) scale. glmer does not report a
  # "Residual" row for a Poisson family (the scale is fixed), but filter it out
  # defensively in case it appears.
  var_components <- model |>
    summary() |>
    extracting("varcor") |>
    as_tibble() |>
    filter(grp != "Residual") |>
    select(
      grp,
      var = vcov
    )

  # population-mean count: lambda = exp(beta0 + 0.5 * sum of RE variances)
  beta0 <- lme4::fixef(model)[["(Intercept)"]]
  lambda <- exp(beta0 + 0.5 * sum(var_components$var))

  # distribution-specific (Poisson sampling) variance on the log scale. The
  # trigamma function is the most exact approximation (Nakagawa et al. 2017); the
  # lognormal approximation log(1 + 1 / lambda) agrees to within ~1% at the
  # moderate counts seen here.
  var_dist <- psigamma(lambda, deriv = 1)
  # var_dist <- log(1 + 1 / lambda)  # lognormal approximation (alternative)

  var_components |>
    # add the Poisson sampling variance as its own (within-individual) component
    bind_rows(
      tibble(grp = "poisson_sampling", var = var_dist)
    ) |>
    mutate(
      partition = case_when(
        grp == "part_age" ~ "between ages",
        grp == "part_gender" ~ "between genders",
        grp == "part_occupation" ~ "between occupations",
        grp %in% c("obs", "poisson_sampling") ~ "within individuals",
        .default = "between individuals (unexplained)",
      )
    ) |>
    # pool the observation-level RE and Poisson sampling into "within
    # individuals"
    group_by(
      partition
    ) |>
    summarise(
      var = sum(var),
      .groups = "drop"
    ) |>
    mutate(
      proportion = var / sum(var)
    ) |>
    select(
      partition,
      var,
      proportion
    )

}

# format a proportion (0-1) as a percentage label: to the nearest integer, but
# to a single decimal place when it would otherwise round to 0%, and "<0.1%"
# when even one decimal place would round to 0
format_percent_label <- function(proportion) {
  pct <- 100 * proportion
  case_when(
    round(pct) >= 1 ~ sprintf("%.0f%%", pct),
    round(pct, 1) >= 0.1 ~ sprintf("%.1f%%", pct),
    .default = "<0.1%"
  )
}

# as format_percent_label() but always to one decimal place (used in the tables,
# where there is room for the extra precision); still "<0.1%" below 0.05%
format_percent_label_1dp <- function(proportion) {
  pct <- 100 * proportion
  ifelse(round(pct, 1) >= 0.1, sprintf("%.1f%%", pct), "<0.1%")
}

make_barplot <- function(partitioning, drop = TRUE, label_threshold = 0.05) {

  # semantic colour palette:
  # - within individuals: light grey - the dominant component, but the one of
  #   least significance in the downstream analyses, so de-emphasised
  # - between individuals (unexplained): green, matching the activity-structured
  #   contact matrix panel in analysis 3
  # - between ages: blue, matching the age-structured contact matrix panel in
  #   analysis 3
  # - between genders / occupations: distinguishable shades of the age blue, as
  #   they are largely confounded with age
  colour_palette <- c(
    "within individuals" = "grey88",
    "between individuals (unexplained)" = "#41AE76",
    "between ages" = "#2171B5",
    "between genders" = "#9ECAE1",
    "between occupations" = "#08306B"
  )

  plot_data <- partitioning |>
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
      text_percent = format_percent_label(proportion),
      text_position = rev(cumsum(rev(proportion))) - proportion / 2
    ) |>
    ungroup() |>
    # make some better names
    rename(
      `Variance explained` = proportion,
      Component = partition
    )

  # split the labels: large slices are labelled inside the bar; small slices
  # that would overlap are labelled outside with a leader line (coloured by
  # component, so which label belongs to which slice is unambiguous)
  labels_inside <- plot_data |>
    filter(`Variance explained` >= label_threshold)
  labels_outside <- plot_data |>
    filter(`Variance explained` < label_threshold)
  # nudge outside labels away from the bars: left for the first study, right for
  # the rest
  n_studies <- length(unique(plot_data$Study))
  # sit the labels tight against the outer edge of each bar: nudge just far
  # enough to clear the bar
  outside_nudge <- ifelse(
    as.integer(factor(labels_outside$Study)) <= n_studies / 2,
    -0.5,
    0.5
  )
  # align labels to the bar edge: the first (left-hand) study's labels are
  # right-aligned (right edge against the bar, text into the margin), the rest
  # left-aligned (left edge against the bar)
  outside_hjust <- ifelse(outside_nudge < 0, 1, 0)

  p <- ggplot(
      plot_data,
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
      data = labels_inside,
      aes(
        label = text_percent,
        y = text_position
      ),
      size = 4
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1)
    ) +
    # drop = FALSE keeps all components in the scale (and hence the legend) even
    # when a panel lacks some, so legends can be combined across panels
    scale_fill_manual(values = colour_palette, drop = drop) +
    # matched colour scale for the leader-line labels (no separate legend)
    scale_colour_manual(values = colour_palette, guide = "none") +
    xlab(NULL) +
    theme_minimal()

  # only add the leader-line layer (and the extra horizontal room it needs) when
  # there are small slices to label this way
  if (nrow(labels_outside) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = labels_outside,
        aes(
          label = text_percent,
          y = text_position,
          colour = Component
        ),
        nudge_x = outside_nudge,
        direction = "y",
        hjust = outside_hjust,
        size = 4,
        min.segment.length = 0,
        segment.colour = "black",
        segment.size = 0.25,
        box.padding = 0.1,
        show.legend = FALSE,
        seed = 2026
      ) +
      # just enough horizontal room for the label text to sit beside the bars
      scale_x_discrete(expand = expansion(add = 1))
  }

  p
}

# load the various French Connection datasets
load_fc_datasets <- function() {
  
  # download the French Connection data, if needed
  maybe_download_fc_data()
  
  # participant info
  fc_participants <- read_csv(
    fc_data_filepath("participant", "common"),
    col_types = cols(
      part_id = col_double(),
      hh_id = col_character(),
      part_age = col_double(),
      part_gender = col_character()
    )
  )
  
  # extra info on participants
  fc_participants_extra <- read_csv(
    fc_data_filepath("participant", "extra"),
    col_types = cols(
      part_id = col_double(),
      wave_part_id = col_double(),
      wave = col_double(),
      questionnaire.type = col_double(),
      childRespondentLink = col_double(),
      childRespondentAge = col_double(),
      childRespondentGender = col_double(),
      transportModeWeek1 = col_double(),
      transportModeWeek2 = col_double(),
      transportModeWeek3 = col_double(),
      transportModeWeekEnd1 = col_double(),
      transportModeWeekEnd2 = col_double(),
      transportModeWeekEnd3 = col_double(),
      ZIP = col_double(),
      participant_education = col_double(),
      participant_occupation = col_double(),
      participant_occ_detail = col_double(),
      enfScolarise = col_double(),
      enfGardeMaison = col_double(),
      enfGardeNounou = col_double(),
      enfNbEnfNounou = col_double(),
      enfNounouScolarise = col_double(),
      enfCreche = col_double(),
      enfNbEnfCreche = col_double(),
      enfFreqGarderie = col_double(),
      class_size = col_double(),
      enfCantine = col_double(),
      enfCentreAere = col_double(),
      enfCentreAereEcole = col_double(),
      enfCentreAereVacance = col_double(),
      work_contacts = col_double(),
      more20ContactPro = col_double(),
      work_contacts_nr = col_double(),
      AgeContactPro1 = col_double(),
      AgeContactPro2 = col_double(),
      AgeContactPro3 = col_double(),
      AgeContactPro4 = col_double(),
      AgeContactPro5 = col_double(),
      NbStudentClassroom = col_double(),
      studentCantina = col_double(),
      commonParticip = col_double()
    )
  )
  
  # contact events
  fc_contacts <- read_csv(
    fc_data_filepath("contact", "common"),
    col_types = cols(
      part_id = col_double(),
      cont_id = col_character(),
      cnt_age_exact = col_double(),
      cnt_age_est_min = col_double(),
      cnt_age_est_max = col_double(),
      cnt_gender = col_character(),
      cnt_home = col_logical(),
      cnt_work = col_logical(),
      cnt_school = col_logical(),
      cnt_transport = col_logical(),
      cnt_leisure = col_logical(),
      cnt_otherplace = col_logical(),
      frequency_multi = col_double(),
      phys_contact = col_double(),
      duration_multi = col_double()
    )
  )
  
  # extr info on contact events
  fc_contacts_extra <- read_csv(
    fc_data_filepath("contact", "extra"),
    col_types = cols(
      cont_id = col_character(),
      wave = col_double(),
      studyDay = col_double()
    )
  )  
  
  # return these as a list
  list(
    fc_participants = fc_participants,
    fc_participants_extra = fc_participants_extra,
    fc_contacts = fc_contacts,
    fc_contacts_extra = fc_contacts_extra
  )
  
}

# load and process the French Connection data for learning the age-structured
# contct matrix, augmented with additional heterogeneity classes
process_fc_data_conmat <- function() {
  
  # get the French Connection datasets
  fc_datasets <- load_fc_datasets()
  
  # get all observed combinations of participant and wave, to pad contacts with
  # 0s
  fc_participants_observed <- fc_datasets$fc_participants_extra |>
    select(part_id,
           wave) |>
    expand_grid(
      studyDay = 1:2
    ) |>
    left_join(
      fc_datasets$fc_participants,
      by = "part_id"
    ) |>
    select(
      -hh_id,
      -part_gender
    )
  
  # collapse contact events to count contacts per participant, per wave, per
  # studyDay, per contact age
  fc_contacts_sry <- fc_datasets$fc_contacts |>
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
    left_join(
      fc_datasets$fc_contacts_extra,
      by = "cont_id"
    ) |>
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
    tidyr::complete(
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
                                   max_quantile = 1 - 1e-5) {
  
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

# ---- Continuous kernel-based contact matrix functions ----------------------
# These implement the continuous-quantile-space models from the MATLAB scripts
# in the Matlab/ directory. They coexist with the discrete-class generate_matrix()
# pipeline above; the two approaches are theoretically distinct.
#
# All four functions accept a p_pop vector of population proportions per bin,
# supporting non-uniform bin sizes. For uniform quantile bins pass
# p_pop = rep(1/n_bins, n_bins).

# Gaussian assortativity kernel g(z) = exp(-alpha * z^2).
# alpha controls mixing sharpness: alpha = 0 gives proportionate mixing, large
# alpha gives near-fully-assortative mixing.
# (translation of Matlab/calcKernel.m)
calc_kernel <- function(z, alpha) {
  exp(-alpha * z^2)
}

# Construct a proportionate mixing contact matrix.
# (translation of Matlab/makeContactMatrix_PM.m)
#
# v     : numeric vector of activity levels in each bin (length n_bins)
# p_pop : numeric vector of population proportions per bin (length n_bins);
#         will be normalised internally
make_proportionate_contact_matrix <- function(v, p_pop) {
  p_pop <- p_pop / sum(p_pop)
  Ev    <- sum(p_pop * v)
  outer(p_pop * v, v) / Ev^2
}

# Construct an assortative contact matrix using Mike's iterative column method.
# (translation of Matlab/makeContactMatrix_AM.m)
#
# v     : numeric vector of activity levels in each bin (length n_bins)
# p_pop : numeric vector of population proportions per bin (length n_bins)
# alpha : non-negative scalar kernel width parameter
make_assortative_contact_matrix <- function(v, p_pop, alpha) {
  n_bins <- length(v)
  p_pop  <- p_pop / sum(p_pop)

  # Quantile midpoints derived from p_pop (matches MATLAB: x = 0.5*(c(1:end-1)+c(2:end)))
  c_cum <- c(0, cumsum(p_pop))
  x     <- 0.5 * (c_cum[1:n_bins] + c_cum[2:(n_bins + 1)])

  X  <- matrix(x, nrow = n_bins, ncol = n_bins, byrow = TRUE)
  Y  <- matrix(x, nrow = n_bins, ncol = n_bins, byrow = FALSE)
  gk <- calc_kernel(Y - X, alpha)
  Ev <- sum(p_pop * v)

  # C[i,j] = p_pop[i] * v[i] * gk[i,j]  (row-wise scaling)
  C <- (p_pop * v) * gk

  # Denominator: column sums of lower triangle of C
  C_lower        <- C
  C_lower[upper.tri(C_lower)] <- 0
  den <- colSums(C_lower)

  M1 <- sweep(C, 2, den, "/")

  # First column: scale by v[1]/Ev
  M1[, 1] <- M1[, 1] * v[1] / Ev

  # Subsequent columns: scale by remaining contact budget, accounting for p_pop
  for (j_col in 2:n_bins) {
    M1[j_col:n_bins, j_col] <- M1[j_col:n_bins, j_col] *
      (v[j_col] / Ev - sum(p_pop[1:(j_col - 1)] * M1[j_col, 1:(j_col - 1)]) / p_pop[j_col])
  }

  # Symmetrise: upper triangle = transpose of strict lower triangle, scaled by p_pop[i]/p_pop[j]
  M_lower <- M1
  M_lower[upper.tri(M_lower)] <- 0

  M_strict_lower <- M1
  M_strict_lower[upper.tri(M_strict_lower, diag = TRUE)] <- 0

  M <- M_lower + t(M_strict_lower) * outer(p_pop, 1 / p_pop)

  # Detailed balance check: p_pop[j] * M[i,j] must be symmetric
  agg_cont <- sweep(M, 2, p_pop, "*")
  stopifnot(max(abs(agg_cont - t(agg_cont))) < 1e-12)
  M
}

# Construct an assortative contact matrix using Tom's iterative normalisation method.
# (translation of Matlab/makeContactMatrix_AM_Tom.m)
#
# v        : numeric vector of activity levels in each bin (length n_bins)
# p_pop    : numeric vector of population proportions per bin (length n_bins)
# alpha    : non-negative scalar kernel width parameter
# rel_fact : relaxation factor for fixed-point iteration (default 0.5)
# tol      : convergence tolerance on sup-norm of relative update (default 1e-10)
make_toms_contact_matrix <- function(v, p_pop, alpha, rel_fact = 0.5, tol = 1e-10) {
  n_bins <- length(v)
  p_pop  <- p_pop / sum(p_pop)

  # Quantile midpoints derived from p_pop
  c_cum <- c(0, cumsum(p_pop))
  x     <- 0.5 * (c_cum[1:n_bins] + c_cum[2:(n_bins + 1)])

  X  <- matrix(x, nrow = n_bins, ncol = n_bins, byrow = TRUE)
  Y  <- matrix(x, nrow = n_bins, ncol = n_bins, byrow = FALSE)
  gk <- calc_kernel(Y - X, alpha)
  Ev <- sum(p_pop * v)

  w         <- rep(1, n_bins)
  converged <- FALSE
  while (!converged) {
    w_saved <- w
    # gk %*% (p_pop * w) is n x 1; as.vector() drops to plain vector
    w <- (1 - rel_fact) * w + rel_fact * v / (Ev * as.vector(gk %*% (p_pop * w)))
    converged <- max(abs(w - w_saved)) / max(abs(w_saved)) < tol
  }

  # M[i,j] = p_pop[i] * w[i] * w[j] * gk[i,j]
  T_mat <- outer(p_pop * w, w) * gk

  # Detailed balance check
  agg_cont <- sweep(T_mat, 2, p_pop, "*")
  stopifnot(max(abs(agg_cont - t(agg_cont))) < 1e-12)
  T_mat
}

# Linear blend of a proportionate and an assortative contact matrix.
# eps = 0 gives pure proportionate mixing; eps = 1 gives pure assortative.
# Operates on plain numeric matrices (not the long-format tibbles used by
# generate_matrix()).
make_blended_matrix <- function(M_PM, M_assort, eps) {
  stopifnot(identical(dim(M_PM), dim(M_assort)), eps >= 0, eps <= 1)
  (1 - eps) * M_PM + eps * M_assort
}

# combine the above to make an activity-structured matrix with n_activity_bins
# classes, heterogeneity parameter sigma, activity-assortativity range parameter
# alpha, and assortativity parameter epsilon. The 'binning_method' can either be
# "gausquad" for uneven-sized bins that contain lots of information inrelatively
# few bins, or "even" for even bin sizes, which make it easier to see the
# assortativity level.
make_activity_matrix <- function(n_activity_bins, sigma, alpha, epsilon,
                                 binning_method = c("gaussquad", "even")) {

  # choose binning function
  binning_method <- match.arg(binning_method)
  
  binner <- switch(binning_method,
                   gaussquad = gaussquad_classes,
                   even = quantile2_classes_mean)
  
  # make bins
  activity_bins <- binner(sigma = sigma,
                          n_classes = n_activity_bins)
  
  # fully assortative matrix
  assortative_activity <- make_toms_contact_matrix(
    v = activity_bins$activity,
    p_pop = activity_bins$fraction,
    alpha = alpha)
  
  # fully proportionate matrix
  proportionate_activity <- make_proportionate_contact_matrix(
    v = activity_bins$activity,
    p_pop = activity_bins$fraction)
  
  # partially-assortative matrix
  activity_matrix <- make_blended_matrix(
    M_PM = proportionate_activity,
    M_assort = assortative_activity,
    eps = epsilon)
  
  # set names
  rownames(activity_matrix) <- activity_bins$class
  colnames(activity_matrix) <- activity_bins$class
  
  # add on the activity bins and parameters as attributes
  attr(activity_matrix, "activity_bins") <- activity_bins
  attr(activity_matrix, "parameters") <- list(sigma = sigma,
                                              alpha = alpha,
                                              epsilon = epsilon)
  
  # and return
  activity_matrix
  
}

# given two matrices (e.g. an age-structured matrix and an activity matrix),
# return a combined matrix tiling the two, with each cell giving the product of
# the cell values of the two matrices. Attach the combined row/column names
tile_matrices <- function(a, b) {
  
  # get the kronecker product  
  ab <- kronecker(a, b, FUN = "*")
  
  # define the names appropriately
  a_names <- rownames(a)
  b_names <- rownames(b)
  n_a <- length(a_names)
  n_b <- length(b_names)
  
  ab_names <- paste(
    rep(a_names, each = n_b),
    rep(b_names, n_a),
    sep = "-"
  )
  
  rownames(ab) <- ab_names
  colnames(ab) <- ab_names
  
  ab
}

# is x almost exactly equal to 1 
almost_one <- function(x, tol = 1e-5) {
  max(abs(x - 1)) < tol
}

# given two vectors of population fractions (e.g. for an age-structured matrix
# and an activity matrix), return a combined vector of population fractions
# tiling the two, with each element giving the product of the elements values of
# the two vectors
tile_fractions <- function(a, b) {
  
  if (!almost_one(sum(a))) {
    stop("population fraction a does not sum to 1")
  }
  if (!almost_one(sum(b))) {
    stop("population fraction b does not sum to 1")
  }
  n_a <- length(a)
  n_b <- length(b)
  
  # tile these in the same way as in tile_matrices()
  a_tiled <- rep(a, each = n_b)
  b_tiled <- rep(b, n_a)
  ab <- a_tiled * b_tiled
  ab
}

# return a ggplot visualisation of a contact matrix - like in conmat but
# agnostic to whether it is age
autoplot_contact_matrix <- function(contact_matrix, palette = 1) {
  
  contact_matrix |>
    conmat::matrix_to_predictions() |>
    rename(
      from = age_group_from,
      to = age_group_to
    ) |>
    ggplot(
      aes(
        x = from,
        y = to,
        fill = contacts
      )
    ) +
    geom_tile() +
    coord_fixed() +
    scale_fill_distiller(
      palette = palette,
      direction = 1, 
      trans = "sqrt"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(
        size = 6,
        angle = 45,
        hjust = 1
      )
    )

}

# Run the downstream "final analysis" of analysis 3 (the two-stage method) on a
# participant-level contact dataset, returning the quantities we want bootstrap
# intervals on: the residual heterogeneity sigma, the two R0s, and the overall
# and by-age final sizes for the three scenarios (age matrix, age matrix
# rescaled to the age/activity R0, and age/activity matrix). This encapsulates
# the non-illustrative part of R/analysis_3_real_world_application.qmd so it can
# be re-run on bootstrap resamples. `contact_data` must have columns age_from,
# age_to, setting, contacts and part_id (as produced by process_fc_data_conmat()
# after renaming). All other arguments are held fixed across bootstraps, since
# they are modelling choices rather than survey data.
fc_final_size_pipeline <- function(contact_data,
                                   population,
                                   beta = 0.1,
                                   alpha = 0.5,
                                   epsilon = 0.5,
                                   n_activity_bins = 10,
                                   max_age = 90) {

  # --- step a: age-structured contact model (aggregated across participants) --
  contact_data_summarised <- contact_data |>
    group_by(age_from, age_to, setting) |>
    summarise(contacts = sum(contacts),
              participants = n(),
              .groups = "drop")

  model_aggregated <- fit_single_contact_model(
    contact_data = contact_data_summarised,
    population = population
  )

  # 5-year age matrix up to max_age, and the population fraction in each bin
  age_breaks <- seq(0, max_age, by = 5)
  age_contact_pred <- predict_contacts(model_aggregated,
                                       population,
                                       age_breaks = age_breaks)
  age_matrix <- conmat::predictions_to_matrix(age_contact_pred)
  pop_bin <- population |>
    filter(lower.age.limit < max_age) |>
    pull(population)
  age_fractions <- pop_bin / sum(pop_bin)

  # --- step b: residual between-participant SD (sigma) via a Poisson GLMM ------
  # predict per-year contacts to use as the model offset
  year_range <- range(c(contact_data$age_from, contact_data$age_to),
                      na.rm = TRUE)
  age_breaks_years <- seq(year_range[1], year_range[2] + 1, by = 1)
  year_age_group_lookup <- tibble(age = age_breaks_years) |>
    mutate(age_group = sprintf("[%s,%s)", age, age + 1))

  age_contact_pred_years <- model_aggregated |>
    predict_contacts(population, age_breaks = age_breaks_years) |>
    left_join(year_age_group_lookup, by = c(age_group_from = "age_group")) |>
    rename(age_from = age) |>
    left_join(year_age_group_lookup, by = c(age_group_to = "age_group")) |>
    rename(age_to = age, predicted_contacts = contacts) |>
    select(age_from, age_to, predicted_contacts)

  contact_data_modelling <- contact_data |>
    left_join(age_contact_pred_years, by = c("age_from", "age_to")) |>
    filter(age_from <= max_age, age_to <= max_age) |>
    mutate(obs = factor(row_number()))

  # nAGQ = 0 for numerical robustness (see the qmd for the rationale)
  model_twostage <- lme4::glmer(
    contacts ~ offset(log(predicted_contacts)) + (1 | part_id) + (1|obs),
    family = stats::poisson,
    data = contact_data_modelling,
    nAGQ = 0
  )

  sigma <- model_twostage |>
    summary() |>
    extracting("varcor") |>
    as_tibble() |>
    filter(grp == "part_id") |>
    pull(sdcor)

  # --- activity-structured matrix and combined age/activity matrix ------------
  activity_matrix <- make_activity_matrix(n_activity_bins = n_activity_bins,
                                          sigma = sigma,
                                          alpha = alpha,
                                          epsilon = epsilon)
  activity_fractions <- attr(activity_matrix, "activity_bins")$fraction

  age_activity_matrix <- tile_matrices(age_matrix, activity_matrix)
  age_activity_fractions <- tile_fractions(age_fractions, activity_fractions)

  # --- R0s --------------------------------------------------------------------
  R0_age <- get_eigenval(age_matrix * beta)
  R0_age_activity <- get_eigenval(age_activity_matrix * beta)

  # --- final sizes for the three scenarios ------------------------------------
  n_age <- length(age_fractions)
  n_age_activity <- length(age_activity_fractions)
  dummy_age <- as.matrix(rep(1, n_age))
  dummy_age_activity <- as.matrix(rep(1, n_age_activity))

  # rescale so that (matrix * fractions) has an eigenvalue of 1 (a requirement
  # of final_size)
  age_matrix_scaled <- age_matrix / get_eigenval(age_matrix * age_fractions)
  age_activity_matrix_scaled <- age_activity_matrix /
    get_eigenval(age_activity_matrix * age_activity_fractions)

  size_age <- final_size(r0 = R0_age,
                         contact_matrix = age_matrix_scaled,
                         demography_vector = age_fractions,
                         susceptibility = dummy_age,
                         p_susceptibility = dummy_age)

  size_age_scaled <- final_size(r0 = R0_age_activity,
                                contact_matrix = age_matrix_scaled,
                                demography_vector = age_fractions,
                                susceptibility = dummy_age,
                                p_susceptibility = dummy_age)

  size_age_activity <- final_size(r0 = R0_age_activity,
                                  contact_matrix = age_activity_matrix_scaled,
                                  demography_vector = age_activity_fractions,
                                  susceptibility = dummy_age_activity,
                                  p_susceptibility = dummy_age_activity)

  # collapse the age/activity final sizes down to age groups
  size_age_activity_by_age <- size_age_activity |>
    mutate(demo_grp = str_split_i(demo_grp, "-", 1)) |>
    group_by(demo_grp) |>
    summarise(p_infected = sum(p_infected * activity_fractions),
              .groups = "drop")

  # overall (population-weighted) final sizes
  overall <- c(
    age          = sum(size_age$p_infected * age_fractions),
    age_scaled   = sum(size_age_scaled$p_infected * age_fractions),
    age_activity = sum(size_age_activity$p_infected * age_activity_fractions)
  )

  # tidy by-age final sizes for the three scenarios, with numeric bin edges
  parse_bins <- function(df) {
    df |>
      mutate(
        age_lower = as.numeric(str_match(demo_grp, "([0-9]+)[^0-9]+([0-9]+)")[, 2]),
        age_upper = as.numeric(str_match(demo_grp, "([0-9]+)[^0-9]+([0-9]+)")[, 3])
      ) |>
      select(age_lower, age_upper, p_infected)
  }

  by_age <- bind_rows(
    age          = parse_bins(size_age),
    age_scaled   = parse_bins(size_age_scaled),
    age_activity = parse_bins(size_age_activity_by_age),
    .id = "scenario"
  ) |>
    arrange(scenario, age_lower)

  list(
    sigma = sigma,
    R0_age = R0_age,
    R0_age_activity = R0_age_activity,
    overall = overall,
    by_age = by_age
  )
}

# Draw a cluster (participant-level) bootstrap resample from a participant-level
# contact dataset: sample the participants with replacement and give each drawn
# participant a fresh unique part_id, so that a participant drawn more than once
# is treated as distinct individuals in the (1 | part_id) random effect (rather
# than as one participant with repeated observations, which would bias sigma
# downward). `canonical_ids` fixes the participant ordering so the returned
# membership indices are comparable across resamples (needed for a later BCa
# acceleration term via jackknife-after-bootstrap). Returns the resampled data
# and the integer indices (into canonical_ids) of the drawn participants.
resample_participants <- function(contact_data, canonical_ids) {
  n <- length(canonical_ids)
  drawn_idx <- sample.int(n, size = n, replace = TRUE)
  key <- tibble(
    new_id = paste0("boot", seq_len(n)),
    part_id = canonical_ids[drawn_idx]
  )
  resampled <- key |>
    left_join(contact_data, by = "part_id", relationship = "many-to-many") |>
    mutate(part_id = new_id) |>
    select(-new_id)
  list(data = resampled, membership = drawn_idx)
}

# Bootstrap confidence interval for a scalar quantity, returning both the plain
# percentile interval and the bias-corrected ("BC") interval. `theta_star` is
# the vector of bootstrap replicate estimates and `theta_hat` the point estimate
# from the full data. The BC interval shifts the percentiles to correct for any
# median bias between theta_hat and the bootstrap distribution, via the
# bias-correction constant z0 = Phi^-1(fraction of replicates below theta_hat).
# The acceleration term of a full BCa interval is set to zero here (see the
# analysis-3 qmd for why it is deferred); z0 is returned so the size of the bias
# correction can be inspected. Following Efron & Tibshirani (1993).
bc_ci <- function(theta_star, theta_hat, level = 0.95) {
  theta_star <- theta_star[is.finite(theta_star)]
  n_star <- length(theta_star)
  a_tail <- (1 - level) / 2
  probs <- c(a_tail, 1 - a_tail)

  # plain percentile interval
  pct <- quantile(theta_star, probs, names = FALSE, type = 7)

  # bias correction: z0 from the fraction of replicates below the point estimate
  prop_below <- mean(theta_star < theta_hat)
  # guard degenerate proportions (0 or 1 would give an infinite z0)
  prop_below <- min(max(prop_below, 1 / (2 * n_star)), 1 - 1 / (2 * n_star))
  z0 <- qnorm(prop_below)

  # adjusted percentiles (acceleration a = 0)
  adj_probs <- pnorm(2 * z0 + qnorm(probs))
  bc <- quantile(theta_star, adj_probs, names = FALSE, type = 7)

  tibble(
    estimate = theta_hat,
    pct_lower = pct[1], pct_upper = pct[2],
    bc_lower = bc[1], bc_upper = bc[2],
    z0 = z0,
    n_boot = n_star
  )
}

# Sample skewness (the standardised third moment g1) of a vector, used to gauge
# how far a bootstrap distribution departs from symmetry. Non-finite values are
# dropped.
sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- sqrt(mean((x - m)^2))
  mean((x - m)^3) / s^3
}

# BCa acceleration constant(s) via jackknife-after-bootstrap (Efron 1992). This
# is the quantity that distinguishes a full BCa interval from the plain BC
# interval above: a = (1/6) * skewness of the empirical influence values. Given a
# matrix (or data frame) of bootstrap replicate estimates -- one column per
# quantity, one row per replicate -- and `membership`, the list of resampled unit
# indices for each replicate (aligned row-for-row with `theta_star`), it
# estimates each column's acceleration without any re-fitting: the
# leave-one-unit-out estimate theta_(j) is approximated by the mean of the
# replicates in which unit j did not appear, and
#   a = sum_j (theta_(.) - theta_(j))^3 / (6 * (sum_j (theta_(.) - theta_(j))^2)^{3/2}),
# with theta_(.) the mean of the leave-one-out estimates. Returns one
# acceleration per column of `theta_star`.
bca_acceleration <- function(theta_star, membership, n_units) {
  theta_star <- as.matrix(theta_star)
  B <- nrow(theta_star)
  # accumulate, per unit, the sum and count of replicates in which it is present
  present_sum <- matrix(0, n_units, ncol(theta_star))
  present_cnt <- integer(n_units)
  for (r in seq_len(B)) {
    j <- unique(membership[[r]])
    present_sum[j, ] <- present_sum[j, ] + rep(theta_star[r, ], each = length(j))
    present_cnt[j]   <- present_cnt[j] + 1L
  }
  # leave-one-out mean = mean over the replicates in which the unit was absent
  absent_sum <- sweep(-present_sum, 2, colSums(theta_star), `+`)
  theta_loo  <- absent_sum / (B - present_cnt)
  apply(theta_loo, 2, function(tl) {
    tl <- tl[is.finite(tl)]        # units present in every replicate contribute nothing
    d <- mean(tl) - tl
    sum(d^3) / (6 * sum(d^2)^1.5)
  })
}

# Full BCa confidence interval given a precomputed acceleration `a` (see
# bca_acceleration()). With a = 0 this reduces exactly to the BC interval of
# bc_ci(); a non-zero `a` additionally corrects for skewness. Returns the
# interval endpoints alongside the acceleration used.
bca_ci <- function(theta_star, theta_hat, a, level = 0.95) {
  theta_star <- theta_star[is.finite(theta_star)]
  n_star <- length(theta_star)
  a_tail <- (1 - level) / 2
  z <- qnorm(c(a_tail, 1 - a_tail))

  # bias correction z0, guarded against degenerate proportions as in bc_ci()
  prop_below <- mean(theta_star < theta_hat)
  prop_below <- min(max(prop_below, 1 / (2 * n_star)), 1 - 1 / (2 * n_star))
  z0 <- qnorm(prop_below)

  adj_probs <- pnorm(z0 + (z0 + z) / (1 - a * (z0 + z)))
  bca <- quantile(theta_star, adj_probs, names = FALSE, type = 7)

  tibble(
    estimate = theta_hat,
    bca_lower = bca[1], bca_upper = bca[2],
    z0 = z0, a = a
  )
}
