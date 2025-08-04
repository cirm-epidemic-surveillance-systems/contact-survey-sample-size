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

generate_proportionate_mixing_matrix <- function(s) {
  sum_degree <- sum(s$means_by_quantile)
  cross_join(s, s) |> 
    mutate(
      product = means_by_quantile.x * means_by_quantile.y
    ) |> 
    select(
      from = quantile.x,
      to = quantile.y, product
    ) |> 
    mutate(
      product = product / sum_degree
    ) |> 
    rename(
      weight = product
    )
}

generate_fully_assortative_mixing_matrix <- function(s) {
  cross_join(s, s) |>
    rename(
      from = quantile.x,
      to = quantile.y
    ) |> 
    mutate(
      weight = ifelse(from == to,
                      means_by_quantile.x,
                      0)
    ) |> 
    select(
      from,
      to,
      weight
    )
}

# wrapper function
# sigma is the normal distribution standard deviation
# alpha is the amount of assortativity
generate_matrix <- function(sigma = 1, number_of_quantiles = 2, assort = 1) {
  
  normal_dist_test.df <- as_tibble(
    x = exp(rnorm(n = 1000000))
  ) |> 
    rename(
      X = value
    )
  
  s <- measure_mean_by_contact_quantile(normal_dist_test.df,
                                        number_of_quantiles)
  
  assortative_matrix <- generate_fully_assortative_mixing_matrix(s)
  proportionate_matrix <- generate_proportionate_mixing_matrix(s)
  
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
  decomp$values[eigen_value]
}

map_to_eigen <- function(assort = 1,
                         number_of_quantiles = 3,
                         eigen_value = 1,
                         beta = 1) {
  M <- generate_matrix(assort = assort,
                       number_of_quantiles = number_of_quantiles)
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