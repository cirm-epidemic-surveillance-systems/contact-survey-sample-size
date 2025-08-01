# load and process the French Connection data for the variance partitioning
# analysis
process_fc_data_varpart <- function() {
  
  # participant info
  fc_participants <- read_csv("data/french_connection/2015_Beraud_France_participant_common.csv")
  fc_participants_extra <- read_csv("data/french_connection/2015_Beraud_France_participant_extra.csv")
  
  # get all observed combinations of participant and wave, to pad contacts with 0s
  fc_participants_observed <- fc_participants_extra %>%
    select(part_id,
           wave,
           participant_occupation) %>%
    rename(
      part_occupation = participant_occupation 
    ) %>%
    expand_grid(
      studyDay = 1:2
    ) %>%
    left_join(
      fc_participants,
      by = "part_id"
    ) %>%
    select(
      -hh_id
    )
  
  # contact events
  fc_contacts <- read_csv("data/french_connection/2015_Beraud_France_contact_common.csv")
  fc_contacts_extra <- read_csv("data/french_connection/2015_Beraud_France_contact_extra.csv")
  
  # collapse contact events to count contacts per participant, per wave, per
  # studyDay
  fc_contacts_sry <- fc_contacts %>%
    # add on wave and day columns
    left_join(fc_contacts_extra,
              by = "cont_id") %>%
    # count contacts per participant/wave/day
    group_by(
      part_id,
      wave,
      studyDay
    ) %>%
    summarise(
      contacts = n(),
      .groups = "drop"
    )
  
  # add observed contacts, fill in others with 0s
  fc_participants_observed %>%
    left_join(fc_contacts_sry,
              by = c("part_id", "wave", "studyDay")) %>%
    mutate(
      contacts = replace_na(contacts, 0)
    ) %>%
    # set factors
    mutate(
      part_id = factor(part_id),
      part_age = factor(part_age),
      part_gender = factor(part_gender),
      part_occupation = factor(part_occupation),
      wave = factor(wave),
      studyDay = factor(studyDay)
    ) %>%
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
  
  filepath <- "data/hongkong/HongKongData.csv"
  
  filepath %>%
    read_csv() %>%
    select(
    part_id = pid,
    part_age = age,
    part_gender = sex,
    contacts = n.contact.total
  ) %>% 
    mutate(
      part_id = as_factor(part_id),
      part_age = as_factor(part_age),
      part_gender = as_factor(part_gender)
    )
  
}

# given a fitted random effects model, perform the variance partitioning
partition_variance_lmer <- function (model) {
  model %>%
    summary() %>%
    `[[`("varcor") %>%
    as_tibble() %>%
    mutate(
      var = sdcor ^ 2,
      proportion = var / sum(var)
    ) %>%
    select(
      grp,
      var,
      proportion
    ) %>%
    mutate(
      grp = case_when(
        grp == "part_age" ~ "between ages",
        grp == "part_gender" ~ "between genders",
        grp == "part_occupation" ~ "between occupations",
        grp == "Residual" ~ "within individuals",
        .default = "between individuals (unexplained)",
      )
    ) %>%
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
  
  partitioning %>%
    # set the factor order for partition to how we want to plot it
    mutate(
      partition = factor(partition,
                         levels = names(colour_palette))
    ) %>%
    # sort by partition (to the factor order), so we can place the labels in the
    # right places
    arrange(
      partition
    ) %>%
    group_by(
      Study
    ) %>%
    # for labels
    mutate(
      text_percent = scales::label_percent()(proportion),
      text_position = rev(cumsum(rev(proportion))) - proportion / 2
    ) %>%
    ungroup() %>%
    # make some better names
    rename(
      `Variance explained` = proportion,
      Component = partition
    ) %>%
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
  
  # participant info
  fc_participants <- read_csv("data/french_connection/2015_Beraud_France_participant_common.csv")
  fc_participants_extra <- read_csv("data/french_connection/2015_Beraud_France_participant_extra.csv")
  
  # get all observed combinations of participant and wave, to pad contacts with 0s
  fc_participants_observed <- fc_participants_extra %>%
    select(part_id,
           wave) %>%
    expand_grid(
      studyDay = 1:2
    ) %>%
    left_join(
      fc_participants,
      by = "part_id"
    ) %>%
    select(
      -hh_id,
      -part_gender
    )
  
  # contact events
  fc_contacts <- read_csv("data/french_connection/2015_Beraud_France_contact_common.csv")
  fc_contacts_extra <- read_csv("data/french_connection/2015_Beraud_France_contact_extra.csv")
  
  # collapse contact events to count contacts per participant, per wave, per
  # studyDay, per contact age
  fc_contacts_sry <- fc_contacts %>%
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
    ) %>%
    # add on wave and day columns
    left_join(fc_contacts_extra,
              by = "cont_id") %>%
    # count contacts per participant/wave/day
    group_by(
      part_id,
      wave,
      studyDay,
      cont_age
    ) %>%
    summarise(
      contacts = n(),
      .groups = "drop"
    )
  
  # get the range of contact ages, and pad the observations with explicit 0s
  possible_ages <- seq(
    from = min(fc_contacts_sry$cont_age),
    to = max(fc_contacts_sry$cont_age),
    by = 1)
  
  fc_contacts_complete <- fc_contacts_sry %>%  
    complete(
      # only the participant/wave/studyDay combinations in the data
      nesting(part_id, wave, studyDay),
      # but for all possible ages
      cont_age = possible_ages,
      # if the value is not int he contacts file, there must be 0 contacts
      fill = list(contacts = 0)
    )
    
  # add observed contacts, fill in others with 0s
  fc_participants_observed %>%
    left_join(fc_contacts_complete,
              by = c("part_id", "wave", "studyDay")) %>%
    mutate(
      contacts = replace_na(contacts, 0)
    ) %>%
    # set factors
    mutate(
      part_id = factor(part_id),
      wave = factor(wave),
      studyDay = factor(studyDay)
    ) %>%
    relocate(
      part_age,
      .before = cont_age
    )
}

#supporting function for quantizing
MeasureMeanByContactQuantile<-function(contacts,N_quantiles){
  contacts |> 
    mutate(MeansByQuantile = ntile(X,N_quantiles)) |> 
    group_by(MeansByQuantile) |> 
    summarise(MeansByQuantile = mean(X)) |> 
    mutate(Quantile = row_number())->tmp
  return(tmp)
}

GenerateProportionateMixingMatrix <-function(s){
  SumDegree = sum(s$MeansByQuantile)
  cross_join(s,s)|> 
    mutate(Product = MeansByQuantile.x*MeansByQuantile.y)|> 
    select(From = Quantile.x, To = Quantile.y,Product)|> 
    mutate(Product = Product/SumDegree)|> 
    rename(Weight = Product)->tmp
  
  return(tmp)
}

GenerateFullyAssortativeMixingMatrix <-function(s){
  
  cross_join(s,s)|>
    rename(From = Quantile.x, To = Quantile.y)|> 
    mutate(Weight = ifelse(From == To, MeansByQuantile.x,0 )) |> 
    select(From,To,Weight)->tmp
  return(tmp)    
}

#wrapper function

GenerateMatrix<-function(sigma = 1,number_of_quantiles = 2, assort=1){
  #sigma is the normal distribution standard deviation
  #alpha is the amount of assortativity
  as_tibble(x=exp(rnorm(n=1000000))) |> 
    rename(X=value)-> normal_dist_test.df
  
  MeasureMeanByContactQuantile(normal_dist_test.df,number_of_quantiles)->s
  
  AssortativeMatrix<-GenerateFullyAssortativeMixingMatrix(s)
  ProportionateMatrix<-GenerateProportionateMixingMatrix(s)
  
  AssortativeMatrix
  ProportionateMatrix
  
  AssortativeMatrix |> 
    rename(Weight_Assortative=Weight) |> 
    mutate(Weight_Proportionate = ProportionateMatrix$Weight) |> 
    mutate(WeightCombined = assort*Weight_Assortative + (1-assort)*Weight_Proportionate) |> 
    select(From,To,WeightCombined) -> CombinedMatrix
  
  return(CombinedMatrix)
}

MatrixToEigenvalue<-function(M,eigen_value = 1){
  #M is the matrix in long format, eigen_value is which eigenvalue to extract 
  #defaults to eigen_value = 1 for the dominant 
  M |> usedist::pivot_to_numeric_matrix(From, To, WeightCombined) |> 
    eigen()->decomp 
  tmp = decomp$values[eigen_value]
  return(tmp)  
}

MapToEigen<-function(assort = 1,number_of_quantiles = 3,eigen_value=1,beta = 1){
  GenerateMatrix(assort = assort,number_of_quantiles = number_of_quantiles)->M 
  M |> mutate(WeightCombined = WeightCombined*beta)
  MatrixToEigenvalue(M,eigen_value=eigen_value)->tmp
  return(tmp)
}

# sanity check contact variance analysis by permuting the participant IDs
permute_sim <- function() {
  
  fc_contact_counts %>%
    mutate(
      contacts = sample(contacts)
    ) %>%
    lmer(
      log1p(contacts) ~ (1|part_age) + (1|part_id),
      data = .
    ) %>%
    partition_variance_lmer() %>%
    select(-var) %>%
    pivot_wider(
      names_from = "partition",
      values_from = "proportion"
    )
  
}

# make a proportionate mixing contact matrix from class contact rates
make_contact_matrix <- function(class_means) {
  class_means %*% t(class_means) / sum(class_means)
}

# get the dominant eigenvalue of a matrix
get_eigenval <- function(matrix) {
  Re(eigen(matrix)$values[1])
}