# Fit a modified conmat model to the French Connection dataset, including
# participant ID as a random effect

library(tidyverse)

source("R/functions.R")

library(conmat)

# Load FC data with the contact ages, and rename things to those expected by
# conmat
fc_contact_data <- process_fc_data_conmat() %>%
  rename(age_from = part_age,
         age_to = cont_age) %>%
  mutate(
    setting = "all"
  )

# aggregate observations across all participants (this won't work for our social
# activity version, but doing this for now to check conmat analysis is working
# on these data)
fc_contact_data_summarised <- fc_contact_data %>%
  group_by(
    age_from,
    age_to,
    setting) %>%
  summarise(
    contacts = sum(contacts),
    participants = n(),
    .groups = "drop"
  )

# get population data for France in 2010
fc_population <- age_population(
  data = socialmixr::wpp_age(),
  location_col = country,
  location = c("France"),
  age_col = lower.age.limit,
  year_col = year,
  year = 2010
)

# run conmat on this and return the age-structured matrices
fc_model_aggregated <- fit_single_contact_model(
  contact_data = fc_contact_data_summarised,
  population = fc_population
)

age_breaks <- c(seq(0, 75, by = 5), Inf)
age_contact_pred <- predict_contacts(fc_model_aggregated,
                                     fc_population,
                                     age_breaks = age_breaks)

age_matrix <- conmat::predictions_to_matrix(age_contact_pred)
autoplot(age_matrix)

# modify the fitting function to add a random effect (create a conmat branch)
fc_contact_data_participant <- fc_contact_data %>%
  select(
    age_from,
    age_to,
    setting,
    contacts,
    part_id
  ) %>%
  mutate(
    participants = 1
  )

# use the modified version of conmat to add this term to the model

# because we can't aggregate over participants, we need to aggregate the ages to
# reduce the dataset size and run time.

# aggregate to the midpoints of the discrete age classes we will be predicting to
# later
max_age <- max(fc_contact_data_participant$age_from,
               fc_contact_data_participant$age_to,
               na.rm = TRUE)

age_agg_lookup <- tibble(
  lower = age_breaks[-length(age_breaks)],
  upper = age_breaks[-1]
  ) %>%
  mutate(
    label = sprintf("[%s,%s)", lower, upper),
    upper = case_when(
      !is.finite(upper) ~ max_age + 1,
      .default = upper),
    midpoint = lower + 2
  ) %>%
  rowwise() %>%
  mutate(
    age = list(seq(lower, upper - 1, by = 1))
  ) %>%
  unnest(age) %>%
  select(
    age,
    label,
    midpoint
  )

# fc_contact_data_participant_lores <- fc_contact_data_participant %>%
#   left_join(
#     age_agg_lookup,
#     by = "age_to"
#   ) %>%
#   group_by(
#     part_id,
#     age_from,
#     setting,
#     age_to_midpoint
#   ) %>%
#   summarise(
#     contacts = sum(contacts),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     participants = 1
#   ) %>%
#   rename(
#     age_to = age_to_midpoint
#   )

# two-stage inference for age structure and leftover between-individual
# variation:

# fit the age-only model using only conmat
# model the residuals from this with a random effects model

# join the age-structured matrix estimates onto this dataset

fc_contact_data_participant %>%
  left_join(
    age_agg_lookup,
    by = c(age_from = "age")
  ) %>%
  rename(
    age_group_from = label,
    age_from_midpoint = midpoint
  ) %>%
  left_join(
    age_agg_lookup,
    by = c(age_to = "age")
  ) %>%
  rename(
    age_group_to = label,
    age_to_midpoint = midpoint
  ) %>%
  left_join(
    age_contact_pred,
    by = c(?)
  )



# fit a random effects poisson regression, with these (logged) as an offset



# # run conmat on this and return the age-structured matrix
# fc_model_participant <- fit_single_contact_model(
#   contact_data = fc_contact_data_participant_lores,
#   population = fc_population
# )
# 
# # save this model!
# saveRDS(fc_model_participant,
#         "temporary/fc_model_participant.RDS")
# 
# # Produce the usual age-structured matrix
# 
# # need to hack this to handle part id!
# age_matrix_participant <- fc_model_participant %>%
#   predict_contacts(
#     population = fc_population,
#     age_breaks = age_breaks
#   ) %>%
#   predictions_to_matrix()
# 
# # output the discretised age matrix
# autoplot(age_matrix_participant)

# extract the standard deviation of the participant random effect

# like this? Need to find the right components with grepl on part_id
sd(fc_model_participant$coefficients[-1])

fc_model_participant$smooth

# create the discretised social activity matrix
# add row and column names

# combine the two matrices together (kronecker product)
contact_matrix <- kronecker(social_activity_matrix,
                            age_matrix_participant,
                            FUN = "*")

