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

# extract the standard deviation of the participant random effect

# output the discretised age matrix and the discretised social activity matrix

# combine the two matrices together (kronecker product)

