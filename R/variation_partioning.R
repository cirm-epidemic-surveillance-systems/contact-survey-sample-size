# Analyse partitioning of variation in contact rates (within/between individuals
# and strata)

library(tidyverse)

# load relevant French Connection data

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
fc_contact_counts <- fc_participants_observed %>%
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

# analyse this with a random effects model to capture the within-individual
# variation, the between-individual variation explained by age, and the residual
# between-individual variation
library(lme4)

# Gaussian random effects model (not overfitted, unreasonable distribution?)
model_lmer <- lmer(
  contacts ~ (1|part_age) +
    (1|part_gender) +
    (1|part_occupation) +
    (1|part_id),
  data = fc_contact_counts
)

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

partioning <- partition_variance_lmer(model_lmer)


# add a barplot of this
colour_table <- tribble(
  ~partition, ~colour,
  "within individuals", "pink",
  "between individuals (unexplained)", scales::alpha("blue", 0.1),
  "between genders", scales::alpha("blue", 0.4),
  "between ages", scales::alpha("blue", 0.5),
  "between occupations", scales::alpha("blue", 0.6)
)

colour_vector <- colour_table %>%
  pivot_wider(
    names_from = partition,
    values_from = colour
  ) %>%
  as.vector() %>%
  do.call(c, .)

partioning %>%
  mutate(
    partition = factor(partition,
                       levels = colour_table$partition),
    study = "France"
  ) %>%
  arrange(
    partition
  ) %>%
  # for labels
  mutate(
    text_percent = scales::label_percent()(proportion),
    text_position = rev(cumsum(rev(proportion))) - proportion / 2
  ) %>%
  rename(
    `Variance explained` = proportion,
    Component = partition
  ) %>%
  ggplot(
    aes(
      x = study,
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
    # vjust = 0.5
  ) +
  scale_y_continuous(
    labels = scales::label_percent()
  ) +
  scale_fill_manual(values = colour_vector) +
  theme_minimal()

# sanity check this analysis by permuting the participant IDs
permute_sim <- function() {
  
  fc_contact_counts %>%
    mutate(
      contacts = sample(contacts)
    ) %>%
    lmer(
      contacts ~ (1|part_age) + (1|part_id),
      data = .
    ) %>%
    partition_variance_lmer() %>%
    select(-var) %>%
    pivot_wider(
      names_from = "partition",
      values_from = "proportion"
    )
  
}

sims_list <- replicate(50,
                  permute_sim(),
                  simplify = FALSE)
sims_df <- do.call(bind_rows, sims_list)

partitions <- colnames(partioning)
par(mfrow = c(length(partitions), 1))

for (partition in partitions) {
  perturbed_props <- sims_df[[partition]]
  estimated_prop <- partioning[[partition]]
  
  xlim <- range(c(perturbed_props, estimated_prop))
  
  hist(perturbed_props,
       xlim = xlim,
       breaks = 100,
       xlab = "Proportion",
       main = partition)
  abline(v = estimated_prop)
}

# add participant characteristics to model to explain remaining between
# individual variance in contacts
