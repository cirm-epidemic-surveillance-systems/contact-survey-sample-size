# Analyse partitioning of variation in contact rates (within/between individuals
# and strata)

# given a fitted random effects model, perfrom the variance partitioning
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
    # for labels
    mutate(
      text_percent = scales::label_percent()(proportion),
      text_position = rev(cumsum(rev(proportion))) - proportion / 2
    ) %>%
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

library(tidyverse)
library(patchwork)

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

# Gaussian random effects model, with only age, and between vs within individual
# effects
fc_model_lmer_none <- lmer(
  contacts ~ (1|part_id),
  data = fc_contact_counts
)

fc_model_lmer_age <- lmer(
  contacts ~ (1|part_age) +
    (1|part_id),
  data = fc_contact_counts
)

# the same but with gender and occupation (confounded)
fc_model_lmer_all <- lmer(
  contacts ~ (1|part_age) +
    (1|part_gender) +
    (1|part_occupation) +
    (1|part_id),
  data = fc_contact_counts
)

# do variance partitioning on each, combining across studies
partitioning_none <- bind_rows(
  France = partition_variance_lmer(fc_model_lmer_none),
  .id = "Study"
)

partitioning_age <- bind_rows(
  France = partition_variance_lmer(fc_model_lmer_age),
  .id = "Study"
)

partitioning_all <- bind_rows(
  France = partition_variance_lmer(fc_model_lmer_all),
  .id = "Study"
)

barplot_none <- make_barplot(partitioning_none) +
  ggtitle(label = "no covariates")
barplot_age <- make_barplot(partitioning_age) +
  ggtitle(label = "age only")
barplot_all <- make_barplot(partitioning_all) +
  ggtitle(label = "age, gender, and occupation")


(barplot_none / barplot_age / barplot_all) +
  plot_annotation(title = "Partitioning variance in daily contact rates")

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

partitions <- partitioning_age$partition
par(mfrow = c(length(partitions), 1))

for (this_partition in partitions) {
  perturbed_props <- sims_df[[this_partition]]
  estimated_prop <- partitioning_age %>%
    filter(partition == this_partition) %>%
    pull(proportion)
  
  xlim <- range(c(perturbed_props, estimated_prop))
  
  hist(perturbed_props,
       xlim = xlim,
       breaks = 100,
       xlab = "Proportion",
       main = this_partition)
  abline(v = estimated_prop)
}
