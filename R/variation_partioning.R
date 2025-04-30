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
    wave = factor(wave),
    studyDay = factor(studyDay)
  ) %>%
  relocate(
    part_age,
    .after = part_id
  )

# fc_contact_counts %>%
#   arrange(
#     part_id, wave, studyDay
#   ) %>%
#   view()


# analyse this with a random effects model to capture the within-individual
# variation, the between-individual variation explained by age, and the residual
# between-individual variation
library(lme4)
# 
# # Variance (normal assumption) version (overfitted)
# model_lm <- lm(
#   contacts ~ part_age + part_id,
#   data = fc_contact_counts
# )
# 
# sry_lm <- anova(model_lm)
# 
# variances <- c(
#   between_age = sry_lm$`Sum Sq`[1],
#   between_other = sry_lm$`Sum Sq`[2],
#   within = sry_lm$`Sum Sq`[3])
# r2_lm <- variances / sum(variances)
# r2_lm
# 
# # Poisson GLM version (overfitted) 
# model_glm <- glm(
#   contacts ~ part_age + part_id,
#   family = stats::poisson,
#   data = fc_contact_counts
# )
# 
# sry_glm <- anova(model_glm)
# deviances <- c(
#   between_age = sry_glm$Deviance[2],
#   between_other = sry_glm$Deviance[3],
#   within = sry_glm$`Resid. Dev`[3])
# r2_glm <- deviances / sum(deviances)
# r2_glm
# 

# Gaussian random effects model (not overfitted, unreasonable distribution?)
model_lmer <- lmer(
  contacts ~ (1|part_age) + (1|part_id),
  data = fc_contact_counts
)

sry_lmer <- summary(model_lmer)
variances_df <- as.data.frame(sry_lmer$varcor)
variances_lmer <- c(
  between_age = variances_df$sdcor[2] ^ 2,
  between_other = variances_df$sdcor[1] ^ 2,
  within = variances_df$sdcor[3] ^ 2)
r2_lmer <- variances_lmer / sum(variances_lmer)
r2_lmer
r2_lmer_between <- r2_lmer[1:2]
r2_lmer_between <- r2_lmer_between / sum(r2_lmer_between)
r2_lmer_between


# log1p-Gaussian random effects model (not overfitted, more reasonable distribution)
model_lmer_log <- lmer(
  log1p(contacts) ~ (1|part_age) + (1|part_id),
  data = fc_contact_counts
)

sry_lmer_log <- summary(model_lmer_log)
variances_df_log <- as.data.frame(sry_lmer_log$varcor)
variances_lmer_log <- c(
  between_age = variances_df_log$sdcor[2] ^ 2,
  between_other = variances_df_log$sdcor[1] ^ 2,
  within = variances_df_log$sdcor[3] ^ 2)
r2_lmer_log <- variances_lmer_log / sum(variances_lmer_log)
r2_lmer_log

r2_lmer_log_between <- r2_lmer_log[1:2]
r2_lmer_log_between <- r2_lmer_log_between / sum(r2_lmer_log_between)
r2_lmer_log_between


# sanity check this analysis by permuting the participant IDs
permute_sim <- function() {
  
  fc_contact_counts_permuted <- fc_contact_counts %>%
    mutate(contacts = sample(contacts))
  
  model_lmer <- lmer(
    contacts ~ (1|part_age) + (1|part_id),
    data = fc_contact_counts_permuted
  )
  
  sry_lmer <- summary(model_lmer)
  variances_df <- as.data.frame(sry_lmer$varcor)
  variances_lmer <- c(
    between_age = variances_df$sdcor[2] ^ 2,
    between_other = variances_df$sdcor[1] ^ 2,
    within = variances_df$sdcor[3] ^ 2)
  r2_lmer <- variances_lmer / sum(variances_lmer)
  r2_lmer
}

sims <- replicate(1000, permute_sim())

hist(sims["between_age", ],
     xlim = c(0, 0.1),
     breaks = 100)
abline(v= r2_lmer["between_age"])

hist(sims["between_other", ],
     xlim = c(0, 0.5),
     breaks = 100)
abline(v= r2_lmer["between_other"])

hist(sims["within", ],
     xlim = c(0, 1),
     breaks = 100)
abline(v= r2_lmer["within"])
