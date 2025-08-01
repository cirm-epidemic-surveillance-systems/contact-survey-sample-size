# Analyse partitioning of variation in contact rates (within/between individuals
# and strata)

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the processed French Connection data
fc_contact_counts <- process_fc_data_varpart()

# and the Hong Kong data
hk_contact_counts <- process_hk_data_varpart()

# analyse this with a random effects model to capture the within-individual
# variation, the between-individual variation explained by age, and the residual
# between-individual variation


# model age first,then add estimates as an offset 

fc_model_lmer_only_age <- lmer(
  log1p(contacts) ~ (1|part_age),
  data = fc_contact_counts
)

age_offset_lookup <- fc_model_lmer_only_age |>
  coef |>
  `$`(part_age) |>
  tibble(
    part_age = rownames(.),
    offset = .[, 1]
  ) |>
  select(-`(Intercept)`)

fc_contact_counts_age_est <- fc_contact_counts |>
  left_join(age_offset_lookup,
            by = "part_age")

fc_model_lmer_part_given_age <- lmer(
  log1p(contacts) ~ (1|part_id) + offset(offset),
  data = fc_contact_counts_age_est
)

fc_model_lmer_only_age
# 0.2154 sd age
# 0.6929 sd other
fc_model_lmer_part_given_age
# (of the other)
# 0.4181 sd participant
# 0.5478 participant


age_frac <- (0.2154^2) / (0.2154^2 + 0.6929^2)
non_age_frac <- 1 - age_frac

part_frac_of_other <- (0.4181^2) / (0.4181^2 + 0.5478^2)
non_part_frac_of_other <- 1 - part_frac_of_other

part_frac <- non_age_frac * part_frac_of_other
within_frac <- non_age_frac * non_part_frac_of_other

within_frac # 58%
part_frac # 33%
age_frac # 9%

# Gaussian random effects model, with only age, and between vs within individual
# effects
fc_model_lmer_none <- lmer(
  log1p(contacts) ~ (1|part_id),
  data = fc_contact_counts
)

fc_model_lmer_age <- lmer(
  log1p(contacts) ~ (1|part_age) +
    (1|part_id),
  data = fc_contact_counts
)

# the same but with gender and occupation (confounded)
fc_model_lmer_all <- lmer(
  log1p(contacts) ~ (1|part_age) +
    (1|part_gender) +
    (1|part_occupation) +
    (1|part_id),
  data = fc_contact_counts
)

# now for Hong Kong

hk_model_lmer_none <- lmer(
  log1p(contacts) ~ (1|part_id),
  data = hk_contact_counts
)

hk_model_lmer_age <- lmer(
  log1p(contacts) ~ (1|part_age) +
    (1|part_id),
  data = hk_contact_counts
)

# the same but with gender (confounded)
hk_model_lmer_all <- lmer(
  log1p(contacts) ~ (1|part_age) +
    (1|part_gender) +
    (1|part_id),
  data = hk_contact_counts
)

# do variance partitioning on each, combining across studies
partitioning_none <- bind_rows(
  France = partition_variance_lmer(fc_model_lmer_none),
  `Hong Kong` = partition_variance_lmer(hk_model_lmer_none),
  .id = "Study"
)

partitioning_age <- bind_rows(
  France = partition_variance_lmer(fc_model_lmer_age),
  `Hong Kong` = partition_variance_lmer(hk_model_lmer_age),
  .id = "Study"
)

partitioning_all <- bind_rows(
  France = partition_variance_lmer(fc_model_lmer_all),
  `Hong Kong` = partition_variance_lmer(hk_model_lmer_all),
  .id = "Study"
)

barplot_none <- make_barplot(partitioning_none) +
  ggtitle(label = "no covariates")
barplot_age <- make_barplot(partitioning_age) +
  ggtitle(label = "age only")
barplot_all <- make_barplot(partitioning_all) +
  ggtitle(label = "additional predictors")


barplot_age +
  ggtitle(
    label = "Partitioning variance in daily (log1p) contact rates"
  )
ggsave("figures/varpart_age.png",
       bg = "white", scale = 0.75)


(barplot_none / barplot_age / barplot_all) +
  plot_annotation(title = "Partitioning variance in daily (log1p) contact rates")

ggsave("figures/varpart_all.png",
       bg = "white", scale = 1)


# sanity check this analysis by permuting the participant IDs

sims_list <- replicate(50,
                  permute_sim(),
                  simplify = FALSE)
sims_df <- do.call(bind_rows, sims_list)

partitions <- partitioning_age$partition
par(mfrow = c(length(partitions), 1))

for (this_partition in partitions) {
  perturbed_props <- sims_df[[this_partition]]
  estimated_prop <- partitioning_age |>
    filter(partition == this_partition) |>
    pull(proportion)
  
  xlim <- range(c(perturbed_props, estimated_prop))
  
  hist(perturbed_props,
       xlim = xlim,
       breaks = 100,
       xlab = "Proportion",
       main = this_partition)
  abline(v = estimated_prop)
}
