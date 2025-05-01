# Download contact survey data


# Download all French Connection study data:

# Béraud G, Kazmercziak S, Beutels P, Levy-Bruhl D, Lenne X, Mielcarek N, et al.
# (2015) The French Connection: The First Large Population-Based Contact Survey
# in France Relevant for the Spread of Infectious Diseases. PLoS ONE 10(7):
# e0133203. https://doi.org/10.1371/journal.pone.0133203

# data on zenodo here: https://zenodo.org/records/3886590
fc_dir <- "data/french_connection"
dir.create(fc_dir, recursive = TRUE, showWarnings = FALSE)
library("socialmixr")
survey <- get_survey("https://doi.org/10.5281/zenodo.3886590")
saveRDS(survey, file.path(fc_dir, "survey.rds"))
