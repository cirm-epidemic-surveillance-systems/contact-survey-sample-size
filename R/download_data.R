# Download contact survey data


# Download all French Connection study data:

# Béraud G, Kazmercziak S, Beutels P, Levy-Bruhl D, Lenne X, Mielcarek N, et al.
# (2015) The French Connection: The First Large Population-Based Contact Survey
# in France Relevant for the Spread of Infectious Diseases. PLoS ONE 10(7):
# e0133203. https://doi.org/10.1371/journal.pone.0133203

# data on zenodo here: https://zenodo.org/records/3886590
fc_dir <- "data/french_connection"
dir.create(fc_dir, recursive = TRUE, showWarnings = FALSE)
fc_filepath <- file.path(fc_dir, "everything.zip")
download.file("https://zenodo.org/api/records/3886590/files-archive",
              fc_filepath)
unzip(fc_filepath,
      exdir = fc_dir)

# Download Hong Kong Survey data 
hk_dir <-"data/hongkong"
dir.create(hk_dir,recursive = TRUE, showWarnings = FALSE)
hk_filepath <-file.path(hk_dir,"HongKongData.csv")
download.file("https://royalsocietypublishing.org/action/downloadSupplement?doi=10.1098%2Frsif.2017.0838&file=rsif20170838supp2.csv",
              hk_filepath)
