library(readxl)
library(tidyverse)
library(usethis)

# Read in from spreadsheet and filter only complete entries

pgs_metadata <- read_excel("//ess01/P471/data/durable/common/pgs_directory/supporting/prs_directory_metadata_safe.xlsx")%>%
  select(Sumstats_filename:runs)

use_data(pgs_metadata, overwrite=TRUE)


# Read in the ldpred2 pipeline Rscript to save as a package data object

ldp_raw = readLines("./data-raw/ldpred.pipl")

use_data(ldp_raw, overwrite=TRUE)
