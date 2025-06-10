library(tidyverse)
#read unprocessed additions
files = list.files("//ess01/P471/data/durable/common/pgs_directory/add_to_dir") 
csvs= paste0("//ess01/P471/data/durable/common/pgs_directory/add_to_dir/",
             files[str_detect(files,".csv")])
all = map(csvs,readr::read_csv, col_types=cols(.default = "c")) %>% 
  purrr::reduce(bind_rows)
#combine in a single file
write_csv(all, file="//ess01/P471/data/durable/common/pgs_directory/add_to_dir/combined.csv")
#add to master file
master=readxl::read_xlsx("//ess01/P471/data/durable/common/pgs_directory/supporting/prs_directory_metadata_safe.xlsx")
new=master %>% bind_rows(all %>% 
                           rename("Sumstats_filename"=filename,`Ref (PMID)`=Pmid)) %>% 
  distinct()
#save and mark for QC
openxlsx::write.xlsx(new, file= "//ess01/P471/data/durable/common/pgs_directory/supporting/prs_directory_metadata_forQC.xlsx")
###
#move all to processed
file.copy(from=csvs,to="//ess01/P471/data/durable/common/pgs_directory/add_to_dir/processed")
file.remove(csvs)
#clean up
file.remove("//ess01/P471/data/durable/common/pgs_directory/add_to_dir/combined.csv")
