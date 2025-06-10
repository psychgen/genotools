#' Get info about PGS from directory
#'
#' \code{get_pgs_log} copies the log files for the PGS creation into a specified directory
#' \code{get_pgs_snps} get lists of all SNPs available after clumping for any set of PGS
#' \code{pgs_nspns_thresh} get tables with number of SNPs at each thresholds for any set of PGS
#'
#' Detailed description...
#'
#' @param pgs_list a string or vector of strings from the Pheno_shortname column of available_pgs()
#' @param thresholds a vector of numeric values corresponding to the p-value thresholds you want
#' @param geno_data Which genotype data files do scores come from? Defaults to
#' "MoBaPsychGen_v1-ec-eur-batch-basic-qc"; also supports "98k-mobagenetics-ieu",
#' "98k-ec-eur-fin-batch-basic-qc", "hrv_rot1" for legacy reasons - these should
#' not be needed for new analyses
#' @param pgs_directory Where are the PGS stored locally? Defaults to
#' "//ess01/P471/data/durable/common/pgs_directory/pgs"
#' @param pgs_software With what software were PGS created? Defaults to "prsice2"
#' @param maf must be a string from pgs directory tree e.g. "0.05"; see all possibilities
#'  using \code{stringr::str_remove_all(unique(available_pgs()$maf), "maf")}
#' @param clump must be a string from pgs directory tree e.g. "500_1_0.25"; see all possibilities
#'  using \code{stringr::str_remove_all(unique(available_pgs()$clump), "clump")}
#' @rdname get_pgs_info
#' @export
#' @importFrom dplyr "%>%"

get_pgs_log <- function(pgs_list=NULL,
                        thresholds= c(5e-08, 5e-07, 5e-06, 5e-05, 5e-04, 0.001, 0.01, 0.005, 0.1, 0.5, 1),
                        geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                        pgs_directory="//ess01/P471/data/durable/common/pgs_directory/pgs",
                        pgs_software="prsice2",
                        maf="0.01",
                        clump="250_1_0.1"){

  dir.create(file.path("./pgs_log_files"), showWarnings = FALSE, recursive=TRUE)


  for(pgs in pgs_list){
  suppressMessages(file.copy(from = paste(pgs_directory,pgs_software,geno_data,
                         paste0("maf",maf), paste0("clump",clump),pgs,
                         paste0(pgs,".log"), sep="/"),
            to = paste0("./pgs_log_files/", pgs, "_", maf, "_", clump,".log" ) ))


  }

  message(paste0("\nRequested log files now available in ", getwd(),"/pgs_log_files"))
}

#'
#' @rdname get_pgs_info
#' @export
#'

get_pgs_snps <- function(pgs_list=NULL,
                         geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                         pgs_directory="//ess01/P471/data/durable/common/pgs_directory/pgs",
                         pgs_software="prsice2",
                         maf="0.01",
                         clump="250_1_0.10"){
  snplists <- list()
  for(pgs in pgs_list){
    tmp <- read.table(paste(pgs_directory,pgs_software,geno_data,
                                            paste0("maf",maf), paste0("clump",clump),pgs, paste0(pgs,".snp"), sep="/"),
                      header = T)
    snplists[[pgs]] <- tmp
  }

  message("\n NB, this function returns (as a list of dataframes) the list of SNPs available after clumping for
constructing each PGS, but this DOES NOT tell you whether or not a SNP is included at a given threshold. You need
to manually filter on the 'P' column for each threshold to get this information")
  return(snplists)
}


#'
#' @rdname get_pgs_info
#' @export
#'


get_pgs_nsnps <- function(pgs_list=NULL,
                          thresholds= c(5e-08, 5e-07, 5e-06, 5e-05, 5e-04, 0.001, 0.01, 0.005, 0.1, 0.5, 1),
                          geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                          pgs_directory="//ess01/P471/data/durable/common/pgs_directory/pgs",
                          pgs_software="prsice2",
                          maf="0.01",
                          clump="250_1_0.10"){

  snp_ns<- list()

  for(pgs in pgs_list){
  tmp <- read.table(paste(pgs_directory,pgs_software,geno_data
                          , paste0("maf",maf), paste0("clump",clump),pgs, paste0(pgs,".prsice"), sep="/"), header = T) %>%
    dplyr::select(Threshold, Num_SNP) %>%
    dplyr::mutate_all(as.numeric) %>%
    dplyr::filter(Threshold %in% thresholds )

  snp_ns[[pgs]] <- tmp

  }

return(snp_ns)
}
