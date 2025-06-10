#' Available PGS
#'
#' \code{available_pgs} returns a dataframe with details of already available
#' polygenic scores (PGS)
#'
#' Detailed description...
#'
#' @param geno_data Which genotype data files do scores come from? Defaults to
#' "MoBaPsychGen_v1-ec-eur-batch-basic-qc"; previous versions also  available for
#' continuity (but not recommended for new projects): "98k-mobagenetics-ieu",
#'  "98k-ec-eur-fin-batch-basic-qc","hrv_rot1"
#' @param pgs_directory Where are the PGS stored locally? Defaults to
#' "//ess01/P471/data/durable/common/pgs_directory/pgs"
#' @param pgs_software With what software were PGS created? Defaults to "prsice2"
#' @param pgs_meta What R object holds the PGS metadata? Defaults to
#' genotools::pgs_metadata
#' @param within_fam Are scores for within-family analyses such as trio-PGS? If so,
#' you may prefer to use a score based on only very well imputed (INFO>0.95) variants (see
#' https://www.medrxiv.org/content/10.1101/2024.10.01.24314703v1.full-text for
#' context on this decision). Defaults to FALSE.
#' @export
#' @importFrom dplyr "%>%"


available_pgs <- function(geno_data="MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                          pgs_directory="//ess01/P471/data/durable/common/pgs_directory/pgs",
                          pgs_software="ldpred2",
                          pgs_meta=genotools::pgs_metadata,
                          within_fam = FALSE){

  message("Scanning PGS directory...")

  if(pgs_software=="prsice2"){

    paths <- list.files(path = paste0(pgs_directory,"/",pgs_software ), include.dirs = TRUE, recursive = TRUE)
    mytable <- stringr::str_split_fixed(paths, pattern = "/", n = stringr::str_count(paths, "/") + 1) %>%
      as.data.frame() %>%
      dplyr::mutate_all(dplyr::na_if,"") %>%
      dplyr::select(1:3)

    colnames(mytable) <- c("maf","clump","Pheno_shortname")

        suppressMessages(mytable <- mytable %>%
                       tidyr::drop_na(Pheno_shortname) %>%
                       dplyr::distinct() %>%
                       dplyr::left_join(pgs_meta %>%
                                          dplyr::select(Pheno_shortname, Phenotype, Outcome_type,Source,`Ref (PMID)`, Sumstats_filename)))



  } else if(pgs_software=="ldpred2"){

    if(within_fam == FALSE){
      paths <- list.files(path = paste0(pgs_directory,"/",pgs_software, "/no_filter" ), include.dirs = TRUE, recursive = TRUE)
    } else if(within_fam == TRUE){
      paths <- list.files(path = paste0(pgs_directory,"/",pgs_software, "/filter_info_95"), include.dirs = TRUE, recursive = TRUE)
      }

    mytable <- stringr::str_split_fixed(paths, pattern = "/", n = ((length(stringr::str_count(paths, "/")) - sum(stringr::str_count(paths, "/"))) + 1)) %>%
      as.data.frame() %>%
      dplyr::filter(stringr::str_detect(V2, "_pred_auto.txt") ) %>%
      dplyr::select("Pheno_shortname" = V1)


    suppressMessages(mytable <- mytable %>%
                       tidyr::drop_na(Pheno_shortname) %>%
                       dplyr::distinct() %>%
                       dplyr::left_join(pgs_meta %>%
                                          dplyr::select(Pheno_shortname, Phenotype, Outcome_type,Source,`Ref (PMID)`, Sumstats_filename)))
  }

  return(mytable)

}

