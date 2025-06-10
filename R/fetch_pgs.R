#' Fetch PGS from directory
#'
#' \code{fetch_pgs} returns a dataframe with requested polygenic scores (PGS)
#'
#' Detailed description...
#'
#' @param pgs_list a string or vector of strings from the Pheno_shortname column of available_pgs()
#' @param thresholds (ignored unless using a clumping and thresholding approach to PGS creation)
#' a vector of numeric values corresponding to the p-value thresholds you want
#' or a single value between 0-1 to return all thresholds lower than this value by using the
#' threshold_range = TRUE option; defaults to c(5e-08, 5e-07, 5e-06, 5e-05, 5e-04, 0.001,
#' 0.01, 0.005, 0.1, 0.5, 1)
#' @param threshold_range (ignored unless using a clumping and thresholding approach to PGS creation)
#' use TRUE with a single value in thresholds to get a range; defaults to FALSE
#' @param geno_data Which genotype data files do scores come from? Defaults to
#' "MoBaPsychGen_v1-ec-eur-batch-basic-qc"; previous versions also available for continuity (but
#'  not recommended for new projects): "98k-mobagenetics-ieu","98k-ec-eur-fin-batch-basic-qc","hrv_rot1"
#' @param pgs_directory Where are the PGS stored locally? Defaults to
#' "//ess01/P471/data/durable/common/pgs_directory/pgs/ldpred2/no_filter"
#' @param pgs_software With what software were PGS created? Defaults to "ldpred2"
#' @param pgs_meta What R object holds the PGS metadata? Defaults to genotools::pgs_metadata
#' @param maf Used for retrieving scores made with PRSice2 based on the clumping and
#' thresholding method; must be a string from pgs directory tree e.g. "0.05"; see all possibilities
#'  using \code{stringr::str_remove_all(unique(available_pgs()$maf), "maf")}
#' @param clump Used for retrieving scores made with PRSice2 based on the clumping and
#' thresholding method; must be a string from pgs directory tree e.g. "500_1_0.25"; see all possibilities
#'  using \code{stringr::str_remove_all(unique(available_pgs()$clump), "clump")}
#' @param within_fam Are scores for within-family analyses such as trio-PGS? If so,
#' you may prefer to use a score based on only very well imputed (INFO>0.95) variants (see
#' https://www.medrxiv.org/content/10.1101/2024.10.01.24314703v1.full-text for
#' context on this decision). Defaults to FALSE.
#' @param zip_file is the file for the PGS saved as .gz, defualt is = TRUE
#' @export
#' @importFrom dplyr "%>%"




fetch_pgs <- function(pgs_list=NULL,
                      thresholds=c(5e-08, 5e-07, 5e-06, 5e-05, 5e-04, 0.001, 0.01, 0.005, 0.1, 0.5, 1),
                      threshold_range=FALSE,
                      geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                      pgs_directory="//ess01/P471/data/durable/common/pgs_directory/pgs/ldpred2/no_filter",
                      pgs_software="ldpred2",
                      pgs_meta=genotools::pgs_metadata,
                      maf="0.05",
                      clump="500_1_0.25",
                      within_fam = FALSE,
                      zip_file = TRUE){

  #SETUP
  #########################

  if(within_fam== TRUE) {

    pgs_directory <- stringr::str_replace_all(pgs_directory, "no_filter", "filter_info_95")
    message("within_fam = TRUE, PGS were constructed on summary statistics that have been filter based on info score > 0.95. If non-filtered scores are needed run with within_fam = FALSE")

  }


  message("\nRetrieving available PGS...")

  if(pgs_software=="ldpred2"){

    #Check pgs_list is provided
    if(is.null(pgs_list))
      stop("\nNo PGS requested. For details of available PGS, run available_pgs(). PGS must be requested using their Pheno_shortname.")
    #Check pgs requested in pgs_list are available and exclude any that are not
    if(length(setdiff(pgs_list,
                      list.files(paste0(pgs_directory,"/")))) != 0)
      warning(paste0("\n\nUnavailable PGS requested. The following PGS will not be retrieved:\n",
                     paste0(setdiff(pgs_list, list.files(paste0(pgs_directory,"/"))), collapse = ", "),
                     "\nPGS must be requested using their Pheno_shortname (see available_pgs())."))
    pgs_list<- pgs_list[!pgs_list %in% setdiff(pgs_list,
                                               list.files(paste0(pgs_directory,"/")))]
    #Throw an error if no available PGS remain in pgs_list
    if(length(pgs_list)==0)
      stop("\n\nNo available PGS remaining. For details of available PGS, run available_pgs(). If your PGS is listed, check you have used the correct Pheno_shortname.")

    #FETCHING SCORES
    ##########################
    all_pgs <- list()



    for(pgs in pgs_list){
      message(paste0("\nFetching ",pgs," scores... (This is #", match(pgs,pgs_list)," of ", length(pgs_list)," available phenotypes requested)"))

      if(zip_file == FALSE){
        suppressMessages(all_pgs[[pgs]] <- data.table::fread(paste0(pgs_directory,"/",pgs,"/",pgs,"_pred_auto.txt"), header = T) %>%
                           dplyr::select("FID"=family.ID,"IID"=sample.ID,"f_id"=paternal.ID,"m_id"=maternal.ID, !!pgs:= final_pred_auto))


      }else{

        suppressMessages(all_pgs[[pgs]] <- data.table::fread(paste0(pgs_directory,"/",pgs,"/",pgs,"_pred_auto.txt.gz"), header = T) %>%
                           dplyr::select("FID"=family.ID,"IID"=sample.ID,"f_id"=paternal.ID,"m_id"=maternal.ID, !!pgs:= final_pred_auto))


      }

        }


    suppressMessages(all_pgs_tbl <- all_pgs %>%
                       purrr::reduce(dplyr::left_join) %>%
                       dplyr::as_tibble() %>%
                       dplyr::rename_with(.cols=c(-FID,-IID,-f_id,-m_id),.fn = ~ paste0(.,"_pgs")))



  } else if (pgs_software=="prsice2"){

    pgs_directory <- stringr::str_replace_all(pgs_directory, "ldpred2/no_filter", "prsice2")

    #CHECKING INPUTS
    ##########################

    #Check MAF and clumping pars are acceptable
    if(!maf %in% stringr::str_remove(list.files(paste0(pgs_directory,"/",pgs_software,"/",geno_data,"/")), "maf"))
      stop(
        paste0("\nScores are not available at this MAF - you will need to create them from scratch. Available MAF thresholds: \n",
               stringr::str_remove_all(
                 paste0(
                   list.files(
                     paste0(pgs_directory,"/",pgs_software,"/",geno_data,"/")
                   )
                   ,collapse = ", ")
                 , "maf")
        )
      )
    #Check MAF and clumping pars are acceptable
    if(!clump %in% stringr::str_remove(list.files(paste0(pgs_directory,"/",pgs_software,"/",geno_data,"/maf",maf,"/")), "clump"))
      stop(
        paste0("\nScores are not available with these clumping parameters - you will need to create them from scratch. Available clumping parameter combinations: \n",
               stringr::str_remove_all(
                 paste0(
                   list.files(
                     paste0(pgs_directory,"/",pgs_software,"/",geno_data,"/maf",maf,"/")
                   )
                   ,collapse = ", ")
                 , "clump")
        )
      )

    #Check pgs_list is provided
    if(is.null(pgs_list))
      stop("\nNo PGS requested. For details of available PGS, run available_pgs(). PGS must be requested using their Pheno_shortname.")
    #Check pgs requested in pgs_list are available and exclude any that are not
    if(length(setdiff(pgs_list,
                      list.files(paste0(pgs_directory,"/",pgs_software,"/",geno_data,"/maf",maf,"/clump",clump,"/")))) != 0)
      warning(paste0("\n\nUnavailable PGS requested. The following PGS will not be retrieved:\n",
                     paste0(setdiff(pgs_list, list.files(paste0(pgs_directory,"/",pgs_software,"/",geno_data,"/maf",maf,"/clump",clump,"/"))), collapse = ", "),
                     "\nPGS must be requested using their Pheno_shortname (see available_pgs())."))
    pgs_list<- pgs_list[!pgs_list %in% setdiff(pgs_list,
                                               list.files(paste0(pgs_directory,"/",pgs_software,"/",geno_data,"/maf",maf,"/clump",clump,"/")))]
    #Throw an error if no available PGS remain in pgs_list
    if(length(pgs_list)==0)
      stop("\n\nNo available PGS remaining. For details of available PGS, run available_pgs(). If your PGS is listed, check you have used the correct Pheno_shortname.")




    #Check that the thresholds settings make sense
    if(threshold_range==TRUE & length(thresholds)>1){
      warning("\nYou have set threshold_range=TRUE and specified multiple values in thresholds.
When threshold_range=TRUE, a single 'thresholds1 value is required as the upper limit of the range (the lower limit is always 0), and the function will return scores at all available thresholds in the range.
Defaulting to threshold_range=FALSE and fetching scores at specified thresholds, if available." )
      threshold_range=FALSE}
    if(threshold_range==TRUE & length(thresholds)==1){
      message(paste0("\nThreshold range selected: ",paste0(range(0,thresholds),collapse=" - ") ))
    }


    ##########################

    #FETCHING SCORES
    ##########################
    all_pgs <- list()

    for(pgs in pgs_list){
      message(paste0("\nFetching ",pgs," scores... (This is #", match(pgs,pgs_list)," of ", length(pgs_list)," available phenotypes requested)"))
      suppressMessages(all_pgs[[pgs]] <- read.table(paste0(pgs_directory,"/",pgs_software,"/",geno_data,"/maf",maf,"/clump",clump,"/",pgs,"/",pgs,".all_score"), header = T) %>%
                         tidyr::gather(threshold, value, -IID,-FID) %>%
                         dplyr::mutate(threshold = stringr::str_remove(threshold,"Pt_"),
                                       threshold = stringr::str_replace(threshold,"e.","e-"),
                                       threshold = as.numeric(threshold)))

      if(threshold_range==TRUE){
        suppressMessages(all_pgs[[pgs]] <- all_pgs[[pgs]] %>%
                           dplyr::filter(threshold<thresholds) %>%
                           tidyr::spread(threshold,value) %>%
                           dplyr::rename_at(dplyr::vars(-IID, -FID) ,function(x){paste0(pgs,"_p<", x)}))
      }else{
        suppressMessages(all_pgs[[pgs]] <- all_pgs[[pgs]] %>%
                           dplyr::filter(threshold%in%thresholds) %>%
                           tidyr::spread(threshold,value) %>%
                           dplyr::rename_at(dplyr::vars(-IID, -FID) ,function(x){paste0(pgs,"_p<", x)}))
      }

      if(dim(all_pgs[[pgs]])[2]<3 )
        stop("\n\nFor some reason, none of the thresholds you specified are available. Check your input or use threshold_range==TRUE and thresholds=1 to get all available thresholds.")

    }



    suppressMessages(all_pgs_tbl <- all_pgs %>%
                       purrr::reduce(dplyr::left_join) %>%
                       dplyr::as_tibble())

    if(threshold_range==TRUE){message("\nFetch complete: here are all available scores from your request, within the requested threshold range.
To get logfiles from the creation of any set of PGS, run get_pgs_log().
To get lists including details of all SNPs available after clumping for any set of PGS, run get_pgs_snps().
To get tables with the number of SNPs included at each threshold for any set of PGS, run get_pgs_nsnps().
\nFor merging with phenotypic data, regressing on batch and PCs and/or retreiving these as covariates for analyses,
and retrieving relevant info to exclude individuals based on relatedness, use process_pgs(). ")
    }else{message("\nFetch complete: here are all available scores from your request, at all available requested thresholds.
To get logfiles from the creation of any set of PGS, run get_pgs_log().
To get lists including details of all SNPs available after clumping for any set of PGS, run get_pgs_snps().
To get tables with the number of SNPs included at each threshold for any set of PGS, run get_pgs_nsnps().
\nFor merging with phenotypic data, regressing on batch and PCs and/or retreiving these as covariates for analyses,
and retrieving relevant info to exclude individuals based on relatedness, use process_pgs().")}

  } else { stop("Input for option 'pgs_software' not recognised; currently only 'prsice2' and 'ldpred2' are valid inputs.")}

  return(all_pgs_tbl)

}

