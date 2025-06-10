#' Process PGS fetched from directory
#'
#' \code{fetch_pgs} takes raw PGS and joins to phenotypic data IDs, whilst also
#' regressing on key covariates and/or returning these covariates for inclusion in
#' analyses.
#'
#' Detailed description...
#'
#' @param fetched_pgs_df input should be data.frame results of fetch_pgs
#' @param regress_pgs do you want PGS regressed on batch and 20PCs? Defaults to TRUE
#' @param geno_data Which genotype data files do scores come from? Defaults to
#' "MoBaPsychGen_v1-ec-eur-batch-basic-qc"; also supports "98k-mobagenetics-ieu",
#' "98k-ec-eur-fin-batch-basic-qc", "hrv_rot1" for legacy reasons - these should
#' not be needed for new analyses
#' @param covs_dir Where is the covariate file for the geno_data? Defaults to
#' "//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1" for p471
#' @param return_covs Do you want batch and PCS returned? Defaults to TRUE
#' @param analysis_type Not yet implemented (for exclusion flag list vars)
#' @export
#' @importFrom dplyr "%>%"



process_pgs <- function(fetched_pgs_df,
                        regress_pgs = TRUE,
                        geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                        covs_dir = "//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1",
                        return_covs=TRUE,
                        analysis_type = NULL){

  ##Get correct cov file and adjust names to be consistent

  if(geno_data == "MoBaPsychGen_v1-ec-eur-batch-basic-qc"){

    covs <- readr::read_tsv(paste0(covs_dir,"/",geno_data,"-cov.txt"), col_types = readr::cols(.default = "c")) %>%
      dplyr::mutate_at(dplyr::vars(dplyr::matches("PC")), as.numeric) %>%
      dplyr::select(id_moba = dplyr::matches("ID_",ignore.case = FALSE), dplyr::everything()) %>%
      dplyr::mutate(preg_id_BARN_NR = ifelse(Role =="Child", id_moba,NA),
                    f_id = ifelse(Role == "Father",id_moba,NA),
                    m_id = ifelse(Role == "Mother",id_moba,NA)) %>%
      tidyr::separate(preg_id_BARN_NR, into=c("preg_id", "BARN_NR"), sep="_")

  } else if(geno_data == "hrv_rot1") {
    covs <- read.table("//ess01/P471/data/durable/data/Linkage files/core_IDs&covars_hrv_njl_v2.txt", header = T) %>%
      dplyr::rename(genotyping_batch = BATCH) %>%
      dplyr::select(preg_id = dplyr::matches("PREG_ID_",ignore.case = FALSE), dplyr::everything())
  } else if(geno_data =="98k-ec-eur-fin-batch-basic-qc") {
    covs <- read.table("//ess01/P471/data/durable/data/Linkage files/98k-ec-eur-fin-batch-basic-qc-cov.txt", header = T, sep="\t") %>%
      dplyr::select(id_moba = dplyr::matches("ID_",ignore.case = FALSE), dplyr::everything()) %>%
      dplyr::mutate(preg_id_BARN_NR = ifelse(Role =="Child", id_moba,NA),
                    f_id = ifelse(Role == "Father",id_moba,NA),
                    m_id = ifelse(Role == "Mother",id_moba,NA)) %>%
      tidyr::separate(preg_id_BARN_NR, into=c("preg_id", "BARN_NR"), sep="_")
  } else if(geno_data =="98k-mobagenetics-ieu") {
    covs <- read.table("//ess01/P471/data/durable/projects/moba_interim_release_post_imp_qc/FINAL QC VERSION/merged_cov_file_ljh.txt", header = T, sep="\t") %>%
      dplyr::select(id_moba = dplyr::matches("ID_",ignore.case = FALSE), dplyr::everything()) %>%
      dplyr::mutate(preg_id_BARN_NR = ifelse(Role =="Child", id_moba,NA),
                    f_id = ifelse(Role == "Father",id_moba,NA),
                    m_id = ifelse(Role == "Mother",id_moba,NA)) %>%
      tidyr::separate(preg_id_BARN_NR, into=c("preg_id", "BARN_NR"), sep="_")
  }else if(geno_data =="release2-ec-eur-batch-basic-qc") {
    covs <- read.table("//ess01/P471/data/durable/data/genetic/Release2/release2-ec-eur-batch-basic-qc-cov.txt", header = T, sep="\t") %>%
      dplyr::select(id_moba = dplyr::matches("ID_",ignore.case = FALSE), dplyr::everything()) %>%
      dplyr::mutate(preg_id_BARN_NR = ifelse(Role =="Child", id_moba,NA),
                    f_id = ifelse(Role == "Father",id_moba,NA),
                    m_id = ifelse(Role == "Mother",id_moba,NA)) %>%
      tidyr::separate(preg_id_BARN_NR, into=c("preg_id", "BARN_NR"), sep="_")
  }else {
    stop("geno_data input not recognised - check available_pgs() for possibilities")
  }

  if(any(stringr::str_detect(names(fetched_pgs_df), "m_id"))){

    suppressMessages(suppressWarnings(processed_pgs <- fetched_pgs_df %>%
                                        dplyr::select(-m_id,-f_id)%>%
                                        dplyr::left_join(covs)))
  } else {

    suppressMessages(suppressWarnings(processed_pgs <- fetched_pgs_df %>%
                                        dplyr::left_join(covs)))
  }

  ##Regress PGS on covariates if requested, otherwise just select vars to return
  if(regress_pgs == TRUE){
    if(geno_data=="MoBaPsychGen_v1-ec-eur-batch-basic-qc") {
      model<- function(y){
        m<- lm( y ~  processed_pgs$PC1 + processed_pgs$PC2 + processed_pgs$PC3 + processed_pgs$PC4 + processed_pgs$PC5 +
                  processed_pgs$PC6 + processed_pgs$PC7 + processed_pgs$PC8 + processed_pgs$PC9 + processed_pgs$PC10 +
                  processed_pgs$PC11 + processed_pgs$PC12 + processed_pgs$PC13 + processed_pgs$PC14 + processed_pgs$PC15 +
                  processed_pgs$PC16 + processed_pgs$PC17 + processed_pgs$PC18 + processed_pgs$PC19 + processed_pgs$PC20 +
                  processed_pgs$genotyping_batch + processed_pgs$imputation_batch,  na.action = na.exclude)
        return( rstandard(m))
      }

      message("Regressing out genotyping batch, imputation batch, and 20PCs from all supplied scores...")

    }else if(geno_data == "hrv_rot1") {
      model<- function(y){
        m<- lm( y ~  processed_pgs$PC1 + processed_pgs$PC2 + processed_pgs$PC3 + processed_pgs$PC4 + processed_pgs$PC5 +
                  processed_pgs$PC6 + processed_pgs$PC7 + processed_pgs$PC8 + processed_pgs$PC9 + processed_pgs$PC10 +
                  processed_pgs$genotyping_batch,  na.action = na.exclude)
        return( rstandard(m))
      }

    }else if(geno_data %in% c("98k-ec-eur-fin-batch-basic-qc")){
      model<- function(y){
        m<- lm( y ~  processed_pgs$PC1 + processed_pgs$PC2 + processed_pgs$PC3 + processed_pgs$PC4 + processed_pgs$PC5 +
                  processed_pgs$PC6 + processed_pgs$PC7 + processed_pgs$PC8 + processed_pgs$PC9 + processed_pgs$PC10 +
                  processed_pgs$genotyping_batch + processed_pgs$imputation_batch,  na.action = na.exclude)
        return( rstandard(m))
      }
    }else if(geno_data %in% c("release2-ec-eur-batch-basic-qc")){
      model<- function(y){
        m<- lm( y ~  processed_pgs$PC1 + processed_pgs$PC2 + processed_pgs$PC3 + processed_pgs$PC4 + processed_pgs$PC5 +
                  processed_pgs$PC6 + processed_pgs$PC7 + processed_pgs$PC8 + processed_pgs$PC9 + processed_pgs$PC10 +
                  processed_pgs$genotyping_batch ,  na.action = na.exclude)
        return( rstandard(m))
      }
    }else if(geno_data == "98k-mobagenetics-ieu"){
      model<- function(y){
        m<- lm( y ~  processed_pgs$PC1 + processed_pgs$PC2 + processed_pgs$PC3 + processed_pgs$PC4 + processed_pgs$PC5 +
                  processed_pgs$PC6 + processed_pgs$PC7 + processed_pgs$PC8 + processed_pgs$PC9 + processed_pgs$PC10 +
                  processed_pgs$PC11 + processed_pgs$PC12 + processed_pgs$PC13 + processed_pgs$PC14 + processed_pgs$PC15 +
                  processed_pgs$PC16 + processed_pgs$PC17 + processed_pgs$PC18 + processed_pgs$PC19 + processed_pgs$PC20 +
                  processed_pgs$genotyping_batch ,  na.action = na.exclude)
        return( rstandard(m))
      }
    }


    if(return_covs==FALSE){
      processed_pgs <- processed_pgs%>%
        dplyr::mutate(dplyr::across(c(tidyr::matches("_p<|_pgs")), list(res = ~model(.)))) %>%
        dplyr::select(IID,FID,tidyr::matches("preg_id"),BARN_NR,tidyr::matches("f_id"),tidyr::matches("m_id"),Role,
                      tidyr::matches("_p<|_pgs"))

    }else{
      if(geno_data%in%c("hrv_rot1","98k-mobagenetics-ieu","release2-ec-eur-batch-basic-qc") ){
        processed_pgs <- processed_pgs%>%
          dplyr::mutate(dplyr::across(c(tidyr::matches("_p<|_pgs")), list(res = ~model(.)))) %>%
          dplyr::select(IID,FID,tidyr::matches("preg_id"),BARN_NR,tidyr::matches("f_id"),tidyr::matches("m_id"),Role,
                        tidyr::matches("_p<|_pgs"), tidyr::matches("PC"),genotyping_batch)
      }else {
        processed_pgs <- processed_pgs%>%
          dplyr::mutate(dplyr::across(c(tidyr::matches("_p<|_pgs")), list(res = ~model(.)))) %>%
          dplyr::select(IID,FID,tidyr::matches("preg_id"),BARN_NR,tidyr::matches("f_id"),tidyr::matches("m_id"),Role,
                        tidyr::matches("_p<|_pgs"), tidyr::matches("PC"),genotyping_batch,imputation_batch)
      }
    }

    message("\nProcessing complete.
PGS variables with _res suffix are standardised residuals from the regression model.
You also now have ID variables (preg_id, F_ID, M_ID) for linkage to phenotypic data.")
  }else{
    if(return_covs==FALSE){
      processed_pgs <- processed_pgs%>%
        dplyr::select(IID,FID,tidyr::matches("preg_id"),BARN_NR,tidyr::matches("f_id"),tidyr::matches("m_id"),Role,
                      tidyr::matches("_p<|_pgs"))
      warning("\nYou have requested neither regression of PGS on covariates nor return of covariates;
it is highly recommended that you account for potential effects of genotyping batch
and population structure in your analyses by adjusting for these covariates.")
    }else{
      processed_pgs <- processed_pgs%>%
        dplyr::select(IID,FID,tidyr::matches("preg_id"),BARN_NR,tidyr::matches("f_id"),tidyr::matches("m_id"),Role,
                      tidyr::matches("_p<|_pgs"), tidyr::matches("PC"),genotyping_batch)
    }

    message("\nProcessing complete.
You also now have ID variables (preg_id, f_id, m_id) for linkage to phenotypic data.")

  }



  ## Section for adding in relatedness flags depending on analysis_type input

  ##

  return(processed_pgs)
}

