#' Get an unrelated sub-sample of MoBa data based on your analysis type
#'
#' \code{unrelate} takes a MoBa dataset with IDs (genetic, phenotypic, or both),
#' and - based on the kind of relatedness you declare as allowed in your analysis -
#' removes individuals to limit additional relatedness while prioritising
#' genotype and phenotype data availability
#'
#' Detailed description...
#'
#' @param input_df a data.frame with at least i) preg_id + BARN_NR;
#' ii) FID and IID; or (if analysis only includes parents) iii) m_id and/or f_id
#' @param unrelated a string describing who must be unrelated to one another in
#' the output; can be "trios" (for unrelated trios), "nuclear" (for unrelated
#' nucleary families - i.e., trios + other children) "children" (for unrelated
#' children, excluding even all-but-one child per nuclear family), "sibs" (for
#' unrelated children outside nuclear family sib pairs), "mothers", "fathers",
#' "parents","mcduos" (for unrelated sets of mother-child duos), or "fcduos"
#' (for unrelated sets of father-child duos)
#' @param complete should only complete sets of the specified unrelated
#' groupings be returned? Defaults to FALSE, meaning that incomplete
#' groupings can also be returned (i.e., mother child duos when unrelated = "trios" -
#' as long as they are not related to an included complete trio)
#' @param phenotypes column names from 'input_df' containing phenotypes you want to
#' prioritise data availability for when excluding related individuals
#' @param geno_data Which genotype data files do scores come from? Defaults to
#' "MoBaPsychGen_v1-ec-eur-batch-basic-qc"; legacy genotype data files not supported
#' by this function (cf. the _pgs functions)
#' @param covs_dir Where is the covariate file for the geno_data? Defaults to
#' "//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1" for p471
#' @param keep_flag Keep related individuals and add a flag variable? Defaults
#' to FALSE - i.e., remove them from the final dataset
#' @export
#' @importFrom dplyr "%>%"



unrelate <- function(input_df=NULL,
                     unrelated = "trios",
                     complete=FALSE,
                     phenotypes=NULL,
                     geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                     covs_dir = "//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1",
                     keep_flag=FALSE){

  #Check inputs
  #Does input_df have correct ID cols?

  #Does unrelated string correspond to available options?


  ###### Read in necessary files

  message("\nReading in files for relatedness checking...")

  link <- readr::read_tsv(paste0(covs_dir,"/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt")) %>%
    dplyr::select(ID_2306, SENTRIXID, FID, IID, Role, SEX,
                  YOB, genotyping_chip,PC1, PC2, PC3, PC4,
                  PC5, PC6, PC7, PC8, PC9, PC10 )

  fam <- readr::read_delim(paste0(covs_dir,"/MoBaPsychGen_v1-ec-eur-batch-basic-qc.fam"),col_names=F) %>%
    `colnames<-`(c("FID", "IID", "PID", "MID", "sex", "Phenotype" )) %>%
    dplyr::select(FID, IID, PID, MID)

  ###### Join link and fam for spine

  spine <- link %>%
    dplyr::select(ID_2306,IID, FID, Role, SEX)%>%
    dplyr::left_join(fam) %>%
    dplyr::mutate(preg_id_BARN_NR= ifelse(Role=="Child", ID_2306, NA),
                  f_id= ifelse(Role=="Father", ID_2306, NA),
                  m_id= ifelse(Role=="Mother", ID_2306, NA)) %>%
    tidyr::separate(preg_id_BARN_NR, into=c("preg_id","BARN_NR"), sep="_")

  ##### Join input_df if exists (reshape to one row per individual if needed)

  if(!is.null(input_df)){
    dat <- spine %>%
      dplyr::right_join(input_df)
  }else{
    dat <- spine
  }

  ##### So, from here, the goal is to include as many of the unique IIDs in dat
  ##### as possible in the final analytic dataset, adhering to the restrictions on
  ##### relatedness and priorities for inclusion set by the user

  #### Read in files with info on relatedness

  genome <- readr::read_delim(paste0(covs_dir,"/MoBaPsychGen_v1-ec-eur-batch-basic-qc-ibd.genome"),col_names=T)

  kin<-  readr::read_delim(paste0(covs_dir,"/MoBaPsychGen_v1-ec-eur-batch-basic-qc-rel.kin"),col_names=T)


  # Summarise N per FID

  dat_fid <- dat %>%
    dplyr::group_by(FID) %>%
    dplyr::summarise(n=dplyr::n())

  ##### UNREALTED TRIOS

  if (unrelated=="trios"){

    message("\nIdentifying unrelated trios...")

    ###### Identify "easy" trios (3xFID in fam;  Role == "Mother", "Father,"Child"; rxp in genome==PO*2)
    trios <- dat %>%
      dplyr::filter(FID %in% c(
        dat_fid %>%
          dplyr::filter(n==3) %>%
          .$FID )) %>%
      dplyr::group_by(FID) %>%
      dplyr::filter(dplyr::n_distinct(Role)==3)%>%
      dplyr::left_join(genome %>%
                         dplyr::select("FID"=FID1,RT)) %>%
      dplyr::distinct() %>%
      dplyr::filter(RT=="PO") %>%
      dplyr::select(FID,IID,Role)

    ###### Identify "harder" trios (>3xFID in fam;  at least 3 distinct Roles ("Mother", "Father,"Child"); rxp in genome==PO*2)
    nucs <- dat %>%
      dplyr::filter(FID %in% c(
        dat_fid %>%
          dplyr::filter(n>3) %>%
          .$FID )) %>%
      dplyr::group_by(FID) %>%
      dplyr::filter(dplyr::n_distinct(Role)==3)

    nuc_rxps <- genome %>%
      dplyr::filter(FID1==FID2,
                    FID1 %in% nucs$FID) %>%  #All unique relationships in these families - for now we only want parent-offspring relationships
      dplyr::filter(RT=="PO") %>%
      dplyr::select(IID1,IID2) %>%
      tidyr::pivot_longer(tidyr::everything(), values_to="IID") %>%
      dplyr::select(IID) %>%
      dplyr::distinct()

    nucs <- nucs %>%
      dplyr::right_join(nuc_rxps) %>%
      dplyr::filter(!(Role=="Child"&!PID%in%IID)) %>%
      dplyr::filter(!(Role=="Child"&!MID%in%IID)) %>%
      dplyr::filter(!(Role=="Father"&!IID%in%PID)) %>%
      dplyr::filter(!(Role=="Mother"&!IID%in%MID) )   #Prioritise complete trios in extended familes by dropping parents whose IDs are not listed for a child in the family

    if(is.null(phenotypes)){
      phenotypes=c("IID")
    }

    set.seed(31409)
    slct_nuc_children <- nucs %>%
      dplyr::filter(Role=="Child") %>%
      dplyr::rowwise() %>%
      dplyr::mutate(na=sum(is.na(dplyr::c_across(dplyr::all_of(phenotypes))))) %>%
      dplyr::group_by(FID) %>%
      dplyr::filter(na==min(na)) %>%
      dplyr::slice_sample(n=1)

    trios_from_nucs <- nucs %>%
      dplyr::filter(IID %in% slct_nuc_children$IID) %>%
      dplyr::bind_rows(nucs %>%
                         dplyr::filter(IID %in% slct_nuc_children$PID)) %>%
      dplyr::bind_rows(nucs %>%
                         dplyr::filter(IID %in% slct_nuc_children$MID)) %>%
      dplyr::select(FID,IID,Role)

    all_complete_trios <- trios %>%
      dplyr::bind_rows(trios_from_nucs) %>%
      dplyr::distinct()

    return_dat <- all_complete_trios
  }

  if( keep_flag==TRUE){
    return_dat <- spine %>%
      dplyr::left_join(return_dat %>%
                         dplyr::mutate(keep=TRUE))
    message("\nYou set keep_flag=TRUE, so your input dataset is returned in full with individuals flagged
(value TRUE on variable 'keep') for inclusion in analysis")
  } else {
    return_dat <- return_dat %>%
      dplyr::left_join(spine)
  }
  return(return_dat)
}

