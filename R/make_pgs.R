#' Make and run ob script for polygenic score creation
#'
#' \code{make_pgs} makes and submits a submission script for
#' running on Colossus to make polygenic scores
#'
#' Detailed description...
#'
#' @param user username for cluster
#' @param host hostname for cluster - defaults to p471 submit node
#' @param software which software to used to make the PGS? Currently "ldpred2" (default) or "prsice2"
#' @param shortname shortname for the PGS - must be distinct from others in available_pgs()$Pheno_shortname
#' @param account the colossus/SLURM account from which CPU useage will be deducted for making the scores;
#' in p471, this can be "p471" (the default) or "p471_tsd" - use 'qsumm' command in bash to see availability
#' @param cpu_time how long is required for the job on Colossus (defaults to 05:00:00)
#' @param memory how much memory per core is required for the job on Colossus (defaults to 8GB)
#' @param cores how many cores should be used? Defaults to 32
#' @param lib (optional) add a path to a common R library with required packages installed
#' @param sumstats_dir directory for sumstats, suggest leaving at default if working in p471 - NB, this filepath
#' will be used by the HPC on linux, so should be specified accordingly ("ess" not "ess01")
#' @param outputs_dir directory for outputs; NB, this filepath
#' will be used by the HPC on linux, so should be specified accordingly ("ess" not "ess01")
#' @param sumstats_filename sumstats filename, including extension
#' @param genotype_dir location of genotype data on cluster; defaults for mobapsychgen_v1
#' @param genotype_data name of genotype data on cluster, formatted as an RDS file containing the subset
#' of the MoBa genetic data that overlaps with HapMap3 SNPs...; defaults set for mobapsychgen_v1
#' @param ldref_dir location of the HAPMAP3+ and ld blocks needed by ldpred2
#' @param case_control whether the GWAS is case control - will try to pull from existing metadata if left blank; defaults to TRUE
#' @param A1 effect allele column name - will try to pull from existing metadata if left blank
#' @param A2 other allele column - will try to pull from existing metadata if left blank
#' @param stat stat column - will try to pull from existing metadata if left blank
#' @param stat_type whether the stat is beta or OR- will try to pull from existing metadata if left blank
#' @param pvalue pval column - will try to pull from existing metadata if left blank
#' @param SNP snp id column - will try to pull from existing metadata if left blank
#' @param CHR CHR column - will try to pull from existing metadata if left blank
#' @param BP bp column - will try to pull from existing metadata if left blank
#' @param SE standard error column - will try to pull from existing metadata if left blank
#' @param N either, if trait is continuous a) a column name with SNP-wise Ns or b) a value
#' representing the overall N, or - if binary - c) a comma-separated pair of column names for SNP-wise
#' N-cases and N-controls or d) a vector of values representing the N-cases and N-controls
#' @param info (optional) info score column name for filtering
#' @param MAF (optional) maf column name for filtering - NB can also provide EAF here
#' @param maf_threshold minor allele frequency threshold - defaults to 0.01
#' @param source Where do the summary stats come from? (ideally a url) - either this or pmid
#' must be non-missing - will try to pull from existing metadata if left blank
#' @param pmid What publication do the summary stats relate to? - either this or source must be
#' non-missing - will try to pull from existing metadata if left blank
#' @param pgs_meta What R object holds the PGS metadata? Defaults to genotools::pgs_metadata
#' @param sd_maf calculate sd from an allele freq in summary stats - default is FALSE
#' @param session (optional) provide an active SSH session to avoid repeated login prompts if
#' mapping or looping this function
#' @param submit_ssh true/false to submit the job to the cluster. Default  is true. If false dir, .pipl file, and .sh file are created and can be edited and submitted manually.

#' @export
#' @importFrom dplyr "%>%"


make_pgs <- function(user,
                     host="p471-hpc-01.tsd.usit.no",
                     software = "ldpred2",
                     shortname,
                     account="p471",
                     cpu_time="05:00:00",
                     memory="8G",
                     cores=32,
                     lib="//ess/p471/data/durable/common/rlibs/4-2-1",
                     sumstats_dir="//ess/p471/data/no-backup/gwas_sumstats/",
                     outputs_dir= NA,
                     sumstats_filename = NA,
                     genotype_dir="//ess/p471/data/durable/data/genetic/MoBaPsychGen_v1/",
                     genotype_data="genoHapMap3plus_N200k.rds",
                     ldref_dir="//ess/p471/data/durable/data/genetic/ldref/",
                     case_control = TRUE,
                     A1 = NA,
                     A2 = NA,
                     stat = NA,
                     stat_type = NA,
                     pvalue = NA,
                     SNP = NA,
                     CHR= NA,
                     BP = NA,
                     SE = NA,
                     N = NA,
                     info  = NA,
                     MAF = NA,
                     maf_threshold = 0.01,
                     source = NA,
                     pmid =NA,
                     sd_maf = FALSE,
                     pgs_meta=genotools::pgs_metadata,
                     session=NULL,
                     submit_ssh=TRUE
){

  ## First check if a validated version of the score exists, and ask if the user wants to stop and use it if it does

  if(length(list.files("//ess01/p471/data/durable/common/pgs_directory/pgs/ldpred2"))>0){
    if(any(stringr::str_detect(list.files("//ess01/p471/data/durable/common/pgs_directory/pgs/ldpred2/no_filter"), shortname))){
      cont = readline(prompt="A validated version of this score seems already to exist in /data/durable/pgs_directory/ldpred2/no_filter/. Do you wish to continue, or abort and use this score instead? [continue/abort]")
      stopinot(cont=="continue")
    }
  }

  ## Try to extract necessary metadata if not provided, throw errors if metadata does not exist

  if(anyNA(c(sumstats_filename,A1,A2,stat,stat_type, case_control,pvalue,SNP,CHR,BP,SE,N))){
    tmp = genotools::pgs_metadata %>% dplyr::filter(Pheno_shortname==shortname)
    if(nrow(tmp)==0){
      stop(paste0("Some essential information (sumstats_filename, A1, A2, stat,stat_type,case_control, pvalue, SNP, CHR, BP, SE, or N) was not provided.
                  \nGenotools looked for metadata in genotools::pgs_metadata, but could not find an entry for ",shortname,".
                  If you beleive the relevant metadata exists, check that your input ",shortname," matches the corresponding Pheno_shortname value in genotools::pgs_metadata. If not, re-run the function providing the missing information as an argument (e.g., \"A1=\"effect_allele\"\""))
    }
    if(is.na(sumstats_filename)){
      sumstats_filename = tmp$Sumstats_filename
    }
    if(is.na(A1)){
      A1 = tmp$A1
    }
    if(is.na(A2)){
      A2 = tmp$A2
    }
    if(is.na(stat)){
      stat = tmp$stat
    }
    if(is.na(stat_type)){
      stat_type = tmp$stat_type
    }
    if(is.na(case_control)){
      case_control = tmp$case_control
    }
    if(is.na(pvalue)){
      pvalue = tmp$pvalue
    }
    if(is.na(SNP)){
      SNP = tmp$SNP
    }
    if(is.na(CHR)){
      CHR = tmp$CHR
    }
    if(is.na(BP)){
      BP = tmp$BP
    }
    if(is.na(SE)){
      SE = tmp$SE
    }
    if(is.na(N)){
      N = tmp$N
    }
  }
  if(anyNA(c(sumstats_filename,A1,A2,stat,stat_type,case_control,pvalue,SNP,CHR,BP,SE,N))){


    ess = c("sumstats_filename","A1","A2","stat","stat_type","case_control", "pvalue","SNP","CHR","BP","SE","N")

    stop(paste0("Some essential information (", ess[which(is.na(c(sumstats_filename,A1,A2,stat, stat_type, case_control,pvalue,SNP,CHR,BP,SE,N)))],"), was not provided and could not be found in the genotools::pgs_metadata entry for ",shortname,".
              Please re-run the function providing the missing information as an argument (e.g., \"A1=\"effect_allele\"\"" ))
  }

  ## Ensure case is correct for stat_type

  stat_type = toupper(stat_type)

  ## Check that either source of PMID are non-missing

  if(all(is.na(c(source,pmid)))){
    tmp = genotools::pgs_metadata %>% dplyr::filter(Pheno_shortname==shortname)
    if(nrow(tmp)==0){
      stop(paste0("Neither source nor pmid was were provided.
                  \ngenotools looked for metadata in genotools::pgs_metadata, but could not find an entry for ",shortname,".
                  If you believe the relevant metadata exists, check that your input ",shortname," matches the corresponding Pheno_shortname value in genotools::pgs_metadata. If not, re-run the function providing the missing information as an argument (e.g., \"A1=\"effect_allele\"\""))
    }
    if(is.na(source)){
      source = tmp$source
    }
    if(is.na(pmid)){
      pmid = tmp$pmid
    }

  }
  source = paste(source, pmid, sep= "; PMID: ")

  if(is.na(MAF)) {
    MAF = "NA"
  }

  if(is.na(info)) {
    info = "NA"
  }


  ## LDPRED2 specific section

  if(software=="ldpred2"){

    ## Make a user/job specific output directory

    outputs_subdir = paste(outputs_dir,user,shortname, sep="/")
    dir.create(paste(stringr::str_replace(outputs_dir, "/ess/", "/ess01/"),user,shortname, sep="/"), showWarnings = F, recursive = T)


    ## All essential info is present - make the R script that runs ldpred2 on the cluster

    fileConn<-file(paste(stringr::str_replace(outputs_dir, "/ess/", "/ess01/"),
                         user,
                         shortname,
                         paste0(shortname,"-ldpred_pipl.R"),
                         sep="/"),
                   "wb")

    ldp_pop = genotools::ldp_raw %>%
      stringr::str_replace_all( "opt\\$lib", lib) |>
      stringr::str_replace_all( "opt\\$out", outputs_subdir) |>
      stringr::str_replace_all( "opt\\$misc", ldref_dir) |>
      stringr::str_replace_all( "opt\\$cores", as.character(cores)) |>
      stringr::str_replace_all( "opt\\$geno", paste0(genotype_dir,genotype_data)) |>
      stringr::str_replace_all( "opt\\$Sdir", sumstats_dir) |>
      stringr::str_replace_all( "opt\\$shortname", shortname) |>
      stringr::str_replace_all( "opt\\$stattype", stat_type) |>
      stringr::str_replace_all( "opt\\$samplesize", as.character(N)) |>
      stringr::str_replace_all( "opt\\$type", as.character(case_control)) |>
      stringr::str_replace_all( "opt\\$sumstats", sumstats_filename) |>
      stringr::str_replace_all( "opt\\$chr", CHR)|>
      stringr::str_replace_all( "opt\\$pos", BP)|>
      stringr::str_replace_all( "opt\\$A1", A1)|>
      stringr::str_replace_all( "opt\\$A2", A2)|>
      stringr::str_replace_all( "opt\\$se", SE)|>
      stringr::str_replace_all( "opt\\$stat", stat)|>
      stringr::str_replace_all( "opt\\$MAFname", MAF)|>
      stringr::str_replace_all( "opt\\$sd_maf", as.character(sd_maf))|>
      stringr::str_replace_all( "opt\\$MAF", as.character(maf_threshold))|>
      stringr::str_replace_all( "opt\\$info", info) |>
      stringr::str_replace_all( "opt\\$source", source)

    writeLines( ldp_pop, fileConn)

    close(fileConn)

    ### Now make the submission script for the job

    fileConn<-file(paste(stringr::str_replace(outputs_dir, "/ess/", "/ess01/"),
                         user,
                         shortname,
                         paste0(shortname,"-sub.sh"),
                         sep="/"),
                   open="wb", encoding = "UTF-8")

    writeLines( paste0(c(
      "#!/bin/bash -l
#SBATCH --job-name="),shortname,c("
#SBATCH --output="),shortname,c(".scoring.LDpred2.out
#SBATCH --error="),shortname,c(".scoring.LDpred2.err
#SBATCH --time="),cpu_time,c("
#SBATCH --cpus-per-task="),cores,c("
#SBATCH --mem-per-cpu="),memory,c("
#SBATCH --account="),account,c("

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

## Set up job enviroment:
source /cluster/bin/jobsetup

set -o errexit

WORKdir="),outputs_subdir,c("

echo I am job $SLURM_JOBID

/cluster/software/EL9/amd/zen/easybuild/software/R/4.4.2-gfbf-2024a/bin/Rscript --vanilla "), paste0(shortname,"-ldpred_pipl.R")), fileConn, useBytes=T)


    close(fileConn)

    if(submit_ssh==TRUE){
      ## Connect to the cluster
      if(is.null(session)){ #Allow for the possibility that session is set once in a wrapper script

        message("Connecting to the cluster, please enter your password when prompted")
        session = ssh::ssh_connect(paste(user, host, sep="@"))

        }


      ## Navigate to the directory where the Rscript and submission script were saved, submit the submission script,
      ## wait and execute the squeue command
      ssh::ssh_exec_wait(session, command = c(paste0('cd ',outputs_subdir ),
                                              paste0('iconv -f "windows-1252" -t "UTF-8" ',shortname,'-sub.sh -o ',shortname,'-utf-sub.sh' ),
                                              paste0('rm -f ',shortname,'-sub.sh'),
                                              paste0('sbatch ',shortname,'-utf-sub.sh' ),'echo "Pausing to allow job to be scheduled, will then check status..." ' , 'sleep 10', 'squeue'))


      ## Other PGS software implementations can be added below
      if(software!="ldpred2"){

        stop("Only ldpred2 is currently supported as a PGS creation software in genotools.")
        }
      }
  }
}
