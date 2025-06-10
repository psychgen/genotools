#' Make a batch of polygenic scores using the \code{make_pgs} function
#'
#' \code{make_pgs_batch} is used to make a batch of polygenic scores based
#' on a data.frame of inputs (see \code{help(make_pgs)} for main inputs)
#'
#'  \code{make_pgs} creates and runs a job script to make a single PGS based on
#'  user-supplied arguments containing information about column names etc.
#'  \code{make_pgs_batch} is a wrapper script to allow a batch of PGS to be made at
#'  the same time, woth
#'
#' @param inputs a dataframe containing metadata to allow PGS to be created. The
#' formatting should should match genotools::pgs_metadata
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
#' @param genotype_dir location of genotype data on cluster; defaults for mobapsychgen_v1
#' @param genotype_data name of genotype data on cluster, formatted as an RDS file containing the subset
#' of the MoBa genetic data that overlaps with HapMap3 SNPs...; defaults set for mobapsychgen_v1
#' @param ldref_dir location of the HAPMAP3+ and ld blocks needed by ldpred2
#' @param ... arguments to pass to internal functions, see ?make_pgs
#' @export
#' @importFrom dplyr "%>%"


make_pgs_batch <- function(inputs,
                           user,
                     host="p471-hpc-01.tsd.usit.no",
                     software = "ldpred2",
                     shortname,
                     account="p471",
                     cpu_time="05:00:00",
                     memory="8G",
                     cores=32,
                     lib="//ess/p471/data/durable/common/rlibs/4-2-1",
                     sumstats_dir="//ess/p471/cluster/p/p471/cluster/common/gwas_sumstats/",
                     outputs_dir="//ess/p471/data/durable/common/pgs_directory",
                     genotype_dir="//ess/p471/data/durable/data/genetic_data/MoBaPsychGen_v1/",
                     genotype_data="genoHapMap3plus_N200k.rds",
                     ldref_dir="//ess/p471/data/durable/data/genetic/ldref/",
                     maf_threshold = 0.01,
                     ...
){

  # Connect SSH
  message("Connecting to the cluster, please enter your password when prompted")
  session = ssh::ssh_connect(paste(user, host, sep="@"))

  # Set up inputs dataframe

  inputs$user = user
  inputs$host = host
  inputs$software = software
  inputs$account = account
  inputs$cpu_time = cpu_time
  inputs$memory = memory
  inputs$cores = cores
  inputs$lib = lib
  inputs$sumstats_dir= sumstats_dir
  inputs$outputs_dir = outputs_dir
  inputs$genotype_dir = genotype_dir
  inputs$genotype_data = genotype_data
  inputs$ldref_dir = ldref_dir
  inputs$maf_threshold = maf_threshold

inputs = inputs %>%
  dplyr::select(user, host, software, shortname, account, cpu_time, memory, cores, lib, sumstats_dir, outputs_dir, sumstats_filename,
                genotype_dir, genotype_data, ldref_dir, case_control, A1, A2, stat, stat_type, pvalue, SNP, CHR, BP, SE, N, info, MAF,
                maf_threshold,source,pmid)

safe_pgs <- purrr::safely(make_pgs) #wraps make_pgs in safely so if one PGS fails it moves on to the next instead of failing completely

# Make the scores

inputs %>%
  purrr::pwalk(safe_pgs,
               pgs_meta=genotools::pgs_metadata,
               session=session )


}
