  # *Changelog for genotools package*

# genotools 0.3.1 (xxxxxxx)



# genotools 0.3.0 (14.11.24)

* This version contains major updates to allow polygenic scores to be made using ldpred2
* Specifically:
         1. the make_prsice and make_prsice_batch functions are replaced by make_pgs, with options for
        which software should be used. In this version, only LDpred is supported here, but make_pgs will be updated to be backwards compatible with PRSice2 in the next version
         2. R code for QCing summary statistics and running LDpred is added as a data object (ldpred.pipl)
         3. available_pgs and fetch_pgs now have default options for retrieving ldpred2 scores; addtionally, both functions refer to a (read-only) validated directory of scores by default, but also have the option to pull scores from a user specified raw-output directory
         4. process_pgs has minor changes to work with ldpred2 scores
         5. the vignette is updated with an overview of the new workflow
         6. the pgs_search helper is added

* Laura Hegemann became the joint developer of genotools on this version

# genotools 0.3.1 (15.04.2025) 


* Updates:
	1. bug fixes to avail_PGS and make_PGS
	2. updates to QC procedures Primarily adding: 
		A) slightly more restrictive QC thresholds (QC1 and QC2 in Privé, Arbel, et al. (2022))
		B) adjusting SD estimates for imputation quality, option to estimate sd from allele frequencies provided in summary statistics (compared to the default which estimates from the reference panel)
		C) Fixed a bug in the creation of allele frequency / MAF column
	3. Additional options added to make_PGS
		A) can now enter two allele frequency / MAF columns for case/controls 
		B) submit_ssh which is a TRUE/FALSE option to submit the job to cluster (or just create the directory, R script, and .sh submission script) 
		C) added the ability to provide a Z statistic as the stat_type
	4. Created 18 PGS for the repository using the new QC pipeline
	5. Restructured the PGS repository. Old ldpred2 scores can be found in ldpred2/archive 
	6. Updates to documentation and outline/sources for QC and creation of scores is provided
	7. R version used on cluster upgraded. 

	


