#' PGS PCA
#'
#' \code{pgs_pca} implements the PGS-PCA method to return the first principle
#' component from PCA of a set of scores at a range of thresholds made using
#' a clumping and thresholding method such as prsice2 (not applicable for
#' scores from ldpred2)
#'
#' See https://doi.org/10.1002/gepi.22339
#'
#' @param geno_dat data frame with polgenic scores (probably from \code{process_pgs})
#' @param indid individual level id variable
#' @param pgs_var_stem string(s) for the stem that is in common for all variables
#' @param return some combination of "pcs", "r2s", and "loadings"
#' representing PGS for a given phenotype at some p-value threshold
#' @export
#' @importFrom dplyr "%>%"



pgs_pca <- function(dat,
                    indid,
                    pgs_var_stem,
                    return ="pcs"){

  # FROM Supplemental R Code for PGS-PCA paper https://doi.org/10.1002/gepi.22339

  ###### FUNCTION TO PERFORM PGS-PCA ##########################################
  # INPUTS:
  # dat = n x (p+1) dataframe of PGSs under different settings with first column as ID
  # x = label for PGS (typically named for phenotype (i.e. MDD))
  # OUTPUTS:
  # list of
  #  - data = dataframe with cols (ID , PGS-PCA1 , PGS-PCA2)
  #  - r2 = variance explained by each PC of the PGS matrix
  #  - loadings = the PC-loadings used to create PGS-PCA1
  pgs.pc <- function(dat,x){
    xo <- scale(as.matrix(dat[,-1]))  ## scale cols of matrix of only PGSs (remove ID)
    g <- prcomp(xo)   ## perform principal components
    pca.r2 <- g$sdev^2/sum(g$sdev^2)    ## calculate variance explained by each PC
    pc1.loadings <- g$rotation[,1];     ## loadings for PC1
    pc2.loadings <- g$rotation[,2]      ## loadings for PC2
    ## flip direction of PCs to keep direction of association
    ## (sign of loadings for PC1 is arbitrary so we want to keep same direction)
    if (mean(pc1.loadings>0)==0){
      pc1.loadings <- pc1.loadings*(-1)
      pc2.loadings <- pc2.loadings*(-1)
    }
    ## calculate PGS-PCA (outputs PC1 and PC2 even though PC1 sufficient)
    pc1 <- xo %*% pc1.loadings
    pc2 <- xo %*% pc2.loadings
    dat[,paste0(x,".pgs.pc")] <- scale(pc1)   ## rescales PGS-PCA1
    dat[,paste0(x,".pgs.pc2")] <- scale(pc2)  ## rescales PGS-PCA2
    return(list(data=dat,r2=pca.r2,loadings=pc1.loadings))    #####THIS LINE MODIFIED LJH AS REF TO DF WAS WRONG (SHOULD BE DAT)
  }


  pgs_pca_res <- dat %>%
    dplyr::select(indid) %>%
    tidyr::drop_na()

  for(pheno in pgs_var_stem){

    pgs_tmp <- dat %>%
      dplyr::select(indid, matches(pheno)) %>%
      tidyr::drop_na()
    tmp <-pgs.pc(pgs_tmp,pheno)

    pgs_pca_res <- pgs_pca_res %>%
      dplyr::left_join(tmp$data %>%
                         dplyr::select(indid, tidyr::matches("pc")) %>%
                         dplyr::select(-tidyr::matches("pc2")))


  }

  #Join back to main dataset

  dat <- dat %>%
    dplyr::left_join(pgs_pca_res)

  if(all(return=="pcs")){
    return(dat)
  }else{
    ret=list()
    if(any(return == "pcs")){
      ret[["pcs"]] <- dat
    }
    if(any(return == "r2s")){
      ret[["r2s"]] <- tmp$r2
    }
    if(any(return == "loadings")){
      ret[["loadings"]] <- tmp$loadings
    }
    return(ret)

  }


}
