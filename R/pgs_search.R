#' PGS search
#'
#' \code{pgs_search} is a helper to find PGS that are available by
#' a string anywhere in the variable name or description
#'
#'
#' Detailed description..
#'
#' @param dataframe Where to look for variables -should be a data frame returned
#' by a call to \code{available_pgs}
#' @param string a vector of strings to search for
#' @param where if you know the column name in which the string should be found,
#' you can supply it here; the default is to search in all columns ("anywhere")
#' @export
#' @importFrom dplyr "%>%"

pgs_search <- function(dataframe, string, where="anywhere"){

  collapsed_string <- paste0(string, collapse="|") %>%  stringr::str_remove_all("\\*")

  if(where=="anywhere"){
    tmp <- dataframe %>%
      dplyr::filter( stringr::str_detect(Pheno_shortname, collapsed_string )|
                       stringr::str_detect(Phenotype, collapsed_string )|
                       stringr::str_detect(Sumstats_filename, collapsed_string )|
                       stringr::str_detect(Source, collapsed_string ))
  }else{
    tmp <- dataframe %>%
      dplyr::filter( stringr::str_detect(where, collapsed_string ))
  }

  return(tmp)

}



