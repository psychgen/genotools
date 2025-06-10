.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0(
    "#########################################################################
    \nThis is genotools version ", packageVersion(pkgname),".\n
##########################################################################"))
}
