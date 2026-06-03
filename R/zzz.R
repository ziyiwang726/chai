.onAttach <- function(libname, pkgname) {
  packageStartupMessage("chai version: ", packageVersion(pkgname))
}
