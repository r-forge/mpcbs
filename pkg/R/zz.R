.First.lib <- function(libname, pkgname) {
  library.dynam(pkgname, package=pkgname, lib.loc=libname);

  pi <- utils::packageDescription(pkgname);
  packageStartupMessage("Loaded ", pkgname, 
                        " v", pi$Version, " (", pi$Date, ") ");
}

