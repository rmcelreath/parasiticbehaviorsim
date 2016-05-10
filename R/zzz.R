.onLoad <- function(libname, pkgname) { }

.onAttach <- function(...) {
  theLib <- dirname(system.file(package = "parasiticbehaviorsim"))
  pkgdesc <- packageDescription("parasiticbehaviorsim", lib.loc = theLib)
  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
  msg <- paste("parasiticbehaviorsim (Version ", pkgdesc$Version, ")\nInvoke ?parasiticbehaviorsim for documentation.", sep = "")
  packageStartupMessage(msg)
}