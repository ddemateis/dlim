###
### R routines for the R package dlim (c)
#
#' @importFrom utils packageDescription

.onAttach <- 
  function(lib, pkg) {
    #
    ################################################################################
    #
    meta <- packageDescription("dlim")
    attachmsg <- paste("This is dlim ",meta$Version,
                       ". For details: help(`dlim-package`) and vignette('dlimOverview').",
                       sep="")
    packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
  }