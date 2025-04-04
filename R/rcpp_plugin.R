# Automatically registers a plugin for use with
# Rcpp::cppFunction

## REGISTER RCppPlugin
##' Make the plug in for inline Cxx
##'
##' Make the plug in for inline Cxx
##' @param ... I have no idea.
inlineCxxPlugin <- Rcpp::Rcpp.plugin.maker(
  include.before = "#include <FLasherEMSRR.h>", 
    ## FIND dll/so file at pkg root
  libs = paste0(find.package('FLasherEMSRR'), "/libs/", .Platform$r_arch,"/FLasherEMSRR", .Platform$dynlib.ext),
  package = "FLasherEMSRR"
)
