# Function to return the location of the FLasher library
# Used by other packages that link to FLasher, for example FLasherTests
# Just need to include in Makevars:
#PKG_LIBS = `$(R_HOME)/bin/Rscript -e "FLasher:::FLasherLib()"` 
FLasherLib <- function(){
    location <- file.path(find.package('FLasherEMSRR'), "libs", paste0("FLasherEMSRR", .Platform$dynlib.ext))
   cat(location) 
}
