## run on maximus!!!

cran_mirror = "http://ftp5.gwdg.de/pub/misc/cran/"

args <- commandArgs(trailingOnly = TRUE)

## DEBUG ONLY
## local_library="/data/manke/repository/scripts/snakemake_workflows/R/x86_64-unknown-linux-gnu-library/3.3/"
## args = c("/home/kilpert/git/snakemake_workflows/shared/tools/install_R_packages.txt")
## .libPaths(c(local_library, .libPaths()[2:length(.libPaths())]))

packages_file = args[[1]]

## ## set local library
## ##local_library = file.path( "~/R", paste(R.Version()$platform, "-library", sep=""), paste(R.Version()$major, ".", unlist(strsplit(R.Version()$minor, "[.]"))[1], sep="") )
## ##.libPaths(local_library) 
## ## OR:
## local_library = .libPaths()[[1]]
## print(paste("Local library:", local_library))

local_library = .libPaths()[[1]]
local_library

## create dir
# if ( ! dir.exists(local_library) ) {
#   dir.create(local_library, recursive=T)
# }

source("https://bioconductor.org/biocLite.R")
biocLite()

sessionInfo()

installed.packages(lib.loc=local_library)
old.packages(repos=cran_mirror, lib.loc=local_library)

user_packages = as.character(unlist(read.table(packages_file)))
user_packages


## ## rebuild all packages
## update.packages(repos=cran_mirror, lib.loc=local_library, checkBuilt=TRUE)


## CRAN ########################################################################
print("Old packages:")
old.packages(repos=cran_mirror, lib.loc=local_library)
update.packages(repos=cran_mirror, lib.loc=local_library)
installed.packages(lib.loc=local_library)


## Bioconductor ################################################################
source("https://bioconductor.org/biocLite.R")
biocLite(lib.loc=local_library)


## Install packages from 1) CRAN or 2) bioconductor
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages(lib=local_library)[, "Package"])]
  if (length(new.pkg)) {
    tryCatch( install.packages(new.pkg, dependencies=TRUE, repos=cran_mirror, lib.loc=local_library),
              finally = biocLite(new.pkg, lib.loc=local_library)
    )
  }
  sapply(pkg, require, character.only=TRUE)
}

ipak(user_packages)

sessionInfo()

installed_packages = installed.packages(lib.loc=local_library)[,c(1)]


################################################################################

missing_packages = setdiff(user_packages, installed_packages)
if ( length(missing_packages) > 0 ) {
  cat("\nMissing packages (NOT installed):\n")
  print( missing_packages )
}
