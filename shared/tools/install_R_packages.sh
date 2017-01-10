R="/package/R-3.3.1/bin/R"
local_library="/data/manke/repository/scripts/snakemake_workflows/R/x86_64-unknown-linux-gnu-library/3.3/"
packages_list="install_R_packages.txt"

[ -d $local_library ] || mkdir -p $local_library

( export R_LIBS_USER=${local_library} && cat install_R_packages.R | $R --vanilla --args $packages_list ) 2>&1 | tee $local_library/LOG
