library(devtools)
library(roxygen2)

#Create package (directory, Rproj, and basic structure)
usethis::create_package("~/path")

#Setup git
usethis::use_git()

#Add fxn
use_r("sigma_plot")

#Test the package build
load_all()
#document()
check()

#add licence
use_gpl3_license()
check()

#Add package dependencies
use_package("dplyr")
use_package("tidyr")
use_package("ggplot2")
check()

#pipes
use_pipe()
check()

#Documentation
load_all()
document()
check()

#Add example date
use_data_raw()
use_data(model_result, overwrite = TRUE)

#Data documentation
use_r("model_result")

#Addtl documentation
use_readme_rmd()
build_manual(pkg=".", path="../")
