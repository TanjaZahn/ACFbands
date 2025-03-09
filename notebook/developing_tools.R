library(roxygen2) # In-Line Documentation for R 
library(devtools) # Tools to Make Developing R Packages Easier
library(testthat) # Unit Testing for R
library(usethis)  # Automate Package and Project Setup


# # Naming the package ---------------------------------------------------------
# 
# library(available) # Check if the Title of a Package is Available, 
# # Appropriate and Interesting
# # Check for potential names
# available::suggest("Easily extract information about your sample")
# available::suggest("Uncertainty quantification in forecast comparisons")
# 
# # Check whether it's available
# available::available("ACFbands", browse = FALSE)
# available::available("UQforecasts", browse = FALSE)
# available::available("UQeval", browse = FALSE)

# Don't show specific files ----------------------------------------------------

use_git_ignore("notebook*", directory = ".")

# Importing packages ----------------------------------------------------------

library(usethis)  # Automate Package and Project Setup
use_package("ggplot2")
use_import_from(package = "mvtnorm", "qmvnorm")
use_import_from(package = "dplyr", "lead")
use_import_from(package = "dplyr", "lag")
use_import_from(package = "stats", "qchisq")
use_import_from(package = "stats", "qnorm")


# Create descriptions for your functions ----------------------------------------

library(roxygen2) # Read in the roxygen2 R package
roxygenise()      # Builds the help files


# Load the package -------------------------------------------------------------

# Load the package
library(devtools)
load_all(".")


# Check the package ------------------------------------------------------------

# The following function runs a local R CMD check
devtools::check()

# # Check for CRAN specific requirements
# rhub::check_for_cran()
# 
# # Check for win-builder
# devtools::check_win_devel()


# Create a readme file --------------------------------------------------------

usethis::use_readme_rmd()

devtools::build_readme()




