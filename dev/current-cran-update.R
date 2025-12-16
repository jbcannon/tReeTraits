### where i left off
## trying to switch this to let users select their own python environment
## below code will check for python
## need to add code or function to check for conda, setup a conda env 3.11,
## and use it.
## need to take care not to setup environments, packages, etc. automatically as cran wont allow it.
## then use check pytlidar again.


# ---------------------------
# Minimal example for users
# ---------------------------

# Install package from GitHub (requires remotes or pak)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("jbcannon/tReeTraits")  # replace with your repo

# Load your package
library(tReeTraits)


check_pytlidar_env()

# Check if Python 3.11 + PyTLidar environment is ready
if (!check_pytlidar_env()) {
  message("Python environment not ready. Run setup:")
  setup_pytlidar_env()
}

# If environment is ready, you can run QSM:
# qsm_result <- run_treeqsm("path/to/tree.laz")
# head(qsm_result$qsm)






