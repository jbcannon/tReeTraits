#' Install PyTLidar in a reticulate-managed Python environment
#'
#' This function installs the PyTLidar package from GitHub into a conda
#' environment managed by reticulate. By default, it creates/uses the
#' "r-reticulate-3.11" environment with Python 3.11.
#'
#' @param envname Name of the conda environment to use/create. Default is "r-reticulate".
#' @param python_version Python version to use in the environment. Default is "3.11".
#' @export
#' @examples
#' \dontrun{
#'   install_PyTLidar()
#' }
#' @export
install_PyTLidar <- function(  envname = "r-reticulate-3.11",
                               method = c("auto", "conda"),
                               conda = "auto") {
  # --- Check if conda is available ---
  conda_path <- tryCatch(
    reticulate::conda_binary(),
    error = function(e) NULL
  )

  if (is.null(conda_path) || !file.exists(conda_path)) {
    stop(
      "Conda installation not found.\n",
      "Please install Miniconda or Anaconda before continuing.\n\n",
      "Download Miniconda (recommended, lightweight): https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions\n",
      "Download Anaconda (full distribution): https://www.anaconda.com/download\n"
    )
  } else {
    message("Found conda at: ", conda_path)
  }

  # Create or use conda environment
  if (!reticulate::condaenv_exists(envname)) {
    message("Creating conda environment ", envname, " with Python 3.11...")
    reticulate::conda_create(envname, python_version = "3.11")
  } else {
    message("Using existing conda environment: ", envname)
  }

  #Activeate Environment
  reticulate::use_condaenv(envname, required = TRUE)

  # Check if PyTLidar is already installed
  if (!reticulate::py_module_available("PyTLidar")) {
    message("Installing PyTLidar from GitHub into ", envname, "...")
    reticulate::py_install(
      packages = "git+https://github.com/Landscape-CV/PyTLidar.git",
      envname = envname,
      method = "conda",
      pip = TRUE
    )
  } else {
    message("PyTLidar is already installed in ", envname)
  }

  # Confirm installation
  py <- reticulate::py_config()
  message("Using Python: ", py$python, " (version ", py$version, ")")
  if (reticulate::py_module_available("PyTLidar")) {
    message("PyTLidar is successfully installed and available.")
  } else {
    warning("PyTLidar installation did not succeed. Please check your Python environment.")
  }
}

install_PyTLidar()
