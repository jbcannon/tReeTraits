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

  # Instally roby and PyTLidar if needed
  if (!reticulate::py_module_available("torch")) {
    reticulate::py_install(packages='torch', envname=envname, method='conda', pip=TRUE)
  }
  if (!reticulate::py_module_available("robpy")) {
    reticulate::py_install(packages='robpy', envname=envname, method='conda', pip=TRUE)
  }
  if (!reticulate::py_module_available("PyTLidar")) {
    message("Installing PyTLidar from GitHub into ", envname, "...")
    reticulate::py_install(packages = "PyTLidar",
                           envname = envname,
                           method = "conda", pip = TRUE
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

qsm <- reticulate::import_from_path("qsm_runner", path = dirname(py_path))


#' Run Quantitative Structure Model (QSM) from PyTLidar
#'
#' This function runs the PyTLidar TreeQSM pipeline on a terrestrial LiDAR LAS file
#' from within R using the reticulate interface. Parameters can be adjusted to
#' control cover set sizes and cylinder fitting behavior.
#'
#' @param file Character string. Path to a LAS file to process.
#' @param ipd Numeric. Initial patch diameter (default = 0.05).
#' @param pdmin Numeric. Minimum patch diameter for second pass (default = 0.03).
#' @param pdmax Numeric. Maximum patch diameter for second pass (default = 0.12).
#' @param brad1 Numeric. Ball radius for first pass (default = 0.06).
#' @param brad2 Numeric. Ball radius for second pass (default = 0.13).
#' @param plot Integer (0 or 1). Flag for plot generation during runtime (default = 0).
#'
#' @return A data frame containing tree structural metrics output from TreeQSM.
#'
#' @examples
#' \dontrun{
#'   # Example LAS file
#'   las_file <- system.file("extdata", "example_pine.las", package = "mypkg")
#'
#'   # Run QSM with default parameters
#'   qsm_results <- run_qsm(las_file)
#'
#'   # Run with custom parameters
#'   qsm_results <- run_qsm(las_file, ipd = 0.08, pdmax = 0.15)
#' }
#'
#' @export
run_qsm <- function(file, ipd=0.05, pdmin=0.03, pdmax=0.12,
                    brad1=0.06, brad2=0.13, plot=0) {

  # find the installed path to the Python file
  py_path <- system.file("python", "qsm_runner.py", package = "tReeTraits")

  # import the module from that path
  qsm <- reticulate::import_from_path("qsm_runner", path = dirname(py_path))

  # run the Python function
  res <- qsm$run_qsm(file=file, ipd=ipd, pdmin=pdmin,
                     pdmax=pdmax, brad1=brad1, brad2=brad2, plot=plot)

  as.data.frame(res)
}

x = run_qsm('inst/extdata/tree_0129.laz')
