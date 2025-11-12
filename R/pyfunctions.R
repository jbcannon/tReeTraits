#' Check for and install Python modules in a reticulate environment
#'
#' This function checks if the specified Python modules are available in a
#' given reticulate-managed Conda environment, and installs any that are missing.
#'
#' @param modules Character vector of Python module names to check/install.
#' @param envname Character string. Name of the reticulate/conda environment to use.
#'   Default is "r-reticulate-3.11".
#' @param method Character string specifying installation method for reticulate::py_install.
#'   Default is "auto".
#' @export
#' @examples
#' \dontrun{
#' require_module(c("torch", "robpy", "PyTLidar"))
#' }
require_module <- function(modules, envname = "r-reticulate-3.11") {

  # --- Input safety checks ---
  if (missing(modules) || length(modules) == 0) {
    stop("Please provide a character vector of Python module names in `modules`.")
  }
  if (!is.character(modules)) {
    stop("`modules` must be a character vector of module names.")
  }
  if (!is.character(envname) || length(envname) != 1) {
    stop("`envname` must be a single character string specifying the conda environment.")
  }

  reticulate::py_install(c("robpy"), envname = "r-reticulate-3.11", pip = TRUE, method='conda')
  reticulate::py_install(c("torch"), envname = "r-reticulate-3.11", pip = TRUE, method='conda')


  # --- Loop through modules and install if missing ---
  for (m in modules) {
    if (!reticulate::py_module_available(m)) {
      message("Installing Python module '", m, "' in environment '", envname, "'...")
      reticulate::py_install(packages = m, envname = envname, method = 'conda', pip = TRUE)
    } else {
      message("Python module '", m, "' is already installed in '", envname, "'.")
    }
  }
}


#' Install PyTLidar and dependencies
#'
#' Sets up a Python conda environment and installs dependencies for PyTLidar if needed.
#' This includes Python 3.11, `torch`, `numpy`, `robpy`, and `PyTLidar` (via pip).
#'
#' @param envname Character. Name of the Conda environment to use or create.
#'   Default is `"r-reticulate-3.11"`.
#' @param python_version Character. Python version to install (default `"3.11"`).
#' @param method Character. Installation method for reticulate, either `"auto"` or `"conda"`.
#' @param conda Character. Path to the conda binary or `"auto"`.
#' @param reinstall Boolean If TRUE, delete conda environment first.
#' @return Invisibly returns the Python configuration used.
#' @export
install_PyTLidar <- function(envname = "r-reticulate-3.11",
                             python_version = "3.11",
                             method = c("auto", "conda"),
                             conda = "auto",
                             reinstall=FALSE) {
  method <- match.arg(method)

  # --- Check if conda is available ---
  conda_path <- tryCatch(reticulate::conda_binary(conda), error = function(e) NULL)
  if (is.null(conda_path) || !file.exists(conda_path)) {
    stop(
      "Conda installation not found.\n",
      "Please install Miniconda or Anaconda before continuing.\n\n",
      "Download Miniconda (recommended, lightweight): https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions\n",
      "Download Anaconda (full distribution): https://www.anaconda.com/download"
    )
  } else {
    message("Found conda at: ", conda_path)
  }

  # --- Delete and reinstall conda environment if needed
  if(reisntall) {
    if (reticulate::condaenv_exists(envname)) {
      message("Deleting conda environment '", envname, "' with Python ", python_version, "...")
      reticulate::conda_remove(envname)
    } else {
      message("Conda environment ", envname, " doesn't exist...")
    }
  }

  # --- Create or activate conda environment ---
  if (!reticulate::condaenv_exists(envname)) {
    message("Creating conda environment '", envname, "' with Python ", python_version, "...")
    reticulate::conda_create(envname, python_version = python_version)
  } else {
    message("Using existing conda environment: ", envname)
  }

  # --- Activate environment ---
  reticulate::use_condaenv(envname, required = TRUE)

  # --- Install required modules ---
  required_modules <- c("torch", "numpy", "robpy", "PyTLidar")
  installed <- vapply(required_modules, reticulate::py_module_available, logical(1))
  if (!all(installed)) {
    missing <- names(installed[!installed])
    message("Installing missing Python packages: ", paste(missing, collapse = ", "))
    reticulate::py_install(missing, envname = envname, method = method, pip = TRUE)
  } else {
    message("All required Python modules are already installed.")
  }

  # --- Confirm setup ---
  py <- reticulate::py_config()
  message("Using Python: ", py$python, " (version ", py$version, ")")

  invisible(py)
}


#' Check PyTLidar installation and environment
#'
#' Verifies that the conda environment and all required Python packages for PyTLidar are available.
#' If the environment or any package is missing, it returns FALSE with a message.
#'
#' @param envname Character. Name of the conda environment to check.
#' @param required_modules Character vector of required Python modules.
#' @return Logical indicating whether the environment is ready.
#' @export
check_PyTLidar_install <- function(envname = "r-reticulate-3.11",
                                   required_modules = c("torch", "numpy", "robpy", "PyTLidar")) {

  conda_path <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)
  if (is.null(conda_path) || !file.exists(conda_path)) {
    message("❌ Conda installation not found.")
    return(FALSE)
  }

  if (!reticulate::condaenv_exists(envname)) {
    message("❌ Conda environment '", envname, "' does not exist.")
    return(FALSE)
  }

  reticulate::use_condaenv(envname, required = TRUE)

  installed <- vapply(required_modules, reticulate::py_module_available, logical(1))
  if (!all(installed)) {
    missing <- names(installed[!installed])
    message("❌ Missing Python packages in '", envname, "': ", paste(missing, collapse = ", "))
    return(FALSE)
  }

  message("✅ PyTLidar environment '", envname, "' is properly installed and configured.")
  return(TRUE)
}


#' Run TreeQSM using PyTLidar
#'
#' This function runs the PyTLidar TreeQSM model from R using a LAS file.
#' It filters, normalizes, and centers the point cloud, then executes the PyTLidar command-line interface.
#'
#' @param file Path to the input LAS or LAZ file.
#' @param intensity_threshold Minimum point intensity to retain.
#' @param resolution Thinning voxel size (m).
#' @param output_dir Directory for results (temporary if NULL).
#' @param patch_diam1,patch_diam2min,patch_diam2max Numeric vectors of patch diameter parameters.
#' @param optimization Character. Optimization method for multi-parameter runs.
#' @param verbose Logical. Whether to print details during processing.
#' @param keep_results Logical. Whether to keep output files.
#' @return A QSM object read by \code{read_qsm_PyTLidar()}.
#' @export
run_treeqsm <- function(
    file,
    intensity_threshold = 40000,
    resolution = 0.02,
    output_dir = NULL,
    patch_diam1 = c(0.05, 0.10),
    patch_diam2min = c(0.04, 0.05),
    patch_diam2max = c(0.12, 0.14),
    optimization = "trunk+1branch_mean_dis",
    verbose = TRUE,
    keep_results = !is.null(output_dir)
) {

    envname <- "r-reticulate-3.11"

  # --- Check PyTLidar installation ---
  if (!check_PyTLidar_install(envname)) {
    message("Attempting to install PyTLidar environment...")
    install_PyTLidar(envname)
  }

  reticulate::use_condaenv(envname, required = TRUE)

  # --- Create output directory ---
  if (is.null(output_dir)) {
    output_dir <- file.path(tempdir(check = TRUE), "QSM_tmp")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  } else if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # --- Read and process LAS ---
  las <- lidR::readLAS(file, filter = paste0("-thin_with_voxel ", resolution))
  las <- lidR::filter_poi(las, Intensity > intensity_threshold)
  las <- tReeTraits::normalize_las(las)
  las <- tReeTraits::recenter_las(las, height = 1)
  las <- lidR::las_update(las)

  # --- Write to temporary file ---
  myTempFile <- tempfile(fileext = ".las", tmpdir = output_dir)
  lidR::writeLAS(las, myTempFile)

  # --- Normalize paths (Windows-safe) ---
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  myTempFile <- normalizePath(myTempFile, winslash = "/", mustWork = TRUE)

  # --- Build system call ---
  args <- c("-m", "PyTLidar.treeqsm", myTempFile,
            "--outputdirectory", output_dir,
            "--custominput",
            "--ipd", paste(patch_diam1, collapse = " "),
            "--minpd", paste(patch_diam2min, collapse = " "),
            "--maxpd", paste(patch_diam2max, collapse = " "))

  if (verbose) args <- c(args, "--verbose")

  multiple_parameters <- length(patch_diam1) > 1 ||
    length(patch_diam2max) > 1 ||
    length(patch_diam2min) > 1

  if (multiple_parameters) {
    args <- c(args, "--optimum", optimization)
    cat("Running QSM with multiple parameters\n",
        "...Patchdiam1 =", patch_diam1,
        "\n...PatchDiam2_min =", patch_diam2min,
        "\n...PatchDiam2_max =", patch_diam2max,
        "\n...Optimization =", optimization, "\n")
  } else {
    cat("Running QSM with single parameter set\n",
        "...Patchdiam1 =", patch_diam1,
        "\n...PatchDiam2_min =", patch_diam2min,
        "\n...PatchDiam2_max =", patch_diam2max, "\n")
  }

  # --- Run PyTLidar ---
  system2("python", args = args, stdout = "", stderr = "")

  # --- Locate output ---
  results_dir <- file.path(output_dir, "results")
  if (!dir.exists(results_dir)) results_dir <- output_dir

  qsm_files <- list.files(results_dir, pattern = "^cylinder.*\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(qsm_files) == 0) {
    warning("No QSM result files found in output directory.")
    return(invisible(NULL))
  }

  # Load QSM and parameters
  latest_file <- qsm_files[which.max(file.info(qsm_files)$mtime)]
  qsm <- read_qsm_PyTLidar(latest_file)
  parameters = extract_patch_params(latest_file, patch_diam1, patch_diam2min, patch_diam2max)

  if (!keep_results) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }

  return(list(qsm_pars=parameters, qsm=qsm))
}


#' Save QSM results and parameters to files
#'
#' Internal helper function to write QSM data and its parameters to text files.
#' Not intended for end-user calls.
#'
#' @param qsm Output from function `run_treeqsm`. A list containing `qsm` (data.frame of cylinder data) and `qsm_pars` (named list of parameters)
#' @param name Character. Base name to use for output files.
#' @param output_dir Character. Directory where files will be written. Defaults to current working directory.
#'
#' @return A list with file paths: `qsm` (QSM output) and `pars` (parameter file).
#' @export
write_qsm <- function(qsm, name, output_dir = getwd()) {

  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Build file paths
  qsm_file   <- file.path(output_dir, paste0(name, "_qsm.txt"))
  pars_file  <- file.path(output_dir, paste0(name, "_qsm_pars.txt"))

  # Write QSM data
  write.table(qsm$qsm, file = qsm_file, row.names = FALSE, sep = "\t", quote = FALSE)

  # Convert named list of parameters to 1-row data.frame
  pars <- qsm$qsm_pars
  df <- as.data.frame(t(unlist(pars)))
  colnames(df) <- names(pars)

  # Write parameters to file
  write.table(df, file = pars_file, row.names = FALSE, sep = "\t", quote = FALSE)

  # Inform user
  cat("Outputs written to:\n", qsm_file, "\n", pars_file, "\n")

  # Return file paths
  invisible(list(qsm = qsm_file, pars = pars_file))
}


#' Read and Convert PyTLidar QSM Cylinder File
#'
#' Reads a PyTLidar-generated QSM cylinder file (tab-delimited text) and converts it
#' into a tidy R data frame with separate start and end coordinates, cylinder IDs,
#' parent/extension mapping, radius, length, and volume, as used in tReeTraits
#'
#' @param cyl_file Character. Path to the PyTLidar QSM cylinder output file (.txt).
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{startX, startY, startZ}{Coordinates of the cylinder start point.}
#'   \item{endX, endY, endZ}{Coordinates of the cylinder end point. Computed as start + direction.}
#'   \item{cyl_ID}{Unique integer ID for each cylinder.}
#'   \item{parent_ID}{Parent cylinder ID (from PyTLidar output).}
#'   \item{extension_ID}{Extension ID (from PyTLidar output).}
#'   \item{radius_cyl}{Cylinder radius in meters.}
#'   \item{length}{Cylinder length in meters.}
#'   \item{volume}{Cylinder volume in cubic meters, computed as π * r² * length.}
#'   \item{branching_order}{Branching order from PyTLidar output.}
#' }
#'
#' @examples
#' \dontrun{
#' cyl <- read_qsm_PyTLidar("tree_QSM.txt")
#' head(cyl)
#' }
#'
#' @export
read_qsm_PyTLidar <- function(cyl_file) {
  # ------------------------------
  # Step 1: Read the header line
  # ------------------------------
  headers <- stringr::str_split(readLines(cyl_file, n = 1), "\t")[[1]]

  # ------------------------------
  # Step 2: Read data without headers
  # ------------------------------
  cyl <- readr::read_delim(cyl_file, skip = 1, col_names = FALSE, delim = "\t")

  # ------------------------------
  # Step 3: Fix column names
  # ------------------------------
  new_headers <- c(
    headers[1:2],                       # radius and length
    paste0("start", c("X","Y","Z")),    # start_point -> startX, startY, startZ
    paste0("dir", c("X","Y","Z")),      # axis_direction -> dirX, dirY, dirZ
    headers[5:length(headers)]          # remaining columns
  )
  colnames(cyl) <- new_headers

  # ------------------------------
  # Step 4: Compute end coordinates, IDs, and derived values
  # ------------------------------
  cyl <- cyl %>%
    dplyr::mutate(
      # End coordinates = start + direction
      endX = startX + dirX * `length (m)`,
      endY = startY + dirY * `length (m)`,
      endZ = startZ + dirZ * `length (m)`,

      # Cylinder identifiers
      cyl_ID       = dplyr::row_number(),
      parent_ID    = parent,
      extension_ID = extension,

      # Cylinder properties
      length      = `length (m)`,
      radius_cyl  = `radius (m)`,
      volume      = pi * radius_cyl^2 * length,
      branching_order = branch_order
    ) %>%
    dplyr::select(
      startX, startY, startZ,
      endX, endY, endZ,
      cyl_ID, parent_ID, extension_ID,
      radius_cyl, length, volume,
      branching_order
    )

  return(cyl)
}

