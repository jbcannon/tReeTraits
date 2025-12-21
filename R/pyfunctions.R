
#' Run TreeQSM on a LAS/LAZ file using PyTLidar
#'
#' Processes a point cloud, filters and normalizes it, then runs the PyTLidar TreeQSM model
#' to reconstruct cylinders representing the tree structure.
#'
#' @param file Path to the input LAS or LAZ file.
#' @param output_dir Directory to save output files; temporary if NULL.
#' @param intensity_threshold Minimum point intensity to retain.
#' @param resolution Thinning voxel size (m).
#' @param patch_diam1 Numeric vector of patch diameter 1 parameters.
#' @param patch_diam2min Numeric vector of minimum patch diameter 2.
#' @param patch_diam2max Numeric vector of maximum patch diameter 2.
#' @param verbose Logical; whether to print details during processing.
#' @return A list with elements:
#'   \describe{
#'     \item{qsm_pars}{Data frame of patch parameters and fit metrics.}
#'     \item{qsm}{Data frame of cylinder-level QSM output.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' run_treeqsm("example_tree.laz")
#' }
run_treeqsm <- function(
    file,
    output_dir = NULL,
    intensity_threshold = 40000,
    resolution = 0.02,
    patch_diam1 = c(0.05, 0.10),
    patch_diam2min = c(0.04, 0.05),
    patch_diam2max = c(0.12, 0.14),
    verbose = TRUE
) {
  if (!check_pytlidar_setup()) {
    stop(
      "PyTLidar environment not ready.\n",
      "Run check_pytlidar_setup() and follow the printed instructions."
    )
  }

  message("✅ Preparing output directory...")
  if (is.null(output_dir)) {
    output_dir <- file.path(tempdir(check = TRUE), "QSM_tmp")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit({
      try(unlink(output_dir, recursive = TRUE, force = TRUE))
    }, add = TRUE)
  } else if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  message("✅ Reading and preprocessing LAS file...")
  las <- lidR::readLAS(file, filter = paste0("-thin_with_voxel ", resolution))
  las <- lidR::filter_poi(las, Intensity > intensity_threshold)
  las <- normalize_las(las)
  las <- recenter_las(las, height = 1)
  las <- lidR::las_update(las)

  # las_points is your normalized, recentered LAS point matrix
  np <- reticulate::import("numpy")
  P <- as.matrix(las@data[, c("X", "Y", "Z")])
  P <- np$array(P)  # This is what you pass to Python

  # Build inputs dictionary for PyTLidar
  define_input <- reticulate::import("PyTLidar.Utils.define_input", delay_load = TRUE)$define_input
  inputs_list <- define_input(P, 1, 1, 1)  # 1 tree, 1 model, 1? (use standard args)
  inputs <- inputs_list[[1]]  # get first dict

  # overwrite just the patch diameters
  patch_diam1 <- c(patch_diam1)
  patch_diam2min <- c(patch_diam2min)
  patch_diam2max <- c(patch_diam2max)
  inputs$PatchDiam1 <- np$array(patch_diam1)
  inputs$PatchDiam2Min <- np$array(patch_diam2min)
  inputs$PatchDiam2Max <- np$array(patch_diam2max)
  inputs$BallRad1 <- np$array(patch_diam1 + 0.01)
  inputs$BallRad2 <- np$array(patch_diam2max + 0.01)

  message("✅ Running PyTLidar TreeQSM...")
  # import necessary modules
  pytlidar <- reticulate::import("PyTLidar", delay_load = TRUE)
  treeqsm <- reticulate::import("PyTLidar.treeqsm", delay_load = TRUE)

  res <- tryCatch(
    treeqsm$treeqsm(P, inputs, results_location = output_dir),
    error = function(e) stop("Error running PyTLidar TreeQSM: ", e$message)
  )

  message("✅ Reading QSM results...")
  results_dir <- file.path(output_dir, "results")
  if (!dir.exists(results_dir)) results_dir <- output_dir

  qsm_files <- list.files(results_dir, pattern = "^cylinder.*\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(qsm_files) == 0) stop("No QSM result files found.")

  qsm_fit_files <- list.files(results_dir, pattern = "^treedata.*\\.txt$", recursive = TRUE, full.names = TRUE)
  qsm_fits <- do.call(rbind, lapply(qsm_fit_files, .parse_patch_params))
  qsm_fits$file <- qsm_fit_files
  qsm_fits <- dplyr::mutate(dplyr::rowwise(qsm_fits),
                            result = mean(c(AverageCylinderPointDistance_Trunk_mm, AverageCylinderPointDistance_BranchOrder1_mm)))
  qsm_fits <- dplyr::arrange(qsm_fits, result)

  best_qsm <- gsub('/treedata_', '/cylinder_', qsm_fits$file[1])
  qsm <- read_qsm(best_qsm)
  parameters <- qsm_fits[1, grep('PatchDiam', names(qsm_fits))]
  parameters$fit_mm <- as.numeric(qsm_fits$result[1])

  list(qsm_pars = dplyr::select(qsm_fits, !file), qsm = qsm)
}



#' Read a PyTLidar QSM cylinder file
#'
#' Reads a PyTLidar-generated cylinder file and converts it into a tidy data frame
#' with start/end coordinates, radius, length, volume, and branching order.
#'
#' @param cyl_file Path to the PyTLidar cylinder output file (.txt).
#' @return A data frame with columns startX, startY, startZ, endX, endY, endZ,
#'   cyl_ID, parent_ID, extension_ID, radius_cyl, length, volume, branching_order.
#' @export
read_qsm <- function(cyl_file) {
  headers <- stringr::str_split(readLines(cyl_file, n = 1), "\t")[[1]]
  cyl <- readr::read_delim(cyl_file, skip = 1, col_names = FALSE, delim = "\t")

  new_headers <- c(
    headers[1:2],                       # radius, length
    paste0("start", c("X","Y","Z")),    # start coordinates
    paste0("dir", c("X","Y","Z")),      # direction vector
    headers[5:length(headers)]
  )
  colnames(cyl) <- new_headers

  cyl <- dplyr::mutate(cyl,
                       endX = startX + dirX * `length (m)`,
                       endY = startY + dirY * `length (m)`,
                       endZ = startZ + dirZ * `length (m)`,
                       cyl_ID = dplyr::row_number(),
                       parent_ID = parent,
                       extension_ID = extension,
                       length = `length (m)`,
                       radius_cyl = `radius (m)`,
                       volume = pi * radius_cyl^2 * length,
                       branching_order = branch_order
  )
  cyl = dplyr::select(cyl,
                      startX, startY, startZ,
                      endX, endY, endZ,
                      cyl_ID, parent_ID, extension_ID,
                      radius_cyl, length, volume,
                      branching_order
  )

  cyl
}



#' Check and guide setup for PyTLidar
#'
#' PyTLidar package is required to create QSMs. This function checks whether a
#' Python 3.11 environment and the PyTLidar Python package
#' are available via \pkg{reticulate}. If anything is missing, prints
#' step-by-step instructions for installing Miniconda, creating a Python 3.11
#' environment, and installing PyTLidar.
#' @return Invisibly \code{TRUE} if the Python environment and PyTLidar are
#'   available and ready. Otherwise, prints instructions and returns
#'   \code{FALSE} invisibly.
#' @export
#' @examples
#' \dontrun{
#' # Check your Python + PyTLidar setup
#' check_pytlidar_setup()
#'
#' # If the environment is not ready, follow the printed instructions:
#' # 1. Install Miniconda if needed
#' #    reticulate::install_miniconda()
#' # 2. Create Python 3.11 environment
#' #    reticulate::conda_create("pytlidar", python = "3.11")
#' # 3. Tell reticulate to use it
#' #    reticulate::use_condaenv("pytlidar", required = TRUE)
#' # 4. Install PyTLidar in that environment
#' #    reticulate::py_install(c("torch", "numpy", "robpy", "PyTLidar"), pip = TRUE)
#' # 5. Restart R and re-run check_pytlidar_setup()
#' }
#' @keywords utilities setup python reticulate
check_pytlidar_setup <- function() {

  cat("🔍 Checking Python configuration...\n\n")

  # Attempt to select the pytlidar env if it exists
  tryCatch(
    reticulate::use_condaenv("pytlidar", required = TRUE),
    error = function(e) {
      cat("❌ Python environment 'pytlidar' not found.\n")
      .print_miniconda_instructions()
      return(invisible(FALSE))
    }
  )

  # ---- Step 1: Ensure Python is initialized ----
  cfg <- tryCatch(
    reticulate::py_config(),
    error = function(e) NULL
  )

  if (is.null(cfg)) {
    cat("❌ No Python detected by reticulate.\n\n")
    .print_miniconda_instructions()
    return(invisible(FALSE))
  }

  # ---- Step 2: Python version ----
  version <- as.character(cfg$version)
  major_minor <- paste(strsplit(version, "\\.")[[1]][1:2], collapse = ".")

  if (!identical(major_minor, "3.11")) {
    cat("❌ Incompatible Python version detected.\n\n")
    cat("Found Python version: ", version, "\n", sep = "")
    cat("Required version: Python 3.11\n\n")
    .print_miniconda_instructions()
    return(invisible(FALSE))
  }

  cat("✅ Python 3.11 detected.\n\n")

  # ---- Step 3: pytlidar availability ----
  ok <- FALSE
  try({
    reticulate::import("PyTLidar", delay_load = TRUE)
    ok <- TRUE
  }, silent = TRUE)

  if (!ok) {
    cat("❌ Python package 'PyTLidar' is not installed.\n\n")
    .print_pytlidar_install_instructions()
    return(invisible(FALSE))
  }

  cat("✅ PyTLidar is installed and importable.\n")
  cat("🎉 Python environment is ready.\n")

  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.print_miniconda_instructions <- function() {
  cat(
    "Recommended setup using reticulate-managed Miniconda:\n",
    "1. Install Miniconda from R:\n",
    "     reticulate::install_miniconda()\n",
    "   ⚠️ IMPORTANT: Restart R immediately after this step!\n",
    "2. In the new R session, create a Python 3.11 environment:\n",
    "     reticulate::conda_create('pytlidar',
    \tpackages = 'python=3.11',
    \tchannels = c('conda-forge'))\n",
    "3. Install required packages into this environment:\n",
    "     reticulate::conda_install('pytlidar', packages = c('torch','numpy','robpy','PyTLidar'), pip = TRUE)\n",
    "4. Tell reticulate to use this environment:\n",
    "     reticulate::use_condaenv('pytlidar', required = TRUE)\n",
    "5. Finally, run:\n",
    "     check_pytlidar_setup()\n",
    sep = ""
  )
}


#' @keywords internal
#' @noRd
.parse_patch_params <- function(file_path) {
  lines <- readLines(file_path)
  lines <- lines[nzchar(lines)]  # remove empty lines

  start_idx <- which(grepl("^PatchDiam1", lines))
  if (length(start_idx) == 0) stop("PatchDiam1 not found in file.")

  lines <- lines[start_idx:length(lines)]
  lines <- lines[!grepl("^AverageSurfaceCoverage_", lines)]

  params <- list()
  for (line in lines) {
    split_line <- strsplit(line, "\\s+")[[1]]
    name <- split_line[1]
    value <- suppressWarnings(as.numeric(split_line[-1]))
    if (all(is.na(value))) value <- paste(split_line[-1], collapse = " ")
    params[[name]] <- value
  }

  as.data.frame(params, stringsAsFactors = FALSE)
}

#' Save QSM results and patch parameters
#'
#' Writes QSM cylinder data and parameter summaries to tab-delimited text files.
#'
#' @param qsm Output from \code{run_treeqsm()}; a list with \code{qsm} and \code{qsm_pars}.
#' @param name Base name to use for output files.
#' @param output_dir Directory where files are written (defaults to current working directory).
#' @return Invisibly returns a list with file paths: \code{qsm} and \code{pars}.
#' @export
write_qsm <- function(qsm, name, output_dir = getwd()) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  qsm_file <- file.path(output_dir, paste0(name, "_qsm.txt"))
  pars_file <- file.path(output_dir, paste0(name, "_qsm_pars.txt"))

  # Write cylinder data
  write.table(qsm$qsm, file = qsm_file, row.names = FALSE, sep = "\t", quote = FALSE)

  # Write parameters
  pars <- qsm$qsm_pars
  df <- as.data.frame(t(unlist(pars)))
  colnames(df) <- names(pars)
  write.table(df, file = pars_file, row.names = FALSE, sep = "\t", quote = FALSE)

  message("Outputs written to:\n", qsm_file, "\n", pars_file)
  invisible(list(qsm = qsm_file, pars = pars_file))
}


