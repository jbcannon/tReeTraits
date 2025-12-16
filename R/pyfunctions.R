# ---------------------------
# Internal: Check required Python modules
# ---------------------------
#' Check if Python modules are available
#'
#' Internal helper that checks if specified Python modules are available in the
#' current Python environment. Issues a warning if any modules are missing.
#'
#' @param modules Character vector of Python module names to check.
#' @return Logical vector indicating availability of each module.
#' @keywords internal
#' @noRd
.check_python_modules <- function(modules) {
  if (missing(modules) || length(modules) == 0) stop("Provide Python module names.")
  available <- vapply(modules, reticulate::py_module_available, logical(1))

  if (!all(available)) {
    warning(
      "Missing Python modules: ", paste(modules[!available], collapse = ", "),
      "\nPlease install them manually. Example:\n",
      "reticulate::py_install(c('torch','numpy','robpy','PyTLidar'), pip = TRUE)"
    )
  }
  return(available)
}

# ---------------------------
# Check if Python 3.11 + PyTLidar environment is ready
# ---------------------------
#' Verify Python 3.11 environment for PyTLidar
#'
#' Internal helper that checks if Python 3.11 is available and whether required
#' packages (`torch`, `numpy`, `robpy`, `PyTLidar`) are installed.
#'
#' @return Logical indicating whether Python 3.11 and all required modules are available.
#' @keywords internal
#' @noRd
check_pytlidar_env <- function() {
  py <- tryCatch(reticulate::py_config(), error = function(e) NULL)
  if (is.null(py)) {
    message("❌ Python not found. Install Python 3.11 and required packages.")
    return(FALSE)
  }

  # Robust version parsing
  py_version <- as.character(py$version)
  version_parts <- suppressWarnings(as.numeric(unlist(strsplit(py_version, "\\."))))
  if (length(version_parts) < 2) {
    message("❌ Could not parse Python version: ", py$version)
    return(FALSE)
  }

  major <- version_parts[1]
  minor <- version_parts[2]

  if (major != 3 || minor != 11) {
    message("❌ Python version must be 3.11. Found: ", py$version)
    return(FALSE)
  }

  required <- c("torch", "numpy", "robpy", "PyTLidar")
  avail <- .check_python_modules(required)

  return(all(avail))
}

# ---------------------------
# Optional helper to guide user
# ---------------------------
#' Guide user to set up PyTLidar Python environment
#'
#' Internal helper that prints instructions to install Python 3.11 and
#' required Python packages for PyTLidar.
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
#' @noRd
setup_pytlidar_env <- function() {
  message(
    "Ensure Python 3.11 is installed and then install required packages:\n",
    "reticulate::py_install(c('torch','numpy','robpy','PyTLidar'), pip = TRUE)"
  )
  invisible(TRUE)
}

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
  if (!check_pytlidar_env()) {
    stop(
      "PyTLidar environment is not ready.\n",
      "Install Python 3.11 and packages: torch, numpy, robpy, PyTLidar.\n",
      "Example:\n",
      "reticulate::py_install(c('torch','numpy','robpy','PyTLidar'), pip = TRUE)"
    )
  }

  # --- Prepare output directory ---
  if (is.null(output_dir)) {
    output_dir <- file.path(tempdir(check = TRUE), "QSM_tmp")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)
  } else if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # --- Process LAS file ---
  las <- lidR::readLAS(file, filter = paste0("-thin_with_voxel ", resolution))
  las <- lidR::filter_poi(las, Intensity > intensity_threshold)
  las <- tReeTraits::normalize_las(las)
  las <- tReeTraits::recenter_las(las, height = 1)
  las <- lidR::las_update(las)

  myTempFile <- tempfile(fileext = ".las", tmpdir = output_dir)
  lidR::writeLAS(las, myTempFile)

  myTempFile <- normalizePath(myTempFile, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  # --- Build PyTLidar command ---
  args <- c("-m", "PyTLidar.treeqsm", myTempFile,
            "--outputdirectory", output_dir,
            "--custominput",
            "--ipd", paste(patch_diam1, collapse = " "),
            "--minpd", paste(patch_diam2min, collapse = " "),
            "--maxpd", paste(patch_diam2max, collapse = " "))

  if (verbose) args <- c(args, "--verbose")

  system2(reticulate::py_config()$python, args = args, stdout = "", stderr = "")

  # --- Read results ---
  results_dir <- file.path(output_dir, "results")
  if (!dir.exists(results_dir)) results_dir <- output_dir

  qsm_files <- list.files(results_dir, pattern = "^cylinder.*\\.txt$", recursive = TRUE, full.names = TRUE)
  if (length(qsm_files) == 0) stop("No QSM result files found.")

  qsm_fit_files <- list.files(results_dir, pattern = "^treedata.*\\.txt$", recursive = TRUE, full.names = TRUE)
  qsm_fits <- do.call(rbind, lapply(qsm_fit_files, parse_pytlidar_patch_params))
  qsm_fits$file <- qsm_fit_files
  qsm_fits <- dplyr::mutate(dplyr::rowwise(qsm_fits),
                            result = mean(c(AverageCylinderPointDistance_Trunk_mm, AverageCylinderPointDistance_BranchOrder1_mm)))
  qsm_fits <- dplyr::arrange(qsm_fits, result)

  best_qsm <- gsub('/treedata_', '/cylinder_', qsm_fits$file[1])
  qsm <- read_qsm_PyTLidar(best_qsm)
  parameters <- qsm_fits[1, grep('PatchDiam', names(qsm_fits))]
  parameters$fit_mm <- as.numeric(qsm_fits$result[1])

  list(qsm_pars = dplyr::select(qsm_fits, !file), qsm = qsm)
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

  cyl <- cyl %>%
    dplyr::mutate(
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
    ) %>%
    dplyr::select(
      startX, startY, startZ,
      endX, endY, endZ,
      cyl_ID, parent_ID, extension_ID,
      radius_cyl, length, volume,
      branching_order
    )

  cyl
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
