## NExt steps for QSM
# extract goodness of fit more smartly
# - address singular matrix issue w/ hagood
# - build help for multiple optimization types
# output/assess fit for multiple parameters


library(reticulate)

install_PyTLidar()
reticulate::use_condaenv("r-reticulate-3.11", required = TRUE)


file = 'inst/extdata/tree_0129.laz'
ipd = NULL; pdmin = NULL; pdmax = NULL; brad1 = NULL; brad2 = NULL; plot = NULL

run_pytlidar <- function(file,
                         envname = "r-reticulate-3.11",
                         ipd = NULL, pdmin = NULL, pdmax = NULL,
                         brad1 = NULL, brad2 = NULL,
                         plot = NULL) {



  if (!file.exists(file)) stop("File does not exist: ", file)

  # Activate the correct Python environment
  reticulate::use_condaenv("r-reticulate-3.11", required = TRUE)

  # Upgrade pip inside the env first (sometimes helps on Windows)
  reticulate::py_run_string("import pip; pip.main(['install', '--upgrade', 'pip'])")

  # Install robpy and dependencies
  reticulate::py_install(c("robpy", "torch"), envname = "r-reticulate-3.11", pip = TRUE)

  reticulate::py_module_available('robpy')


  # Upgrade pip inside the env first (sometimes helps on Windows)
  py_run_string("import pip; pip.main(['install', '--upgrade', 'pip'])")

  # Install robpy and dependencies
  py_install(c("robpy", "torch"), envname = "r-reticulate-3.11", pip = TRUE)

  py <- reticulate::py_config()$python

  # Build arguments list for Python module
  args <- c()
  if (!is.null(ipd))   args <- c(args, "--ipd", ipd)
  if (!is.null(pdmin)) args <- c(args, "--pdmin", pdmin)
  if (!is.null(pdmax)) args <- c(args, "--pdmax", pdmax)
  if (!is.null(brad1)) args <- c(args, "--brad1", brad1)
  if (!is.null(brad2)) args <- c(args, "--brad2", brad2)
  if (!is.null(plot))  args <- c(args, "--plot", plot)
  args <- c(args, file)  # add the LAS/LAZ file last

  # Build command
  cmd <- c("-m", "PyTLidar.treeqsm", file)
  system2(py, args = cmd)

  # Run the Python module
  reticulate::py_run_module("PyTLidar.treeqsm", args = args)

  message("PyTLidar TreeQSM finished processing: ", file)

  invisible(TRUE)










py_path <- system.file("python", "qsm_runner.py", package = "tReeTraits")
qsm <- reticulate::import_from_path("qsm_runner", path = dirname(py_path))
qsm$run_qsm('inst/extdata/tree_0129.laz')

