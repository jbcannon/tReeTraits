library(reticulate)

install_PyTLidar()

py_path <- system.file("python", "qsm_runner.py", package = "tReeTraits")
qsm <- reticulate::import_from_path("qsm_runner", path = dirname(py_path))
qsm$run_qsm('inst/extdata/tree_0129.laz')

# I'd like this to output a csv or other info representing the QSM. could save it to
# temp file as well.
