#' This function installs PyTLidar into the active Python environment
#' managed by reticulate. By default it will create or use a conda environment.
#'
#' @param method Installation method ("auto", "conda", "virtualenv").
#' @param conda Path to conda binary, or "auto".
#'#' @examples
#' \dontrun{
#'   install_pytlidar()
#' }
#' @export
install_PyTLidar <- function(  method = c("auto", "virtualenv", "conda"),
                               conda = "auto") {
  reticulate::py_install("PyTLidar", method = method, conda = conda)
}
install_PyTLidar()
