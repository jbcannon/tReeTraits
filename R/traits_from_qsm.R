#' Get center of mass from QSM (X, Y, Z)
#'
#' This function calculates the center of mass/volume from a QSM
#' by estimating the centroid of cylinder locations, each weighted
#' by their volume. For meaningful results, you may need to filter
#' out trunk sections only (e.g., `branching_order == 0`).
#' For center of mass, assumes constant density
#' within segments.
#' @param qsm qsm object loaded from `[load_qsm]`.
#' @importFrom dplyr mutate summarize
#' @export
get_center_of_mass = function(qsm) {
  qsm = dplyr::mutate(qsm,
                      X.mid = (startX + endX)/2,
                      Y.mid = (startY + endY)/2,
                      Z.mid = (startZ + endZ)/2)
  centroid = dplyr::summarize(qsm,
                              com.X = weighted.mean(X.mid, volume),
                              com.Y = weighted.mean(Y.mid, volume),
                              com.Z = weighted.mean(Z.mid, volume))
  return(centroid)
}
