

# Preparation
normalize_and_recenter_las = function(las, quantile=c(0.001)) {
  #normalize
  ground_level = quantile(las$Z, quantile)
  las@data[, Z := las$Z - ground_level]

  #recenter
  bole = filter_poi(las, Z < 3)
  centroid = c(mean(bole$X), mean(bole$Y))
  las@data[, Y := las$Y - mean(bole$Y)]
  las@data[, X := las$Y - mean(bole$X)]
  las = lidR::las_update(las)
  return(las)
}
