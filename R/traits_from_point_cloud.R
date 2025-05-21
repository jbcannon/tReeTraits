#' Extract height from  `LAS` object representing segmented tree.
#'
#' Function to extract height from LAS object. Function calculates
#' difference between two specified quantiles from the `Z` attribute.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param quantiles Z quantiles at which ground level and highest point
#' are measured. Values in the interval (0,1) are recommended to trim
#' random noise.
#' @examples
#' library(lidR)
#' las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
#' las = clean_las(las)
#' print(get_height(las))
#' @importFrom stats quantile
#' @export
get_height = function(las, quantiles = c(0.001, 0.999)) {
  stopifnot(class(las) == 'LAS')
  stopifnot(length(quantiles) == 2)
  stopifnot(quantiles >= 0 & quantiles <= 1)
  stopifnot(quantiles[1]<quantiles[2])
  x = stats::quantile(las$Z, quantiles)
  x = as.numeric(diff(x))
  return(c(height=x))
}

#' Extract width from  `LAS` object representing segmented tree.
#'
#' Function to extract width from LAS object. Function calculates
#' difference between two specified quantiles from the `X` and `Y` attributes
#' and returns both widths and their average.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param quantiles Z quantiles at which widths are measured
#' are measured. Values in the interval (0,1) are recommended to trim
#' random noise.
#' @examples
#' library(lidR)
#' las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
#' las = clean_las(las)
#' print(get_width(las))
#' @importFrom stats quantile
#' @export
get_width = function(las, quantiles = c(0.001, 0.999)) {
  x = as.numeric(diff(stats::quantile(las$X, quantiles)))
  y = as.numeric(diff(stats::quantile(las$Y, quantiles)))
  wd = mean(c(x,y))
  return(c(mean_width = wd, x_width = x, y_width=y))
}

#' Extract diameter at breast height from  `LAS` object representing segmented tree.
#'
#' Function to extract diameter at breast height (1.37 m) from LAS object.
#' Function filters LAS keeping only points with Intensity greater than
#' specified threshold. Function calculates verticality eigenvalue and filters
#' based on verticality threshold. Last, diameter is calculated using a
#' RANSAC cylinder fitting algorithm.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param intensity_threshold numeric - filter value for Intensity to help remove vegetation
#' @param verticality_threshold numeric - filter value for Verticality threshold to remove
#' horizontal branches.
#' @param select_n numeric - number of points selected on every RANSAC iteration.
#' @examples
#' library(lidR)
#' las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
#' las = clean_las(las)
#' print(get_dbh(las))
#' @importFrom lidR filter_poi add_lasattribute
#' @importFrom spanner eigen_metrics
#' @importFrom TreeLS shapeFit
#' @importFrom magrittr `%>%`
#' @export
get_dbh = function(las, intensity_threshold=41000,
                   select_n = 10, verticality_threshold=0.9) {
  if(is.null(intensity_threshold)) {
    warning('`get_intensity_treshold() not yet implemented, using default value 41000')
    intensity_treshold = 41000
    #message('No intensity threshold specfied. Estimating itensity using get_intensity_threshold')
    #intensity_threshold = get_intensity_threshold(las)$threshold
  }
  stopifnot(!is.null(verticality_threshold))
  stopifnot(verticality_threshold >= 0 & verticality_threshold < 1)
  bole = lidR::filter_poi(las, Intensity > intensity_threshold & Z < 3)
  #filter to speed up processing.
  bole = lidR::decimate_points(bole, lidR::random_per_voxel(res=0.025))
  eigen = spanner::eigen_metrics(bole)$Verticality
  bole = lidR::add_lasattribute(bole, eigen, 'Verticality', 'Verticality')
  bole = lidR::filter_poi(bole, Z> 1 & Z < 2 & Verticality > verticality_threshold)
  cyl = spanner::cylinderFit(bole, method='ransac', n = select_n, n_best = 3)
  dbh = cyl$radius*2
  return(c(dbh = as.numeric(dbh)))
}

#' Estimate Crown base height of `LAS` representing segmented tree.
#'
#' This function estimates the crown base height by analyzing the vertical
#' profile of the tree using [get_area_profile()] which breaks the profile
#'into segments of height `segment_height`. The function estimates
#' segments exceed a threshold specified by `threshold` which must be exceeded
#' `sustain` times.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param threshold numeric - threshold width at which crown becomes apparent.
#' Recommend a width ~2X greater than anticipated DBH.
#' @param sustain numeric - number of segments in a row that treshold must be
#' exceeded before identifying start of crown. This is to exclude small segments
#' of crown isolated from larger continuous canopy.
#' @param segment_height numeric - height of each segment in which to calculate area
#' @param quantile numeric - quantile at which width is measured
#' Values in the interval approaching 0 (e.g., 0.001) are recommended to
#' trim random noise
#' @importFrom stats quantile
#' @export
get_crown_base = function(las, threshold = 0.5, sustain = 2,
                           segment_height = 0.25, quantile = 0.01) {
  area_profile = get_area_profile(las, segment_height=segment_height, quantile=quantile)
  x=area_profile$width
  area_profile$difference = c(sapply(2:length(x), function(i) x[i+1] - x[i]), NA)
  # Which segments have a width increase greater than threshold width?
  area_profile$wide = area_profile$difference > threshold | area_profile$width > 3*threshold
  #when is width sustained? i.e., greater than 2 segments
  x = which(area_profile$wide)
  d = sapply(2:length(x), function(i) x[i] - x[i-1])
  d = d == 1
  cbh = area_profile$bottom[x[d][1]]
  cbh = cbh - segment_height/2 # return the midpoint
  return(c(crown_base_height=cbh))
}

#' Generate area estimates of tree profile in segments
#'
#' This function calculates the area of the tree profile
#' by breaking it into segments of height `segment_height` and estimating
#' the width of each segement. Area profiles are useful for
#' caluclating total area, but also used to detect canopy base
#' height.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param segment_height numeric - height of each segment in which to calculate area
#' @param quantile numeric - quantile at which width is measured
#' Values in the interval approaching 0 (e.g., 0.001) are recommended to
#' trim random noise
#' @param rotation_angle numeric - angle at which to rotate the point cloud prior
#' to estimating area. Useful in a loop if quantifying mulitple angles
#' @importFrom dplyr filter
#' @importFrom stats quantile
#' @importFrom tibble tibble
#' @export
get_area_profile = function(las, segment_height=0.25, quantile = c(0.001), angle = 0) {
  pc = rotate_las_z(las, angle)
  height = get_height(las)
  ground = 0
  myfun = function(slice_min) {
    slice_max = slice_min + segment_height
    slice = dplyr::filter(pc, z > slice_min & z <= slice_max)
    width = diff(stats::quantile(slice$x, probs = c(quantile,1-quantile)))
    area = width * segment_height
    return(tibble::tibble(bottom = slice_min, top = slice_max, width = width, area = area))
  }
  output = do.call(rbind, lapply(seq(ground,height,by=segment_height), function(i) myfun(i)))
  output$angle = angle
  return(output)
}

#' Segment tree crown of `LAS` representing segmented tree.
#'
#' This function labels all points with $Z > `crown_base_height`$
#' and returns a labled LAS. If `crown_base_height` is not specified, it
#' is estimated with [get_crown_base()] using default parameters.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param crown_base_height numeric - height of canopy base for segmentation.
#' `NULL`, it is estimated with [get_crown_base()] using default parameters.
#' @importFrom lidR add_lasattribute
#' @export
segment_crown = function(las, crown_base_height = NULL) {
  if(is.null(crown_base_height)) {
    warning('no crown_base_height specified. Detected automatically with `get_crown_base()`')
    crown_base_height = get_crown_base(las)
  }
  crown = as.numeric(las$Z > crown_base_height)
  las = lidR::add_lasattribute(las, crown, 'Crown', 'Crown points')
  return(las)
}



#' Estimate crown volume by voxelization
#'
#' This function volume of a `LAS` object by thinning to a resolution specified
#' by `resolution`, and estimating volume using the equation
#' \deqn{Volume_{crown} = N_{occupied voxel} * Volume_{voxel}}
#' @param las `LAS` object from `lidR` package representing
#' the CROWN of a tree. Crowns can be segmented using [segment_crown()].
#' @param resolution numeric - resolution of voxelization
#' @importFrom lidR voxelize_points filter_poi
#' @export
get_crown_volume_voxel = function(las, resolution = 0.1) {
  if(!'Crown' %in% colnames(las@data)) stop('las does not contain column called `Crown` use `segment_crown()`')
  crown = lidR::filter_poi(las, Crown == 1)
  vox = lidR::voxelize_points(crown, res = resolution)
  volume = nrow(vox) * resolution^3
  return(c(crown_volume_vox = volume))
}

#' Estimate crown volume by  alpha shape volume
#'
#' This function volume of a `LAS` object by thinning to a resolution specified
#' by `resolution`, and estimating volume by fitting a alpha shape volume.
#' Crowns must be segmented using [segment_crown()].
#' @param resolution numeric - resolution of initial voxelization to increase speed
#' @param alpha numeric - alpha for the computation of the 3D alpha-shape of the point cloud.
#' See [ITSMe::alpha_volume_pc()].
#' @importFrom lidR voxelize_points filter_poi
#' @importFrom ITSMe alpha_volume_pc
#' @export
get_crown_volume_alpha = function(las, resolution = 0.1, alpha=0.5) {
  if(!'Crown' %in% colnames(las@data)) stop('las does not contain column called `Crown` use `segment_crown()`')
  crown = lidR::filter_poi(las, Crown == 1)
  vox = lidR::voxelize_points(crown, res = resolution)
  vol = suppressWarnings(ITSMe::alpha_volume_pc(vox@data[, c('X','Y','Z')], alpha=alpha))
  return(c(crown_volume_alpha = vol$av))
}
