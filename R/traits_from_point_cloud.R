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
get_height = function(las, quantiles = c(0, 1)) {
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
#' of crown isolated from larger continuous crown.
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
#' caluclating total area, but also used to detect crown base
#' height.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree, with the crown labeled.
#' @param segment_height numeric - height of each segment in which to calculate area
#' @param quantile numeric - quantile at which width is measured
#' Values in the interval approaching 0 (e.g., 0.001) are recommended to
#' trim random noise
#' @param angle numeric - angle at which to rotate the point cloud prior
#' to estimating area. Useful in a loop if quantifying mulitple angles
#' @importFrom dplyr filter
#' @importFrom stats quantile
#' @importFrom tibble tibble
get_area_profile = function(las, segment_height=0.25, quantile = c(0.001), angle = 0) {
  pc = rotate_las_z(las, angle)@data[, c('X', 'Y', 'Z')]
  heights = quantile(las$Z, probs = c(quantile, 1-quantile))
  myfun = function(slice_min) {
    slice_max = slice_min + segment_height
    slice = dplyr::filter(pc, Z > slice_min & Z <= slice_max)
    width = diff(stats::quantile(slice$X, probs = c(quantile,1-quantile)))
    area = width * segment_height
    return(tibble::tibble(bottom = slice_min, top = slice_max, width = width, area = area))
  }
  output = do.call(rbind, lapply(seq(heights[1],heights[2],by=segment_height), function(i) myfun(i)))
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
#' @param crown_base_height numeric - height of crown base for segmentation.
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
#' @examples
#' las = lidR::readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
#' las = clean_las(las)
#' cbh = get_crown_base(las, threshold=0.25, sustain=2)
#' las = segment_crown(las, cbh)
#' get_crown_volume_voxel(las)
#' get_crown_volume_alpha(las)
#' st_area(convex_hull_2D(las)) #profile area, convex hull
#' st_area(voxel_hull_2D(las)) #profile area, voxel hull
#' get_lacunarity(las)
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
#' @examples
#' las = lidR::readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
#' cbh = get_crown_base(las, threshold=0.25, sustain=2)
#' las = segment_crown(las, cbh)
#' get_crown_volume_voxel(las)
#' get_crown_volume_alpha(las)
#' st_area(convex_hull_2D(las)) #profile area, convex hull
#' st_area(voxel_hull_2D(las)) #profile area, voxel hull
#' get_lacunarity(las)
get_crown_volume_alpha = function(las, resolution = 0.1, alpha=0.5) {
  if(!'Crown' %in% colnames(las@data)) {
    stop('las does not contain column called `Crown` use `segment_crown()`')
  }
  crown = lidR::filter_poi(las, Crown == 1)
  vox = lidR::voxelize_points(crown, res = resolution)
  vol = suppressWarnings(ITSMe::alpha_volume_pc(vox@data[, c('X','Y','Z')], alpha=alpha))
  return(c(crown_volume_alpha = vol$av))
}


#' Returns the convex hull representing vertical crown area
#'
#' This function generates an `sf` object representing th vertical crown area
#' of a `LAS` object based using the convex hull of a 2D vertical projection.
#' @param las `LAS` object from `lidR` package representing
#' the CROWN of a tree. Crowns must be segmented using [segment_crown()].
#' @param angle numeric - in degrees, rotation angle about Z axis.
#' @importFrom lidR filter_poi
#' @importFrom sf st_as_sf st_convex_hull st_union
#' @export
#' @examples
#' las = lidR::readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
#' cbh = get_crown_base(las, threshold=0.25, sustain=2)
#' las = segment_crown(las, cbh)
#' get_crown_volume_voxel(las)
#' get_crown_volume_alpha(las)
#' st_area(convex_hull_2D(las)) #profile area, convex hull
#' st_area(voxel_hull_2D(las)) #profile area, voxel hull
#' get_lacunarity(las)
convex_hull_2D = function(las, angle = 0) {
  if(!'Crown' %in% colnames(las@data)) {
    stop('las does not contain column called `Crown` use `segment_crown()`')
  }
  #vertical projection on X-Z plane
  las = lidR::filter_poi(las, Crown == 1)
  las = rotate_las_z(las, angle)
  las = las@data[,c('X', 'Z')]
  las = sf::st_as_sf(las, coords = c('X', 'Z'))
  las = sf::st_convex_hull(sf::st_union(las))
  las = sf::st_as_sf(las)
  return(las)
}


#' Returns an `sf`representing the vertical crown area from voxelization
#'
#' This function generates an `sf` object representing th vertical crown area
#' of a `LAS` object by voxelizing a 2D vertical projection.
#' @param las `LAS` object from `lidR` package representing
#' the CROWN of a tree. Crowns must be segmented using [segment_crown()].
#' @param resolution numeric - resolution of voxelization
#' @param angle numeric - in degrees, rotation angle about Z axis.
#' @importFrom lidR filter_poi voxelize_points rasterize_density
#' @importFrom sf st_as_sf st_as_sf st_crs
#' @importFrom terra as.polygons
#' @export
#' @examples
#' las = lidR::readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
#' cbh = get_crown_base(las, threshold=0.25, sustain=2)
#' las = segment_crown(las, cbh)
#' get_crown_volume_voxel(las)
#' get_crown_volume_alpha(las)
#' st_area(convex_hull_2D(las)) #profile area, convex hull
#' st_area(voxel_hull_2D(las)) #profile area, voxel hull
#' get_lacunarity(las)
voxel_hull_2D = function(las, resolution = 0.1, angle = 0) {
  if(!'Crown' %in% colnames(las@data)) {
    stop('las does not contain column called `Crown` use `segment_crown()`')
  }
  crown = filter_poi(las, Crown == 1)
  if(angle != 0) las = rotate_las_z(crown, angle)
  crown@data[, Y := crown$Z]
  flat_las = lidR::voxelize_points(crown, resolution)
  proj = lidR::rasterize_density(flat_las, res = resolution)
  proj = proj > 0
  proj[proj == 0] = NA
  proj = sf::st_as_sf(terra::as.polygons(proj))
  sf::st_crs(proj) = NA
  return(proj)
}


#' Calculate crown lacunarity from a tree crown
#'
#' This function calculates the lacunarity or "porosity" of a tree crown
#' by comparing the ratio of a voxelized crown hull and a convex hull.
#' See `voxel_hull_2D()` and `convex_hull_2D()`
#' @param las `LAS` object from `lidR` package representing
#' the CROWN of a tree. Crowns must be segmented using [segment_crown()].
#' @param resolution numeric - resolution of voxelization
#' @param angle numeric - in degrees, rotation angle about Z axis.
#' @examples
#' las = lidR::readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
#' cbh = get_crown_base(las, threshold=0.25, sustain=2)
#' las = segment_crown(las, cbh)
#' get_crown_volume_voxel(las)
#' get_crown_volume_alpha(las)
#' st_area(convex_hull_2D(las)) #profile area, convex hull
#' st_area(voxel_hull_2D(las)) #profile area, voxel hull
#' get_lacunarity(las)
get_lacunarity = function(las, res = 0.1, angle = 0) {
  # Make a voxel hull, and get the ratio of its area relative to same
  # hull with holes filled.
  if(!'Crown' %in% colnames(las@data)) stop('las does not contain column called `Crown` use `segment_crown()`')
  las = lidR::filter_poi(las, Crown == 1)
  if(angle != 0) las = rotate_las_z(las, angle)
  voxel_hull_area = sf::st_area(voxel_hull_2D(las, res = res))
  convex_hull_area = sf::st_area(convex_hull_2D(las))
  lacunarity = voxel_hull_area/convex_hull_area
  return(lacunarity)
}

#' Calculate crown leverage from point cloud
#'
#' This function calculates the lever arm of canopies. The function
#' #' is a simple wrapper for `get_area_profile()`. It calculates the crown area
#' in segments defined by `segement_height`, multiplies the area of each of those
#' segments by their height, and then returns the sum of all segments. This
#' is proportaionl to drag calculations on the tree assuming  windspeed is
#' invariant with height.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree, with the crown labeled. See `segment_crown()`
#' @param segment_height numeric - height of each segment in which to calculate area
#' @param quantile numeric - quantile at which width is measured
#' Values in the interval approaching 0 (e.g., 0.001) are recommended to
#' trim random noise
#' @param angle numeric - angle at which to rotate the point cloud prior
#'
#' to estimating area. Useful in a loop if quantifying mulitple angles
#' @importFrom dplyr filter
get_crown_lever_arm = function(las, segment_height=0.25, quantile = c(0.001), angle=0) {
  if(!'Crown' %in% colnames(las@data)) {
    stop('las does not contain column called `Crown` use `segment_crown()`')
  }
  las = lidR::filter_poi(las, Crown == 1)
  profile = get_area_profile(las, segment_height=segment_height, quantile=quantile, angle=angle)
  profile = dplyr::mutate(profile, midpt = (bottom+top)/2, lever = area*midpt)
  lever = sum(profile$lever, na.rm=TRUE)
  return(lever)
}
