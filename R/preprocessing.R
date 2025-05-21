#' Normalize `LAS` object representing segmented tree.
#'
#' Function to normalize LAS object. Function calculates
#' ground level based on the parameter specified by `quantile`,
#' subtracts it from all `Z`.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param quantile Z quantile at which grown level is specified since
#' ground points may not be identifiable with common algorithms if
#' ground points are removed during segmentation
#' @examples
#' library(lidR)
#' las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
#' # view histogram of Z values ranging from  -18 to -7 m
#' hist(las$Z)
#' las = normalize_las(las)
#' # view histogram of Z values now ranging from  0 to 11 m
#' hist(las$Z)
#' @importFrom lidR las_update filter_poi
#' @importFrom stats quantile
#' @import data.table
#' @export
normalize_las = function(las, quantile=c(0.001)) {
  #normalize
  ground_level = stats::quantile(las$Z, quantile)
  las@data[, Z := las$Z - ground_level]
  #update LAS header
  las = lidR::las_update(las)
  return(las)
}

#' Recenter `LAS` object representing segmented tree based on the bole
#' location
#'
#' Function calculates the tree location using points below specified
#' `height` and recenters on `X=0 Y=0`
#' @param las `LAS` object from `lidR` package representing
#' @param height consider only points where Z < height, if specified.
#' Useful for considering only the tree bole, for centering.
#' individually segmented tree. Set `height = NULL` to recenter
#' using all points.
#' @examples
#' library(lidR)
#' las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
#' # view histogram of original X/Y values
#' par(mfrow=c(1,2))
#' hist(las$X)
#' hist(las$Y)
#' las = recenter_las(las)
#' # view histogram of X/Y values centered on 0,0
#' par(mfrow=c(1,2))
#' hist(las$X)
#' hist(las$Y)
#' @importFrom lidR las_update
#' @import data.table
#' @export
recenter_las = function(las, height = 1) {
  if(is.null(height)) {
    centroid = apply(las@data[, c('X', 'Y')],2,mean)
  } else {
    bole = lidR::filter_poi(las, Z < height)
    centroid = apply(bole@data[, c('X', 'Y')],2,mean)
  }
  x_offset = centroid[1]
  y_offset = centroid[2]
  las@data[, X := las$X - x_offset]
  las@data[, Y := las$Y - y_offset]
  las = lidR::las_update(las)
  return(las)
}

#' Load, Recenter, and Remove low vegetation from `LAS` object representing
#' segmented tree
#'
#' Function to normalize, remove noise, remove vegetation, and recenter `LAS`
#' representing segmented tree. Vegetation cleaning is accomplished by
#' identifying stem points (TreeLS::StemPoints) and removing all but the Stem
#' below the `z.treshold`.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param Z.threshold numeric, height threshold below which all stem points
#' can be considered vegetation.
#' @param quantile See `normalize_las`. Z quantile at which grown level is specified since
#' ground points may not be identifiable with common algorithms if
#' ground points are removed during segmentation#'
#' @param height See `recenter_las`. Consider only points where Z < height, if specified.
#' Useful for considering only the tree bole, for centering.
#' individually segmented tree. Set `height = NULL` to recenter
#' using all points.
#' @examples
#' library(lidR)
#' las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
#' las_cleaned = clean_las(las)
#' plot(las)
#' plot(las_cleaned)
#' @importFrom lidR filter_duplicates classify_noise ivf
#' @importFrom TreeLS stemPoints
#' @import data.table
#' @export
clean_las = function(las, Z.threshold = 1, bole_height=1, quantile=0.001) {
  las = normalize_las(las)
  #identify points that are part of the stem and remove them
  las = suppressMessages(TreeLS::stemPoints(las))
  las = lidR::filter_poi(las, Z >  Z.threshold | (Z <  Z.threshold & Stem == TRUE))
  #recenter on the bole.
  las = recenter_las(las)
  las = lidR::filter_duplicates(las)
  las = lidR::classify_noise(las, lidR::ivf(res=0.1,n=3))
  return(las)
}


# Function to set path of Tree qsm to use in other functions
set_TREEQSM_PATH = function(path){
  if(!dir.exists(path)) stop('invalid path')
  assign(TREEQSM_PATH, path, envir = .GlobalEnv)
}

# Function to set path of  Matlab to use in other functions
set_MATLAB = function(path){
  if(!dir.exists(path)) stop('invalid path')
  assign(MATLAB_PATH, path, envir = .GlobalEnv)
}


# Function to generate TreeQSM matlab script with specified parameters
treeQSM = function() {stop('matlab download is down')}

#' Rotate `LAS` object about the `Z` axis
#'
#' Rotate `LAS` object about the `Z` axis for specified angle.
#' @param las `LAS` object from `lidR` package representing
#' individually segmented tree
#' @param angle numeric - in degrees, rotation angle about Z axis.
#' @examples
#' library(lidR)
#' las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
#' plot(las)
#' las = rotate_las_z(las, 90)
#' plot(las)
#' @importFrom recexcavAAR rotate
#' @export
rotate_las_z = function(las, angle) {
  pc = las@data[,c('X','Y','Z')]
  pc = recexcavAAR::rotate(pc$X, pc$Y, pc$Z, degrx = 0, degry = 0, degrz = angle)
  return(pc)
}


