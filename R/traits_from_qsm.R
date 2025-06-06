#' Get center of mass from tree bole segments QSM (X, Y, Z)
#'
#' This function calculates the center of mass/volume from a QSM
#' by estimating the centroid of cylinder locations, each weighted
#' by their volume. Only trunk sections are included (e.g., `branching_order == 0`).
#' For center of mass, assumes constant density
#' within segments.
#' @param qsm qsm object loaded from `[load_qsm]`.
#' @importFrom dplyr mutate summarize
#' @examples
#' # example code
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' print(get_center_of_mass)
#' @export
get_center_of_mass = function(qsm) {
  qsm = dplyr::filter(qsm, branching_order == 0)
  qsm = dplyr::mutate(qsm,
                      X.mid = (startX + endX)/2,
                      Y.mid = (startY + endY)/2,
                      Z.mid = (startZ + endZ)/2)
  centroid = dplyr::summarize(qsm,
                              X = weighted.mean(X.mid, volume),
                              Y = weighted.mean(Y.mid, volume),
                              Z = weighted.mean(Z.mid, volume))
  return(centroid)
}


#' Get internode distances between primary branches from a QSM
#'
#' This function estimates the internode distance between primary branches
#' from a QSM. It filters out all primary branches `branching_order == 1`
#' calculates their attachment points (Z) to the trunk, and then
#' returns the distances betweent the branches
#' @param qsm a QSM loaded using `[load_qsm()]`.
#' @param min_diam numeric - minimum diameter (in cm) to include branch
#' @importFrom dplyr filter mutate arrange
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' inodes = internode_distances(qsm)
#' print(inodes)
#' median(inodes)
#' @export
internode_distances = function(qsm, min_diam = 2) {
  x = get_primary_branches(qsm)
  x = dplyr::arrange(x, ht_m)
  x = dplyr::filter(x, diam_cm >= min_diam)
  x = dplyr::mutate(x, diff = NA)
  if(nrow(x) < 2) return(NA) #requires 2 branches to calculate
  for(i in 2:nrow(x)) {x$diff[i] = x$ht_m[i]-x$ht_m[i-1]}
  return(x$diff[-1])
}


#' Get Branch size distribution from a QSM
#'
#' This function outputs the volume (in mL) distribution of branches across
#' different branch diameter classes (in cm). Outputs a table that can be used
#' in other functions to find branch_skewness or branch_volume_weighted_stats()
#' @param qsm -- a QSM loaded using `[load_qsm()]`.
#' @param breaks numeric -- a vector of diameter classes (in cm) by which to
#'   summarize branch volume. If `NULL` the branch of branch sizes will be
#'   distributed across 1 cm bins.
#' @param plot boolean -- indicates whether the branch diameter distribution
#'   should be plotted as a histogram.
#' @importFrom dplyr summarize group_by
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' branch_distribution = branch_size_distribution(qsm, plot=TRUE)
#' print(branch_distribution)
#'
#' #volume-weighted mean
#' branch_volume_weighted_stats(qsm, FUN = function(x) mean(x))
#'
#' #volume-weighted median
#' branch_volume_weighted_stats(qsm, FUN = function(x) median(x))
#'
#' # volume-weighted skewness
#' branch_volume_weighted_stats(qsm, FUN = function(x) 3*(mean(x) - median(x)) / sd(x)))
#' @export
branch_size_distribution = function(qsm, breaks = NULL, plot=TRUE) {
  # Summarize primary branch attachment points
  qsm = dplyr::filter(qsm, branching_order > 0)
  if(nrow(qsm) < 2) return(NA) # need multiple branches to get a distribution
  qsm$volume_mL = qsm$volume * 100*100*100
  qsm$diam_cm = qsm$radius_cyl*200
  qsm$midpoint = floor(qsm$diam_cm) + 0.5
  if(is.null(breaks)) {breaks = seq(0, floor(max(qsm$diam_cm)) + 1, by = 1)}
  qsm$diameter_cm = cut(qsm$diam_cm, breaks= breaks)
  qsm = dplyr::summarize(dplyr::group_by(qsm,diameter_cm), midpoint = midpoint[1], volume_mL = sum(volume_mL))
  if(plot) barplot(volume_mL ~ midpoint, data=qsm, x='Branch diameter (cm)', y='Total volume (mL)')
  return(qsm)
}


#' Calculate volume weighted branch statistics
#'
#' This function calculates statistics on branch diameters weighted by the
#' volume of branches of that size based on outputs from
#'  `branch_size_distribution()`. The user defined function `FUN` can take
#'  any form of f(x) where x is a vector of diameters of length 1 for every
#'  mL of volume of that branch size class. See Details for recommended values
#'  for `FUN`.
#'  #Details
#'  Values of central tendency are recommended, but not variance since the
#'  weighted means are simulated.
#'
#'  Recommended values of `FUN` are:
#'
#'  Mean
#'  FUN = function(x) mean(x)
#'
#'  Median
#'  FUN = function(x) median(x)
#'
#'  Skewness
#'  FUN = function(x) 3*(mean(x) - median(x)) / sd(x)
#' @param qsm a QSM loaded using `[load_qsm()]`.
#' @param breaks numeric -- a vector of diameter classes (in cm) by which to
#' summarize branch volume. If `NULL` the branch of branch sizes will be
#' distributed across 1 cm bins.
#' @param FUN function -- central tendency function to be weighted based on
#' branch volume.
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' branch_distribution = branch_size_distribution(qsm, plot=TRUE)
#' print(branch_distribution)
#'
#' #volume-weighted mean
#' branch_volume_weighted_stats(qsm, FUN = function(x) mean(x))
#'
#' #volume-weighted median
#' branch_volume_weighted_stats(qsm, FUN = function(x) median(x))
#'
#' # volume-weighted skewness
#' branch_volume_weighted_stats(qsm, FUN = function(x) 3*(mean(x) - median(x)) / sd(x))
#' @export
branch_volume_weighted_stats = function(qsm, breaks=NULL, FUN = function(x) mean(x)) {
  stopifnot(is.function(FUN))
  x = branch_size_distribution(qsm, plot=FALSE, breaks=breaks)
  x = lapply(1:nrow(x), function(i) rep(x$midpoint[i], round(x$volume_mL[i])))
  x = do.call('c', x)
  return(FUN(x))
}

#' Extract primary branches from a QSM
#'
#' Extract primary branches from a QSM by filtering out the the trunk,
#' and identifying all cylinders where `branching_order == 1` and
#' are attached to the trunk. Returns a tibble containing the
#' basal diameter and height of attachment points.
#' @param qsm a QSM loaded using `[load_qsm()]`.
#' @importFrom dplyr filter reframe
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' primary_branches = get_primary_branches(qsm)
#'
#' #number of primary branches
#' nrow(primary_branches)
#' @export
get_primary_branches = function(qsm)
{
  # Summarize primary branch attachment points
  trunk_qsm = dplyr::filter(qsm, branching_order == 0)
  branch_qsm = dplyr::filter(qsm, branching_order == 1)
  primary_branches = dplyr::filter(branch_qsm, parent_ID %in% trunk_qsm$segment_ID)
  branches = dplyr::reframe(primary_branches,
                            section = 'branches',
                            diam_cm = radius_cyl*200,
                            ht_m = startZ*1,
                            volume = NA)
  return(branches)
}


#' Volume distribution from QSM
#'
#' This function estimates tree volume and its vertical distribution
#' from a QSM. The function separates the QSM into (1) trunk sections
#' (2) terminus (top of trunk < 4 cm dbh), and (3) primary branches.
#' The function divides trunk into segments defined by `segment_size`, calculates QSM volume,
#' For tree portions identified as branches the function only returns
#' the diameter. Both of these can be used in mass-volume equations
#' as needed.
#' @param qsm a QSM loaded using `[load_qsm()]`.
#' @param terminus_diam_cm numeric - trunk diameter at which it is treated as a branch.
#' @param segment_size numeric length of trunk segments in which to summarize volume.
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' volume = qsm_volume(qsm)
#' print(volume)
#' @importFrom dplyr filter group_by summarize reframe
#' @export
qsm_volume_distribution = function(qsm, terminus_diam_cm = 4, segment_size=0.5) {
  # Identify Trunk sections
  trunk_qsm = dplyr::filter(qsm, branching_order == 0)
  terminal_branch.ht = max(dplyr::filter(trunk_qsm, radius_cyl > terminus_diam_cm/200)$startZ)
  trunk = dplyr::filter(trunk_qsm, startZ <= terminal_branch.ht)
  #Divide trunk into  sections/chunks and estimate mass of each.
  trunk$ht = trunk$startZ - trunk$startZ[1]
  trunk$midpt = apply(trunk[,c('startZ', 'endZ')],1, mean)
  section_hts = seq(0,max(trunk$ht),by=segment_size)
  if(max(trunk$ht) > max(section_hts)) section_hts = c(section_hts, tail(section_hts,1)+segment_size)
  trunk$section = as.numeric(cut(trunk$ht, section_hts, right=FALSE))
  trunk = dplyr::summarize(dplyr::group_by(trunk, section),
                           diam_cm = weighted.mean(radius_cyl, volume)*200,
                           ht_m = weighted.mean(midpt, volume),
                           volume = sum(volume))
  trunk$section = 'trunk'

  #Summarize terminal section
  terminus = dplyr::filter(trunk_qsm, startZ > terminal_branch.ht)
  terminus$midpt = apply(terminus[,c('startZ', 'endZ')],1, mean)
  terminus = dplyr::summarize(terminus, section = 'terminus',
                              diam_cm = radius_cyl[which.min(startZ)]*200,
                              ht_m = weighted.mean(midpt, volume),
                              volume = sum(volume))


  # Summarize primary branch attachment points
  branch_qsm = dplyr::filter(qsm, branching_order == 1)
  primary_branches = dplyr::filter(branch_qsm, parent_ID %in% trunk_qsm$segment_ID)
  branches = dplyr::reframe(primary_branches,
                            section = 'branches',
                            diam_cm = radius_cyl*200,
                            ht_m = startZ*1,
                            volume = NA)

  output = rbind(trunk, terminus, branches)
  return(output)
}

#' Fit taper equation to QSM
#'
#' This function fits a taper equation to trunk sections of a QSM using
#' Kozak the Model (2002, 2007).
#'
#' $ d(h)/D = a0 (h/H) + a1 (h/H) + a2 (h/H)^2 + a3 (h/H)^3 $
#'
#' The function groups QSM cylinders into segments of `segment_size` up to
#' `terminus_diam` which is the maximum diameter at which the taper equation
#' ends.
#' @param qsm a QSM loaded using `[load_qsm()]`.
#' @param dbh numeric -- required to fit Kozak model, not calculated from QSM, so
#' as not to conflict with other more accurate means of measurement e.g., `get_DBH'
#' @param terminus_diam_cm numeric -- the trunk diameter at which is no longer
#' considered trunk
#' @param segment_size numeric -- the length of segments that QSM cylinders are
#' grouped into
#' @param plot boolean -- indicates whether model output should be plotted. Plots
#' are found in the output list as object$plot, regardless of this setting.
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot aes geom_point theme_bw labs
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' fit_taper(qsm, dbh = 13.8)
#' @export
#'
fit_taper_Kozak = function(qsm, dbh, terminus_diam_cm = 4, segment_size=0.25, plot=TRUE) {
  taper = qsm_volume_distribution(qsm, terminus_diam_cm = 2, segment_size=0.25)
  taper = dplyr::select(dplyr::filter(taper, section == 'trunk'), ht_m, diam_cm)
  H = max(qsm$endZ)
  # Kozak function
  d = diam_cm ~ (a0 + a1*(ht_m/H) + a2*(ht_m/H)^2 + a3*(ht_m/H)^3)/dbh
  f = function(a0, a1, a2, a3, ht_m) (a0 + a1*(ht_m/H) + a2*(ht_m/H)^2 + a3*(ht_m/H)^3)/dbh
  mod = nls(d, data=taper, start = list(a0=1, a1=-1.3, a2=3, a3=-5))
  r2 = cor(taper$diam_cm, predict(mod))^2
  rmse = sqrt(mean(resid(mod)^2))
  results = as.data.frame(t(coef(mod)))
  results$r2 = r2
  results$rmse = rmse
  Hs = seq(0,H,by=0.1)
  modfit = data.frame(ht_m = Hs,diam_cm=f(results$a0, results$a1, results$a2, results$a3, Hs))
  myPlot = suppressWarnings(ggplot2::ggplot(taper, ggplot2::aes(x=ht_m, y=diam_cm)) + ggplot2::geom_point() +
    ggplot2::geom_line(data=modfit) + ggplot2::lims(y=c(0,max(modfit$diam_cm)), x=c(0,H)) +
    ggplot2::theme_bw() + ggplot2::labs(x='Height (m)', y='Diameter (cm)') + theme(panel.grid = ggplot2::element_blank()))
  if(plot) suppressWarnings(print(myPlot))
  output = list(data=taper, plot=myPlot, results=results)
  return(output)
}

#' Horizontal offset of center of mass from QSM
#'
#' This function extracts the horizontal distance between the base of a tree and
#' the center of a tree from a QSM. The function takes the coordinate of the
#' lowest QSM segment, and the center of mass, and finds the horizontal
#' distance between them.
#' @param qsm a QSM loaded using `[load_qsm()]`.
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' print(get_center_of_mass(qsm))
#' print(get_com_offset(qsm))
#' @export
get_com_offset = function(qsm) {
  #get XY of lowest segmented
  coord_base = qsm[which.min(qsm$startZ),][, c('startX', 'startY', 'startZ')]
  coord_com = get_center_of_mass(qsm)
  #return horizontal distance between points
  coords = matrix(c(coord_base[1,1:2], coord_com[1,1:2]), nrow=2, byrow = TRUE)
  offset = as.numeric(dist(coords))
  return(offset)
}


#' Calculate tree sweep from straight line from QMS
#'
#' This function calculates tree sweep from a QSM. Starting with an idealized
#' vector of a straight tree (straight line from top to bottom QSM segment)
#' the function caclucates devations of points along the trunk from the
#' idealized vector. It resturns sweep from each QSM segment so that
#' summary statistics can be computed by the user. See also `get_stem_deflection()`
#' @param qsm QSM loaded using `[load_qsm()]`.
#' @param terminus_diam_cm numeric -- the trunk diameter at which is no longer
#' considered trunk
#' @param plot boolean -- indicates whether graph of sweep should be plotted.
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' print(get_center_of_mass(qsm))
#' print(get_com_offset(qsm))
get_stem_sweep = function(qsm, terminus_diam_cm = 4, plot=TRUE) {
  bole = dplyr::filter(qsm, branching_order==0 & radius_cyl > terminus_diam_cm/200)

  # Setup endpoints for tilt and sweep reference line
  endpoint_segments = c(which.min(bole$startZ), which.max(bole$endZ))
  bottom = as.numeric(bole[endpoint_segments[1], c('startZ','startY', 'startZ')])
  top = as.numeric(bole[endpoint_segments[2], c('endX','endY', 'endZ')])
  endpoints = matrix(c(bottom, top), nrow=2, byrow=TRUE)
  colnames(endpoints) = c('X','Y','Z')

  # Calculate sweep
  trunk_pts = dplyr::mutate(bole, X = (startX + endX)/2,
                            Y = (startY + endY)/2,
                            Z = (startZ + endZ)/2, .keep='none')

  # function to get vector length
  length = function(x) sqrt(sum(x^2))
  # calculate projection of trunk point onto idealized vector, and get distance
  # between them "sweep"
  #https://www.youtube.com/watch?v=_XcDo2FNRe0

  sweep_vector = function(p) {
    a = apply(endpoints, 2, diff) #idealized vector
    b = p - endpoints[1,] #vector from tree base to point p
    proj = sum(a*b)/length(a) #projection of b onto a
    sweep = sqrt(length(b)^2 - proj^2) #distance from projected line to trunk point
    return(sweep)
  }
  sweep = apply(trunk_pts, 1, sweep_vector)
  sweep = data.frame(Height = as.numeric(trunk_pts$Z), sweep = sweep)
  if(plot) plot(sweep ~ Height, data = sweep, type='l', xlab = 'Sweep (m)', ylab = 'Height (m)')
  return(sweep)
}


#' Get tree tilt from QSM
#'
#' This function calculates tilt of a tree from a QSM. The function identifies
#' the upper and lower extreme segments of the QMS (trunk sections only) and
#' computes a vector between them, and returns the devation of that angle
#' from directly vertical.
#' @param qsm a QSM loaded using `[load_qsm()]`.
#' @param terminus_diam_cm numeric - trunk diameter at which it is treated as a branch.
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' get_stem_tilt(qsm)
#' @export
get_stem_tilt = function(qsm, terminus_diam_cm = 4) {
  #extract tree bole
  bole = dplyr::filter(qsm, branching_order==0 & radius_cyl > terminus_diam_cm/200)
  # Setup endpoints for tilt and sweep reference line
  endpoint_segments = c(which.min(bole$startZ), which.max(bole$endZ))
  bottom = as.numeric(bole[endpoint_segments[1], c('startZ','startY', 'startZ')])
  top = as.numeric(bole[endpoint_segments[2], c('endX','endY', 'endZ')])
  endpoints = matrix(c(bottom, top), nrow=2, byrow=TRUE)
  colnames(endpoints) = c('X','Y','Z')
  #get straight line tilt
  tilt = as.vector(acos(diff(endpoints[,'Z'])/ dist(endpoints)) * 180/pi)
  return(tilt)
}
