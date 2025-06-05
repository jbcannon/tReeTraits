





#' Make 3-panel plot of tree point cloud to check for errors
#'
#' Plots 2 profiles X, Y, and and overhead Z view of a point
#' cloud to allow users to identify stray points, or errors in
#' segmentations.
#' @param las `LAS` object from `lidR` package representing
#' the CROWN of a tree. Crowns can be segmented using [segment_crown()].
#' @param res numeric - resolution of voxelization to speed up plotting
#' @examples
#' # example code
#' las = lidR::readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
#' las = clean_las(las)
#' plot_tree(las)
#' @importFrom lidR decimate_points random_per_voxel
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @export
plot_tree = function(las, res = 0.05, plot=TRUE) {
  las_thin = lidR::decimate_points(las, lidR::random_per_voxel(res=res))
  ggsettings = function(x) {
    x = x + ggplot2::geom_point(size=0.1)  +
      ggplot2::coord_fixed() +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     legend.position = 'none') +
      ggplot2::scale_color_viridis_c(direction = 1, option='G') +
      ggplot2::labs(y=ggplot2::element_blank(), x=ggplot2::element_blank()) +
      ggplot2::theme(axis.text = ggplot2::element_text(colour = '#FFFFFF'),
                     axis.ticks = ggplot2::element_line(color='#FFFFFF'))
    return(x)
  }
  v1 = ggsettings(ggplot2::ggplot(las_thin@data, ggplot2::aes(x=X,y=Z, color=Y)))
  v2 = ggsettings(ggplot2::ggplot(las_thin@data, ggplot2::aes(x=Y,y=Z, color=X)))
  #v3 = ggsettings(ggplot2::ggplot(dplyr::filter(las@data,Z<3), ggplot2::aes(x=X,y=Z, color=Y)))
  v4 = ggsettings(ggplot2::ggplot(las_thin@data, ggplot2::aes(x=X,y=Y, color=-Z)))
  profiles = ggpubr::ggarrange(v1, v2, v4, nrow=1, widths = c(0.3,0.3,0.4), labels=LETTERS)
  if(plot) plot(profiles)
  return(profiles)
}

#' Plot QSM in base R
#'
#' Simple function to create a diagnostic plot to view QSMs colored
#' by branching order.
#' @param qsm  a QSM loaded using `[load_qsm()]`.
#' @param scale a factor by which to multiply the `radius_cyl` column to
#' give line segments the appearance of volume
#' @examples
#' qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
#' qsm = load_qsm(qsm_file)
#' plot_qsm(qsm)
#' @export
plot_qsm = function(qsm, scale = 150, rotation=TRUE) {
  if(rotation) par(mfrow = c(1,2))
  par(mar=c(2,2,1,1)+0.5)
  ylim = range(qsm[, c('startZ', 'endZ')])
  plot(NA, NA, xlim = c(-1, 1), ylim = ylim ,asp=1, xlab='', ylab='', axes=FALSE)
  axis(1);axis(2);
  with(qsm, arrows(startX, startZ, endX, endZ, length=0, code=2, col=branching_order+1,lwd=radius_cyl*scale))
  if(rotation) {
    plot(NA, NA, xlim = c(-1, 1), ylim = ylim ,asp=1, xlab='', ylab='', axes=FALSE)
    axis(1);
    with(qsm, arrows(startY, startZ, endY, endZ, length=0, code=2, col=branching_order+1,lwd=radius_cyl*scale))
  }
  par(mfrow=c(1,1))
}


# Diagonistic plot to display basic tree  measurements
basics_diagnostic_plot = function(las, height, cbh, crown_width, dbh, res=0.1) {
  las_thin = lidR::decimate_points(las, lidR::random_per_voxel(res=res))
  marker_df = data.frame(
    x = c(crown_width/-2, crown_width/2, 0, 0, min(las$X)+2, min(las$X)+2, min(las$X)+3, min(las$X)+3, dbh/-2, dbh/2),
    y = c(height-(height-cbh)/2, height - (height-cbh)/2, min(las$Z), cbh, min(las$Z), height, min(las$Z), cbh, 1.5, 1.5),
    name = sort(rep(letters[1:5],2)),
    color = c(rep('red',4), rep('black', 4), rep('blue',2)))

  crown_met = ggplot2::ggplot() +
    ggplot2::geom_point(data = las_thin@data, ggplot2::aes(x=X, y=Z), size=0.05, color=grey(0.4))  +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() + ggplot2::geom_line(data = marker_df, ggplot2::aes(x=x,y=y,group=name, color=color),
                                             linewidth=1, arrow = ggplot2::arrow(length=ggplot2::unit(.1,"cm"), angle=90, ends="both")) +
    ggplot2::theme(legend.position='none') + ggplot2::labs(x=ggplot2::element_blank(), y=ggplot2::element_blank()) +
    ggplot2::theme(#plot.background = ggplot2::element_rect(fill = "black", colour = "black"),
      #panel.background = ggplot2::element_rect(fill = "black", colour = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = 'none')
  return(crown_met)
}

# Diagonistic plot to display generated hulls
hull_diagnostic_plot = function(las, res=0.1) {
  if(!'Crown' %in% colnames(las@data)) stop('las does not contain column called `Crown` use `segment_crown()`')
  crown = lidR::filter_poi(las, Crown == 1)
  convex_hull = convex_hull_2D(crown)
  voxel_hull = voxel_hull_2D(crown, res = res)
  x = ggplot2::ggplot() + ggplot2::geom_sf(data=voxel_hull, fill='chartreuse4') +
    ggplot2::geom_sf(data=convex_hull, fill=NA, linewidth=1, linetype='dashed', color=grey(0.2)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  return(x)
}


# Diagnostic plot to view taper function. Wrapper for fit_taper_Kozak which
# has a plot output.
taper_diagnostic_plot = function(qsm, dbh) {
  return(fit_taper_Kozak(qsm, dbh, plot=FALSE)$plot)
}

# diagnostic plot to view branch diameter distribution. Simpe wrapper
# for `branch_size_distribution` which does most of the heavy lifting.
branch_distribution_plot = function(qsm) {
  branches = branch_size_distribution(qsm, plot=FALSE)
  myPlot = ggplot2::ggplot(branches, ggplot2::aes(x=midpoint, y=volume_mL)) +
    ggplot2::geom_col() + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = 'Branch diameter (cm)', y = 'Total volume (mL)')
  return(myPlot)
}


#' Generate a diagnostic plot to assess basic metrics and QSM output
#'
#' @param las `LAS` object from `lidR` package representing
#' the CROWN of a tree. Crowns can be segmented using [segment_crown()]
#' @param qsm  a QSM loaded using `[load_qsm()]`.
#' @param res numeric - resolution of voxelization to speed up plotting
#' @export
full_diagnostic_plot = function(las, qsm, height, cbh, crown_width, dbh, res=0.1) {
  r = plot_tree(las)
  s = basics_diagnostic_plot(las, height, cbh, crown_width, dbh, res)
  t = hull_diagnostic_plot(las, res=res)
  u = ggplotify::as.ggplot(~plot_qsm(qsm, rotation=FALSE, scale=50))
  v = taper_diagnostic_plot(qsm, dbh)
  w = branch_distribution_plot(qsm)
  x = ggpubr::ggarrange(r, s, t, widths = c(0.5, 0.2, 0.3), ncol=3, labels=c(NA, 'D', 'E'))
  y = ggpubr::ggarrange(u, v, w, ncol=3, widths= c(0.2, 0.4,0.4), labels=c('F', 'G', 'H'))
  z = ggpubr::ggarrange(x, y, nrow=2)
  return(z)
}
