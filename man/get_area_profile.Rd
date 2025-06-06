% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traits_from_point_cloud.R
\name{get_area_profile}
\alias{get_area_profile}
\title{Generate area estimates of tree profile in segments}
\usage{
get_area_profile(las, segment_height = 0.25, quantile = c(0.001), angle = 0)
}
\arguments{
\item{las}{`LAS` object from `lidR` package representing
individually segmented tree, with the crown labeled.}

\item{segment_height}{numeric - height of each segment in which to calculate area}

\item{quantile}{numeric - quantile at which width is measured
Values in the interval approaching 0 (e.g., 0.001) are recommended to
trim random noise}

\item{angle}{numeric - angle at which to rotate the point cloud prior
to estimating area. Useful in a loop if quantifying mulitple angles}
}
\description{
This function calculates the area of the tree profile
by breaking it into segments of height `segment_height` and estimating
the width of each segement. Area profiles are useful for
caluclating total area, but also used to detect crown base
height.
}
