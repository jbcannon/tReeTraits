% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{plot_qsm}
\alias{plot_qsm}
\title{Plot QSM in base R}
\usage{
plot_qsm(qsm, scale = 150, rotation = TRUE)
}
\arguments{
\item{qsm}{a QSM loaded using `[load_qsm()]`.}

\item{scale}{a factor by which to multiply the `radius_cyl` column to
give line segments the appearance of volume}

\item{rotation}{boolean - indicates whether the plot should display the
tree from 2 angles TRUE, or just one FALSE.}
}
\description{
Simple function to create a diagnostic plot to view QSMs colored
by branching order.
}
\examples{
qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
qsm = load_qsm(qsm_file)
plot_qsm(qsm)
}
