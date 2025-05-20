rm(list=ls())
library(lidR)

las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
las = normalize_las(las)
las = recenter_las(las)
plot(las)


filename = system.file("extdata", "tree_0723.laz", package="tReeTraits")
las = tReeTraits::clean_las(filename)

plot(las)
