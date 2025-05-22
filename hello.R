rm(list=ls())
library(lidR)
library(tReeTraits)

las = lidR::readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
plot(las)
las = clean_las(las, bole_height=2)
plot(las)

# Get basic tree attributes
get_height(las)
get_width(las)
get_dbh(las, select_n=20)

# Crown traits from point cloud
cbh = get_crown_base(las, threshold=0.25, sustain=2); print(cbh)
las = segment_crown(las, cbh)
plot(las, color='Crown')

#2D Crown Structure
C_hull = convex_hull_2D(las, angle=0)
V_hull = voxel_hull_2D(las, angle=0, resolution = 0.1)

library(ggplot2)
ggplot() + theme_bw() +
  geom_sf(data=C_hull, fill='yellow', alpha=0.5) +
  geom_sf(data=V_hull, fill='blue', alpha=0.5)

# 3D crown structure
#Crown volume from alpha shape
V_a = get_crown_volume_alpha(las, resolution = 0.1)
print(V_a)
# Crown volume from voxelization
V_v = get_crown_volume_voxel(las, resolution = 0.1)
print(V_v)

#show how data gets cleaned
filename = system.file("extdata", "tree_0723.laz", package="tReeTraits")
las = readLAS(filename)
las_cleaned = tReeTraits::clean_las(las)

plot(las)
plot(las_cleaned)


las = lidR::readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
las = clean_las(las)
# remove bole
las = TreeLS::stemPoints(las)
eigen = spanner::eigen_metrics(las)
#add eigens
for(i in c('Verticality', 'Sphericity', 'Linearity')) {
  las = add_lasattribute(las, eigen[[i]], i, i)
}
crown_score = (as.numeric((las$Sphericity > 0.3)) +  as.numeric(las$Verticality < 0.9) + as.numeric(las$Linearity < 0.55) + as.numeric(las$Intensity < 41000))
crown_score = crown_score * (!las$Stem)
crown_score = crown_score>=3
las = add_lasattribute(las, as.numeric(crown_score), 'Crown', 'Crown')
plot(las, color='Crown', legend=TRUE)

library(ggplot2)
library(tidyverse)
las@data %>% filter(!is.na(`Crown Points`)) %>%
  ggplot(aes(x=Intensity, alpha=0.1, fill=as.factor(`Crown Points`))) + geom_density()
plot(las, color='Crown')


crown = filter_poi(las, Crown == 1)
stems = filter_poi(las, Crown == 0)
