## next steps

#qsm metrics.
# update documentation
# update readme.md

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


## Make a QSM
las = readLAS(system.file("extdata", "tree_0723.laz", package="tReeTraits"))
las = filter_poi(las, Intensity > 44000) # remove foliage returns
las = clean_las(las)
tree_mat = las_to_mat(las)
inputs = lapply(default_qsm_inputs(), function(x) x[2])
run_qsm(tree_mat = tree_mat, unique_id = 'tree_0723', parameter_inputs = inputs, overwrite=TRUE,
    output_results = 'R:/landscape_ecology/projects/canopy-traits/qsm-results/',
    TreeQSM_directory = 'R:/landscape_ecology/projects/canopy-traits/docs/TreeQSM/')

qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
qsm = load_qsm(qsm_file)

qsm_bole = dplyr::filter(qsm, branching_order == 0)

get_center_of_mass(qsm_bole)

