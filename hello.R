
library(lidR)
library(tReeTraits)

#show how data gets cleaned
filename = system.file("extdata", "tree_0723.las", package="tReeTraits")
las = lidR::readLAS(filename)
plot(las)
las = clean_las(las, bole_height=3)
plot(las)

# auto generate to check segementation
plot_tree(las)

# Get basic tree attributes
get_height(las)
get_width(las)
get_dbh(las, select_n=30)

# Crown traits from point cloud
cbh = get_crown_base(las, threshold=0.25, sustain=2)
print(cbh)
las = segment_crown(las, cbh)
plot(las, color='Crown')

#2D Crown Structure
C_hull = convex_hull_2D(las, angle=0)
V_hull = voxel_hull_2D(las, angle=0, resolution = 0.1)

library(ggplot2)
ggplot() + theme_bw() +
  geom_sf(data=C_hull, fill='yellow', alpha=0.5) +
  geom_sf(data=V_hull, fill='blue', alpha=0.5)

# crown area and lacunarity
library(sf)
st_area(convex_hull_2D(las))
st_area(voxel_hull_2D(las))
get_lacunarity(las)
get_crown_lever_arm(las)

# get area statistics from multiple angles/crown rotations
angles = c(0, 30, 60, 90)
output = lapply(angles, function(a) {
  data.frame(
    convex_area = st_area(convex_hull_2D(las, angle=a)),
    voxel_area = st_area(voxel_hull_2D(las, angle=a)),
    lacunarity = get_lacunarity(las,angle=a))
})
output = do.call(rbind, output) # combine rows
output
apply(output, 2, mean) # get means

# 3D crown structure
#Crown volume from alpha shape
V_a = get_crown_volume_alpha(las, resolution = 0.1)
print(V_a)
# Crown volume from voxelization
V_v = get_crown_volume_voxel(las, resolution = 0.1)
print(V_v)

## Make a QSM
las = readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
las = filter_poi(las, Intensity > 44000) # remove foliage returns
las = clean_las(las)
tree_mat = las_to_mat(las)
inputs = lapply(default_qsm_inputs(), function(x) x[2])
run_qsm(tree_mat = tree_mat, unique_id = 'tree_0723', parameter_inputs = inputs, overwrite=TRUE,
    output_results = 'R:/landscape_ecology/projects/canopy-traits/qsm-results/',
    TreeQSM_directory = 'R:/landscape_ecology/projects/canopy-traits/docs/TreeQSM/')

# Get traits from a QSM
qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
qsm = load_qsm(qsm_file)
plot_qsm(qsm)

# Fit taper equation Kozak
fit_taper_Kozak(qsm, dbh = 15.2)$results

# Get information on primary branches
primary_branches = get_primary_branches(qsm)
nrow(primary_branches) #Number of primary branches
median(primary_branches$diam_cm)# median primary branch size

# Get median internode distance
median(internode_distances(qsm))

# Branch distribution, volume weighted statistics, and skewness
branch_dist = branch_size_distribution(qsm, plot=TRUE)
print(branch_dist)

#volume-weighted mean
branch_volume_weighted_stats(qsm, FUN = function(x) mean(x))

#volume-weighted median
branch_volume_weighted_stats(qsm, FUN = function(x) median(x))

# volume-weighted skewness and kurtosis
library(moments)
branch_volume_weighted_stats(qsm, FUN = function(x) skewness(x))
branch_volume_weighted_stats(qsm, FUN = function(x) kurtosis(x))

# Get volume distribution
volume = qsm_volume_distribution(qsm)
View(volume)

# Center of mass and its horizontal offset
get_center_of_mass(qsm)
get_com_offset(qsm)

# Stem deflection (tilt and sweep)
get_stem_tilt(qsm)
sweep = get_stem_sweep(qsm)
    max(sweep$sweep) #maximum
    mean(sweep$sweep) #mean deviation
    sqrt(mean((sweep$sweep)^2)) #RMSE sweep

#need to get QSM FIT statistics



filename = system.file("extdata", "tree_0723.las", package="tReeTraits")
las = lidR::readLAS(filename)
las = clean_las(las, bole_height=3)

# Get basic tree attributes
height = get_height(las)
crown_width = get_width(las)[1]
dbh = get_dbh(las, select_n=30)

# Crown traits from point cloud
cbh = get_crown_base(las, threshold=0.25, sustain=2)
las = segment_crown(las, cbh)

qsm_file = system.file("extdata", "tree_0723_qsm.mat", package='tReeTraits')
qsm = load_qsm(qsm_file)

#diagnostic
diagnostic = full_diagnostic_plot(las, qsm, height, cbh, crown_width, dbh)
plot(diagnostic)
ggsave('diagnostic.jpg', diagnostic, height=8, width=10)
