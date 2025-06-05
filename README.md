---
editor_options: 
  markdown: 
    wrap: 72
---

# tReeTraits

An R Package to generate data on tree architecture from terrestrial
lidar scans

# Summary

`tReeTraits` is an `R` package designed to help extract traits on tree
architecture from terrestrial lidar data representing individually
segmented trees, especially traits related to windfirmness (e.g., crown
area, volume, stem taper, branch size distribution, etc). The package
itself draws heavily from other available software especially
<a href='https://github.com/InverseTampere/TreeQSM/'>TreeQSM</a>,
<a href=https://r-lidar.github.io/lidRbook/>lidR</a>, and
<a href=https://github.com/lmterryn/ITSMe>ITSMe</a>. The package brings
these elements together into one package following the methods of Cannon
et al. XXXX.

The package offers functions to (1) pre-process individually segmented
trees, (2) generate a QSM from R by calling a MATLAB program, and (3)
calculate all traits from point clouds and QSMs.

-   **Preprocessing**: Functions to load, recenter, normalize and rotate
    trees, as well as remove vegetation from the vicinity of the bole.
-   **QA/QC**: Function to plot all trees to identify mistakes in
    segmentation
-   **Calculate tree traits from point cloud**: Calculate basic traits
    derived from the point cloud including height, crown width, area,
    and volume.
-   **Make Quantitative Structure Models**: Set parameters and run
    TreeQSM in Matlab without leaving R. This allows calculating
    additional traits on branching architecture
-   **Calculate tree trats from QSMs**: With input from TreeQSM,
    calculate additional traits related to trunk taper, tilt, biomass,
    and branching architecture.

# Package installation

This package has dependencies not on CRAN that must be installed
including `spanner`, `ITSMe`, and `TreeLS`

```{r}
install.packages('lidR')
install.packages('remotes') #allows installation of github packages
remotes::install_github('bi0m3trics/spanner')
remotes::install_github('Imterryn/ITSMe')
remotes::install_github('tiagodc/TreeLS')
```

Get the latest released version of \`tReeTraits from github

```{r}
remotes::install_github('jbcannon/tReeTraits')
```

If you plan to use the QSM aspects of the package, you will also need to
install and license the newest version of Matlab which is requred to
create QMSs.

<a href=https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html>
**Download and Install Matlab**</a>

Next, you must install the the
<a href=https://www.mathworks.com/matlabcentral/answers/4707-how-can-i-download-parallel-computing-toolbox>Parallel
Computing Toolbox</a> and the
<a href=https://www.mathworks.com/products/statistics.html>Statistics
and Machine Learning Toolbox</a> from Matlab.

# Workflow

## Pre-processing

`tReeTraits` contains functions to load, recenter, normalize and rotate
trees, as well as remove vegetation from the vicinity of the bole.

### `clean_las()`

```{r}
library(tReeTraits)
library(lidR)
# load and view example data
las = readLAS(system.file("extdata", "tree_0723.las", package="tReeTraits"))
plot(las)
# clean_las will automatically recenter, normalize and remove vegetation
las = clean_las(las, bole_height=2)
plot(las)

```

![](img/clean_las_ex.jpg){width="300"}

Pine tree with vegetation around bole removed.

### `normalize_las()`

You can normalize and recenter trees with individual functions to
`normalize_las` and `recenter_las`

```{r}
library(lidR)
library(tReeTraits)
las = readLAS(system.file("extdata", "tree_0129.laz", package="tReeTraits"))
las_norm = normalize_las(las)
# view histogram of Z values which now ranging from  0 to 11 m
par(mfrow = 1:2)
hist(las$Z)
hist(las_norm$Z)
```

![](img/normalize_histogram.jpg){width="300"}

Histogram of Z values from lidar scan of pine tree before and after
normalization.

### `recenter_las()`

```{r}
library(lidR)
library(tReeTraits)
las = readLAS(system.file("extdata", "tree_0744.laz", package="tReeTraits"))
# view histogram of original X/Y values
par(mfrow=c(2,2))
hist(las$X)
hist(las$Y)
las = recenter_las(las)
# view histogram of X/Y values centered on 0,0
hist(las$X)
hist(las$Y)
```

![Histograms of X and Y values before (top) and after (bottom)
recentering](img/recenter_las.jpg){width="600"}

### `rotate_las_z`

The function `rotate_las` will rotate a tree about its Z axis which is
useful for analyzing crown traits from multiple angles (although this is
also handled automatically for some functions.

```{r}
library(lidR)
library(tReeTraits)
las = readLAS(system.file("extdata", "tree_0129.laz", package="tReeTraits"))
las_norm = normalize_las(las)
# view histogram of Z values which now ranging from  0 to 11 m
par(mfrow = 1:2)
hist(las$Z)
hist(las_norm$Z)
```

### `plot_tree()`

You may want to plot the tree you have loaded to validate proper
segmentation or look for other errors. `plot_tree()` generates a 3 panel
figure showing two profiles and a birds-eye-view of the tree.

```{r}
library(lidR)
library(tReeTraits)
las = readLAS(system.file("extdata", "tree_0129.laz", package="tReeTraits"))
plot_tree(las)
```

![Figure illustrating 3 views of
tree_0129](img/plot_tree.jpg){width="400"}

## Traits from point cloud

### Basic point cloud measurements

Once the point cloud is loaded you can calculate many tree traits
directly from the point cloud including `get_height()`, `get_width()`,
`get_dbh()`, `get_crown_base()`,

```{r}
library(lidR)
library(tReeTraits)
las = readLAS(system.file("extdata", "tree_0129.laz", package="tReeTraits"))
las = clean_las(las)
height = get_height(las)
crown_width = get_width(las)[1]
dbh = get_dbh(las, select_n=30)
cbh = get_crown_base(las, threshold=0.25, sustain=2)
#undocumented function, see `full_diagnostic_plot()`
basics_diagnostic_plot(las, height, cbh, crown_width, dbh)
```

![Figure illustrating basic canopy measurements of
tree_0129](img/basic_tree_meas.jpg){width="150"}

### Crown area and lever arm

Several canopy measurement metrics related to profile area, volume, and
lacunarity are availalble through the functions `convex_hull_2D()` ,
`voxel_hull_2D()`, `get_lacunarity()` , `get_crown_lever_arm()`.

`convex_hull_2D()`, `voxel_hull_2D()`, `get_lacunarity()` generate
verticle profiles of the crown and measure its area as a convex hull, or
voxelize the crown at a specified resolution to estimate area.
Lacunarity is represents the proportion of "holes" in the canopy and is
defined as $1 - \frac{A_{voxel}}{A_{convex}}$, where A presents the
crown area estimated as a voxel or convex hull, respectively.

```{r}
library(lidR)
library(tReeTraits)
library(sf)
las = readLAS(system.file("extdata", "tree_0129.laz", package="tReeTraits"))
las = clean_las(las)
cbh = get_crown_base(las, threshold=0.25, sustain=2)
las = segment_crown(las, crown_base_height = cbh)
st_area(convex_hull_2D(las)) # ~42 m2
st_area(voxel_hull_2D(las)) # ~33 m2
get_lacunarity(las) # ~ 20%
hull_diagnostic_plot(las) #undocumented function, see `full_diagnostic_plot()`

```

![Illustration of convex hull (dashed line) and voxel hull (green
feature) from tree_0129. The proportion of whitespace within the convex
hull represents lacunarity (\~20%)](img/get_lacunarity.jpg){width="300"}

The function `get_crown_lever_arm()` is relevant to windfirmness. The
function measures the crown area in different segments, and multiplies
each times the height at which they occur. This integration of crown
area and the height at which it occurs is related to total wind forces
acting on a tree. The value is proportional to wind drag on a assuming a
constant wind velocity and drag.

```{r}
library(lidR)
library(tReeTraits)
las = readLAS(system.file("extdata", "tree_0129.laz", package="tReeTraits"))
las = clean_las(las)
cbh = get_crown_base(las, threshold=0.25, sustain=2)
las = segment_crown(las, crown_base_height = cbh)
get_crown_lever_arm(las)
```

### Crown volume

Two methods are available `tReeTraits` for measuring crown volume.
Volume can be measured as an alpha volume using the `ITSMe` package.
Alternatively, crown pixels can be voxelized at a defined resolution and
volume can be estimated as the $$ n_{voxels} \times V_{voxels} $$ where
$n$ is the number of voxels, and $V$ is the volume of the defined voxels
based on their resolution.

```{r}
library(lidR)
library(tReeTraits)
las = readLAS(system.file("extdata", "tree_0129.laz", package="tReeTraits"))
las = clean_las(las)
cbh = get_crown_base(las, threshold=0.25, sustain=2)
las = segment_crown(las, crown_base_height = cbh)
get_crown_volume_alpha(las) # may take a few seconds (~81 m3)
get_crown_volume_voxel(las) # (~31 m3) # this volume will be smaller and approach zero as resolution gets smaller.

# You can visualize the alpha volume with the `ITSMe` package.
library(ITMSe)
crown = filter_poi(las, Crown==1)
vox = voxelize_points(crown, res=0.1)
alpha_volume_pc(vox@data[,c('X','Y','Z')], alpha=0.5, plot=TRUE)
# plotting of the voxelized crown is not currently implemented.

```

![Figure illustrating crown volume of tree_0129 estimated using an alpha
volume from the `ITSMe`
package.](img/crown_alpha_volume.jpg){width="200"}

## Traits From QSM

## Generating QSMs via Matlab:Treeqsm integration

### Set QSM/Matlab paths

### Generate parameters

### Run QSM

### Make QSM JPG
