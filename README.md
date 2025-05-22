---
title: "tReeTraits"
output:
  html_document:
    toc: true
---
An R Package to generate data on tree architecture from terrestrial lidar scans

# Summary

`tReeTraits` is an `R` package designed to help extract traits on tree architecture
from terrestrial lidar data representing individually segmented trees, especially
traits related to windfirmness (e.g., crown area, volume, stem taper, branch size
distribution, etc). The package itself draws heavily from other available software
especially <a href='https://github.com/InverseTampere/TreeQSM/'>TreeQSM</a>, 
<a href=https://r-lidar.github.io/lidRbook/>lidR</a>,
and <a href=https://github.com/lmterryn/ITSMe>ITSMe</a>. The package brings
these elements together into one package following the methods of Cannon et
al. XXXX.

The package offers functions to (1) pre-process individually segmented trees, (2)
generate a QSM from R by calling a MATLAB program, and (3) calculate all traits
from point clouds and QSMs.

* **Preprocessing**: Functions to load, recenter, normalize and rotate trees, as 
well as remove vegetation from the vicinity of the bole.
* **QA/QC**: Function to plot all trees to identify mistakes in segmentation
* **Calculate tree traits from point cloud**: Calculate basic traits derived from
the point cloud including height, crown width, area, and volume.
* **Make Quantitative Structure Models**: Set parameters and run TreeQSM in Matlab
without leaving R. This allows calculating additional traits on branching
architecture
* **Calculate tree trats from QSMs**: With input from TreeQSM, calculate additional
traits related to trunk taper, tilt, biomass, and branching architecture.

# Workflow

## Package installation

## Pre-processing

## CLeanup las

-show side by size of tree-0623

### Normalize and Recenter

-show side by size of tree-0623

### Remove vegetation

-show side by size of tree-0623

## Traits from point cloud

##  Traits From QSM

## Generating QSMs via Matlab:Treeqsm integration

### Set QSM/Matlab paths

### Generate parameters

### Run QSM

### Make QSM JPG





