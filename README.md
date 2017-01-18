trans-diagnostic clustering of brain activity during the imitate / observe task
-------------------------------------------------------------------------------

**cluster_gridsearch.m**

Uses a clustering-stability method to termine the model order (number of clusters) to be used when partioning the data (spoiler: 2 or 3 is the correct answer).

This also could be used to compare distance functions, but that has little influence on the output. This does not do a good job at comparing linkage functions, as our data has some extreme outliers. Therefore, complete, single, average, centroid clustering all produce clusters with single individuals, and a second cluster with the remainder of the data. This produces extremely stable, but uninstering, solutions.

Also, when subdividing the 'oddball' group, the clusters become quite small, as the individuals are quite different from one another.

Outputs `cluster_grid_search.mat`, and `cluster_stability.fig`.

**export.m**

Takes the settings found in `cluster_gridsearch.m`, and runs a GLM on each cluster to test whether each voxel is has a mean value > 0. The outputs of these GLM runs are exported a `.nii` files.

Also plots a dendrogram + distance matrix of the selected cluster settings.

Outputs `cluster_final.mat` and `stat_hist.fig`.

**project-to-surface.sh**

Relies on the connectome-workbench. Not actually functional, but a reminder of how to convert the output `.nii` files from `export.m` to `dscalar.nii` files that can be opened in connectome workbench.

