#!/bin/bash

wb_command -cifti-convert -from-nifti \
    /projects/colin/SPINS_hcp/all_3solution_euclidean_ward3/EM_IMI_OBS/spmT_0001.nii \
    /projects/jviviano/data/cluster-glm/assets/cifti_template.nii \
    /projects/jviviano/data/cluster-glm/all_3solution_euclidean_ward_3.dscalar.nii
