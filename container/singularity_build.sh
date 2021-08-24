#!/bin/sh

singularity build -f intel2021.2_netcdfc4.7.4_ubuntu.sif Singularity.intel_netcdf
singularity build -f cm4_2021.02_ubuntu_intel.sif Singularity.intel_cm4
