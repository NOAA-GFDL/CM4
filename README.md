# GFDL CM4 Model
[Geophysical Fluid Dynamics Laboratory
(GFDL)](https://www.gfdl.noaa.gov).

The layout of this package includes the following directories:

* src - The source code for the CM4 model
* exec - The build directory with Makefiles for building the model executable
* run - Sample run script

## Cloning Instructions

This repository uses [git
submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to
point to other repositories.  Thus, care should be taken when cloning,
and updating the source to ensure all source.  To obtain all source,
use the following git command

```
git clone --recursive https://github.com/NOAA-GFDL/CM4.git
```

The `--recursive` option to `git clone` instructs git to recursively
clone all submodules.  In the event the repository was not cloned
using the `--recursive` option, the following step must be taken to
obtain all sources:

```
# From within the CM4 parent directory
git submodule update --init --recursive
```

## Source Code

All model source is contained in the [src](src) directory.  GFDL
tracks code using the git version control system.  This package
includes a single version of the following GFDL model components.  The
git hash listed corresponds to the commit hash in the internal GFDL
git repository.

Component | Commit Hash
--------- | -----------
atmos_cubed_sphere | b8b05bf650c0d3293b538bdaceb894ba0fd6910b
atmos_drivers | 3be6ed406de2db29766746a69115fd6a47048692
atmos_param | 4fe4ca54a0224ef5c4cf9ebf1010d5b869930a3f
atmos_shared | 2e9d8b770cdb2d70d8d9264e4b2de24213ae21bd
land_lad2 | 154bd2b4bf523f3e699de5017679b156242ec13f 



The following components are available in the
[NOAA-GFDL](https://github.com/NOAA-GFDL) github organization:

* [MOM6](https://github.com/NOAA-GFDL/MOM6)
* [SIS2](https://github.com/NOAA-GFDL/SIS2)
* [GFDL_atmos_cubed_sphere](https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/tree/AM4.0) (as [atmos_cubed_sphere](src/atmos_cubed_sphere))
* [icebergs](https://github.com/NOAA-GFDL/icebergs)
* [ice_param](https://github.com/NOAA-GFDL/ice_param)
* [ocean_BGC](https://github.com/NOAA-GFDL/ocean_BGC)
* [coupler](https://github.com/NOAA-GFDL/FMScoupler)
* [FMS](https://github.com/NOAA-GFDL/FMS) (as [shared](src/shared))
* [mocsy](https://github.com/NOAA-GFDL/mocsy)

## Building CM4

### Containers

The [container folder](container) provides example Dockerfiles and Signularity
definition files to use to build AM4 containers using either GCC/GFORTAN or
Intel oneAPI. There is a script that can be used to build the intel
singularity containers, and the first step of this script can be used with the
other GFDL climate models.

### From source

The [exec](exec) directory contains Makefiles that can be used to
build the CM4 executable.  These Makefiles were generated using the
[Make Makefile (mkmf)](https://github.com/NOAA-GFDL/mkmf) program.
Included in the exec direcgtory is a sample make template file for the
Intel compilers ([intel.mk](exec/templates/intel.mk)).  This make
template can be used on any system with a relatively recent version of
the Intel compilers, the netCDF 4 library and the MPICH2 MPI library.
Included in the [intel.mk](exec/templates/intel.mk) file are
additional settings that can be modified during the build.  


To run the default build (-O3 -msse2), go to the exec directory and
enter the command
```
make HDF_INCLUDE=-I/path/to/hdf5/include
```
Where */path/to/hdf5/include* is the path to your HDF5 include folder where hdf5.mod
is. 

If you would like to change some of the compiler options, there are several different
options to add to the make command.  For example
```
make ISA=-xhost REPRO=on
```
will replace -msse with -xhost and -O3 with -O2.  The three options for 
building are  
`PROD=on` (-O3) Default
`REPRO=on` (-O2)    
`DEBUG=on` (-O0 and other traps)  
All of the make line options can be
found in the [intel.mk](exec/templates/intel.mk) file.

To build with GNU compilers, add `gcc=on` to the `make` line. The make line
options can be found in the [gnu.mk](exec/templates/gnu.mk) file.
## Obtaining the input data

The input data required for running the CM4 model can be found on
[GFDL's data
portal](http://data1.gfdl.noaa.gov/nomads/forms/cm4/) .

The file `CM4_runDir.tar.gz` contains a configured run directory to run a
sample experiment of the CM4 model.  Included in the tar file is a
README.CM4 with more instructions on how to configure the CM4 run
directory.

## Running CM4

Included in the run directory is a sample run script for reference.
To run the CM4 sample experiment, first download the data file
mentioned in [Obtaining the Input data](#obtaining-the-input-data)
section.  Modify the variables in the configuration section in the
sample run script, and then run the script.

The sample data and run script are configured to run on a total of 8127
processors (864 cores 4 threads for the atmosphere and 4671 ocean cores).  
To run on a different number of processors, or modify the
experiment, refer to the `README.CM4` file included in the CM4
data tarball.

Note: The `input.nml` file (found in the CM4 data tarball) contains
Fortran namelists and namelist variables that modify, at run time, the
model.  To learn more about the settings in the `input.nml` file,
please refer to source code where the namelist/variable are defined.


## Model output and Other References

Please refer to the [CM4 data and code
site](http://data1.gfdl.noaa.gov/nomads/forms/cm4/) for details
about where to find model and OBS data used in the papers.

For all analysis figures and pertaining data, please use the CM4
documentation papers as the original reference.

Please direct your questions and feedback to
gfdl.climate.model.info@noaa.gov

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an 'as is' basis and the user assumes responsibility for
its use.  DOC has relinquished control of the information and no
longer has responsibility to protect the integrity, confidentiality,
or availability of the information.  Any claims against the Department
of Commerce stemming from the use of its GitHub project will be
governed by all applicable Federal law.  Any reference to specific
commercial products, processes, or services by service mark,
trademark, manufacturer, or otherwise, does not constitute or imply
their endorsement, recommendation or favoring by the Department of
Commerce.  The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

This project code is made available through GitHub but is managed by
NOAA-GFDL at https://gitlab.gfdl.noaa.gov.
