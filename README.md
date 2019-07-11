# GFDL CM4 Model
[Geophysical Fluid Dynamics Laboratory
(GFDL)](https://www.gfdl.noaa.gov).

The layout of this package includes the following directories:

* src - The source code for the AM4 model
* exec - The build directory with Makefiles for building the model executable
* run - Sample run script
## * analysis - Sample analysis scripts 

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
# From within the AM4 parent directory
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
* [icebergs](https://github.com/NOAA-GFDL/icebergs)
* [ice_param](https://github.com/NOAA-GFDL/ice_param)
* [ocean_BGC](https://github.com/NOAA-GFDL/ocean_BGC)
* [coupler](https://github.com/NOAA-GFDL/FMScoupler)
* [FMS](https://github.com/NOAA-GFDL/FMS) (as [shared](src/shared))
* [mocsy](https://github.com/NOAA-GFDL/mocsy)

## Building CM4

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
make
```
If you would like to change some of the compiler options, there are several different
options to add to the make command.  For example
```
make ISA=-xhost BLD_TYPE=REPRO
```
will replace -msse with -xhost and -O3 with -O2.  The three options for 
`BLD_TYPE` are  
`PROD` (-O3)  
`REPRO` (-O2)    
`DEBUG` (-O0 and other traps)  
All of the make line options can be
found in the [intel.mk](exec/templates/intel.mk) file.

## Running CM4

Included in the run directory is a sample run script for reference.
##To run the CM4 sample experiment, first download the data file
##mentioned in [Obtaining the Input data](#obtaining-the-input-data)
##section.  Modify the variables in the configuration section in the
##sample run script, and then run the script.
##
##The sample data and run script are configured to run on 
##processors.  To run on a different number of processors, or modify the
##experiment, refer to the `README.CM4_run` file included in the CM4
##data tarball.

##Note: The `input.nml` file (found in the CM4 data tarball) contains
##Fortran namelists and namelist variables that modify, at run time, the
##model.  To learn more about the settings in the `input.nml` file,
##please refer to source code where the namelist/variable are defined.


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
