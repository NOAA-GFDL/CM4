# Makefile for Experiment 'cm4p12_warsaw'

BUILDROOT = $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
SRCROOT = $(abspath $(BUILDROOT)../src)/ # Final '/' added to denote directory

MK_TEMPLATE = $(BUILDROOT)templates/intel.mk

# Override default compile options. See templates/intel.mk for more compile options.
# BLD_TYPE configures the type of build, PROD (production), REPRO (regression) or
# DEBUG (debug options)
BLD_TYPE = PROD
OPENMP = t # Build with OpenMP enabled

include $(MK_TEMPLATE)

SUBMAKEFLAGS = BUILDROOT=$(BUILDROOT) SRCROOT=$(SRCROOT) MK_TEMPLATE=$(MK_TEMPLATE) BLD_TYPE=$(BLD_TYPE)

fms_CM4p0.x: coupler/libcoupler.a sis2/libsis2.a atmos_dyn/libatmos_dyn.a land_lad2/libland_lad2.a atmos_phys/libatmos_phys.a mom6/libmom6.a icebergs/libicebergs.a fms/libfms.a
	$(LD) $^ $(LDFLAGS) -o $@ $(STATIC_LIBS)

fms/libfms.a:  FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)fms $(@F)

atmos_phys/libatmos_phys.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)atmos_phys $(@F)

atmos_dyn/libatmos_dyn.a: atmos_phys/libatmos_phys.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)atmos_dyn $(@F)

sis2/libsis2.a: icebergs/libicebergs.a mom6/libmom6.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)sis2 $(@F)

icebergs/libicebergs.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)icebergs $(@F)

land_lad2/libland_lad2.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)land_lad2 $(@F)

# MOM6 is currently unable to run with OpenMP enabled
mom6/libmom6.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)mom6 $(@F)

coupler/libcoupler.a: atmos_dyn/libatmos_dyn.a sis2/libsis2.a atmos_phys/libatmos_phys.a mom6/libmom6.a land_lad2/libland_lad2.a icebergs/libicebergs.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) OPENMP=$(OPENMP) --directory=$(BUILDROOT)coupler $(@F)

FORCE:

clean:
	echo fms
	$(MAKE) --directory=$(BUILDROOT)fms clean
	echo atmos_phys
	$(MAKE) --directory=$(BUILDROOT)atmos_phys clean
	echo atmos_dyn
	$(MAKE) --directory=$(BUILDROOT)atmos_dyn clean
	echo sis2
	$(MAKE) --directory=$(BUILDROOT)sis2 clean
	echo land_lad2
	$(MAKE) --directory=$(BUILDROOT)land_lad2 clean
	echo icebergs
	$(MAKE) --directory=$(BUILDROOT)icebergs clean
	echo mom6
	$(MAKE) --directory=$(BUILDROOT)mom6 clean
	echo coupler
	$(MAKE) --directory=$(BUILDROOT)coupler clean

clean_all: clean
	$(RM) *.x

