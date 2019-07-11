# Makefile created by mkmf 19.3.0

CPPDEFS = -DINTERNAL_FILE_NML -g 

OTHERFLAGS = -I$(BUILDROOT)/sis2 -I$(BUILDROOT)/atmos_dyn -I$(BUILDROOT)/land_lad2 -I$(BUILDROOT)/icebergs -I$(BUILDROOT)/atmos_phys -I$(BUILDROOT)/mom6 -I$(BUILDROOT)/fms

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libcoupler.a
atm_land_ice_flux_exchange.o: $(SRCROOT)FMScoupler/atm_land_ice_flux_exchange.F90 surface_flux.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/atm_land_ice_flux_exchange.F90
coupler_main.o: $(SRCROOT)FMScoupler/coupler_main.F90 flux_exchange.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/coupler_main.F90
flux_exchange.o: $(SRCROOT)FMScoupler/flux_exchange.F90 atm_land_ice_flux_exchange.o land_ice_flux_exchange.o ice_ocean_flux_exchange.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/flux_exchange.F90
ice_ocean_flux_exchange.o: $(SRCROOT)FMScoupler/ice_ocean_flux_exchange.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/ice_ocean_flux_exchange.F90
land_ice_flux_exchange.o: $(SRCROOT)FMScoupler/land_ice_flux_exchange.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/land_ice_flux_exchange.F90
surface_flux.o: $(SRCROOT)FMScoupler/surface_flux.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)FMScoupler/surface_flux.F90
./atm_land_ice_flux_exchange.F90: $(SRCROOT)FMScoupler/atm_land_ice_flux_exchange.F90
	cp $(SRCROOT)FMScoupler/atm_land_ice_flux_exchange.F90 .
./surface_flux.F90: $(SRCROOT)FMScoupler/surface_flux.F90
	cp $(SRCROOT)FMScoupler/surface_flux.F90 .
./flux_exchange.F90: $(SRCROOT)FMScoupler/flux_exchange.F90
	cp $(SRCROOT)FMScoupler/flux_exchange.F90 .
./ice_ocean_flux_exchange.F90: $(SRCROOT)FMScoupler/ice_ocean_flux_exchange.F90
	cp $(SRCROOT)FMScoupler/ice_ocean_flux_exchange.F90 .
./land_ice_flux_exchange.F90: $(SRCROOT)FMScoupler/land_ice_flux_exchange.F90
	cp $(SRCROOT)FMScoupler/land_ice_flux_exchange.F90 .
./coupler_main.F90: $(SRCROOT)FMScoupler/coupler_main.F90
	cp $(SRCROOT)FMScoupler/coupler_main.F90 .
SRC = $(SRCROOT)FMScoupler/coupler_main.F90 $(SRCROOT)FMScoupler/atm_land_ice_flux_exchange.F90 $(SRCROOT)FMScoupler/flux_exchange.F90 $(SRCROOT)FMScoupler/ice_ocean_flux_exchange.F90 $(SRCROOT)FMScoupler/land_ice_flux_exchange.F90 $(SRCROOT)FMScoupler/surface_flux.F90
OBJ = coupler_main.o atm_land_ice_flux_exchange.o flux_exchange.o ice_ocean_flux_exchange.o land_ice_flux_exchange.o surface_flux.o
OFF = $(SRCROOT)FMScoupler/atm_land_ice_flux_exchange.F90 $(SRCROOT)FMScoupler/surface_flux.F90 $(SRCROOT)FMScoupler/flux_exchange.F90 $(SRCROOT)FMScoupler/ice_ocean_flux_exchange.F90 $(SRCROOT)FMScoupler/land_ice_flux_exchange.F90 $(SRCROOT)FMScoupler/coupler_main.F90
clean: neat
	$(RM) -f .libcoupler.a.cppdefs $(OBJ) libcoupler.a *.mod

neat:
	$(RM) -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libcoupler.a: $(OBJ)
	$(AR) $(ARFLAGS) libcoupler.a $(OBJ)

