# Makefile created by mkmf 19.3.0

CPPDEFS = -DINTERNAL_FILE_NML 

OTHERFLAGS = -I$(BUILDROOT)fms/build/.mods -I$(SRCROOT)FMS/include

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libicebergs.a
icebergs.o: $(SRCROOT)/icebergs/icebergs.F90 icebergs_framework.o icebergs_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/icebergs/icebergs.F90
icebergs_framework.o: $(SRCROOT)/icebergs/icebergs_framework.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/icebergs/icebergs_framework.F90
icebergs_io.o: $(SRCROOT)/icebergs/icebergs_io.F90 icebergs_framework.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/icebergs/icebergs_io.F90
./icebergs.F90: $(SRCROOT)/icebergs/icebergs.F90
	cp $(SRCROOT)/icebergs/icebergs.F90 .
./icebergs_framework.F90: $(SRCROOT)/icebergs/icebergs_framework.F90
	cp $(SRCROOT)/icebergs/icebergs_framework.F90 .
./icebergs_io.F90: $(SRCROOT)/icebergs/icebergs_io.F90
	cp $(SRCROOT)/icebergs/icebergs_io.F90 .
SRC = $(SRCROOT)/icebergs/icebergs_framework.F90 $(SRCROOT)/icebergs/icebergs.F90 $(SRCROOT)/icebergs/icebergs_io.F90
OBJ = icebergs_framework.o icebergs.o icebergs_io.o
OFF = $(SRCROOT)/icebergs/icebergs.F90 $(SRCROOT)/icebergs/icebergs_framework.F90 $(SRCROOT)/icebergs/icebergs_io.F90
clean: neat
	$(RM) -f .libicebergs.a.cppdefs $(OBJ) libicebergs.a *.mod

neat:
	$(RM) -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libicebergs.a: $(OBJ)
	$(AR) $(ARFLAGS) libicebergs.a $(OBJ)

