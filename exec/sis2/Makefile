# Makefile created by mkmf 19.3.0


CPPDEFS = -DINTERNAL_FILE_NML -g 

OTHERFLAGS = -I$(BUILDROOT)icebergs -I$(BUILDROOT)mom6 -I$(BUILDROOT)fms

include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libsis2.a
SIS2_ice_thm.o: $(SRCROOT)/SIS2/src/SIS2_ice_thm.F90 $(SRCROOT)MOM6/src//framework/version_variable.h ice_thm.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)MOM6/src//framework	$(SRCROOT)/SIS2/src/SIS2_ice_thm.F90
SIS_continuity.o: $(SRCROOT)/SIS2/src/SIS_continuity.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_hor_grid.o ice_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_continuity.F90
SIS_ctrl_types.o: $(SRCROOT)/SIS2/src/SIS_ctrl_types.F90 SIS_dyn_trans.o SIS_fast_thermo.o SIS_slow_thermo.o SIS_hor_grid.o ice_grid.o SIS_types.o SIS_optics.o SIS_diag_mediator.o SIS_sum_output_type.o SIS_tracer_flow_control.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/SIS_ctrl_types.F90
SIS_debugging.o: $(SRCROOT)/SIS2/src/SIS_debugging.F90 $(SRCROOT)MOM6/src//framework/version_variable.h SIS_hor_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)MOM6/src//framework	$(SRCROOT)/SIS2/src/SIS_debugging.F90
SIS_diag_mediator.o: $(SRCROOT)/SIS2/src/SIS_diag_mediator.F90 $(SRCROOT)MOM6/src//framework/version_variable.h SIS_hor_grid.o ice_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)MOM6/src//framework	$(SRCROOT)/SIS2/src/SIS_diag_mediator.F90
SIS_dyn_bgrid.o: $(SRCROOT)/SIS2/src/SIS_dyn_bgrid.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_debugging.o SIS_hor_grid.o ice_ridge.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_dyn_bgrid.F90
SIS_dyn_cgrid.o: $(SRCROOT)/SIS2/src/SIS_dyn_cgrid.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_debugging.o SIS_hor_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_dyn_cgrid.F90
SIS_dyn_trans.o: $(SRCROOT)/SIS2/src/SIS_dyn_trans.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_debugging.o SIS_sum_output.o SIS_types.o SIS_utils.o SIS_hor_grid.o ice_grid.o SIS2_ice_thm.o SIS_dyn_bgrid.o SIS_dyn_cgrid.o SIS_tracer_flow_control.o SIS_transport.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_dyn_trans.F90
SIS_fast_thermo.o: $(SRCROOT)/SIS2/src/SIS_fast_thermo.F90 $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_debugging.o SIS_optics.o SIS_types.o ice_boundary_types.o SIS_hor_grid.o ice_grid.o SIS2_ice_thm.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)MOM6/src//framework	$(SRCROOT)/SIS2/src/SIS_fast_thermo.F90
SIS_fixed_initialization.o: $(SRCROOT)/SIS2/src/SIS_fixed_initialization.F90 $(SRCROOT)MOM6/src//framework/version_variable.h SIS_debugging.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)MOM6/src//framework	$(SRCROOT)/SIS2/src/SIS_fixed_initialization.F90
SIS_get_input.o: $(SRCROOT)/SIS2/src/SIS_get_input.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/SIS_get_input.F90
SIS_hor_grid.o: $(SRCROOT)/SIS2/src/SIS_hor_grid.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_hor_grid.F90
SIS_optics.o: $(SRCROOT)/SIS2/src/SIS_optics.F90 $(SRCROOT)MOM6/src//framework/version_variable.h ice_shortwave_dEdd.o SIS2_ice_thm.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)MOM6/src//framework	$(SRCROOT)/SIS2/src/SIS_optics.F90
SIS_slow_thermo.o: $(SRCROOT)/SIS2/src/SIS_slow_thermo.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_sum_output.o SIS_types.o SIS_utils.o SIS_hor_grid.o ice_grid.o ice_spec.o SIS_tracer_flow_control.o SIS2_ice_thm.o SIS_optics.o SIS_transport.o SIS_tracer_registry.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_slow_thermo.F90
SIS_sum_output.o: $(SRCROOT)/SIS2/src/SIS_sum_output.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_types.o SIS_hor_grid.o ice_grid.o SIS2_ice_thm.o SIS_sum_output_type.o SIS_tracer_flow_control.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_sum_output.F90
SIS_sum_output_type.o: $(SRCROOT)/SIS2/src/SIS_sum_output_type.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/SIS_sum_output_type.F90
SIS_tracer_advect.o: $(SRCROOT)/SIS2/src/SIS_tracer_advect.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_hor_grid.o ice_grid.o SIS_tracer_registry.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_tracer_advect.F90
SIS_tracer_flow_control.o: $(SRCROOT)/SIS2/src/SIS_tracer_flow_control.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o ice_grid.o SIS_tracer_registry.o SIS_hor_grid.o ice_age_tracer.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_tracer_flow_control.F90
SIS_tracer_registry.o: $(SRCROOT)/SIS2/src/SIS_tracer_registry.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_debugging.o SIS_diag_mediator.o SIS_hor_grid.o ice_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_tracer_registry.F90
SIS_transcribe_grid.o: $(SRCROOT)/SIS2/src/SIS_transcribe_grid.F90 SIS_hor_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/SIS_transcribe_grid.F90
SIS_transport.o: $(SRCROOT)/SIS2/src/SIS_transport.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_tracer_registry.o SIS_tracer_advect.o SIS_continuity.o SIS_hor_grid.o ice_grid.o ice_ridge.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_transport.F90
SIS_types.o: $(SRCROOT)/SIS2/src/SIS_types.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h SIS_hor_grid.o ice_grid.o SIS2_ice_thm.o SIS_diag_mediator.o SIS_debugging.o SIS_sum_output_type.o SIS_tracer_registry.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/SIS_types.F90
SIS_utils.o: $(SRCROOT)/SIS2/src/SIS_utils.F90 SIS_diag_mediator.o SIS_debugging.o SIS_hor_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/SIS_utils.F90
combined_ice_ocean_driver.o: $(SRCROOT)/SIS2/src/combined_ice_ocean_driver.F90 ice_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/combined_ice_ocean_driver.F90
ice_age_tracer.o: $(SRCROOT)/SIS2/src/ice_age_tracer.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_diag_mediator.o SIS_tracer_registry.o SIS_hor_grid.o SIS_utils.o ice_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/ice_age_tracer.F90
ice_albedo.o: $(SRCROOT)ice_param/ice_albedo.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ice_albedo.F90
ice_boundary_types.o: $(SRCROOT)/SIS2/src/ice_boundary_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/ice_boundary_types.F90
ice_grid.o: $(SRCROOT)/SIS2/src/ice_grid.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/ice_grid.F90
ice_model.o: $(SRCROOT)/SIS2/src/ice_model.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h SIS_debugging.o SIS_diag_mediator.o SIS_get_input.o SIS_sum_output.o SIS_transcribe_grid.o ocean_albedo.o ocean_rough.o ice_type.o ice_boundary_types.o SIS_ctrl_types.o SIS_types.o SIS_utils.o SIS_hor_grid.o SIS_fixed_initialization.o ice_grid.o ice_spec.o SIS_tracer_registry.o SIS_tracer_flow_control.o SIS_dyn_trans.o SIS_slow_thermo.o SIS_fast_thermo.o SIS_optics.o SIS2_ice_thm.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/ice_model.F90
ice_ridge.o: $(SRCROOT)/SIS2/src/ice_ridge.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_memory_macros.h SIS_diag_mediator.o SIS_hor_grid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/SIS2/config_src/dynamic_symmetric -I$(SRCROOT)MOM6/src//framework -I$(SRCROOT)/SIS2/src	$(SRCROOT)/SIS2/src/ice_ridge.F90
ice_shortwave_dEdd.o: $(SRCROOT)/SIS2/src/ice_shortwave_dEdd.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/ice_shortwave_dEdd.F90
ice_spec.o: $(SRCROOT)/SIS2/src/ice_spec.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/ice_spec.F90
ice_thm.o: $(SRCROOT)/SIS2/src/ice_thm.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/ice_thm.F90
ice_type.o: $(SRCROOT)/SIS2/src/ice_type.F90 SIS_hor_grid.o ice_grid.o SIS2_ice_thm.o SIS_debugging.o SIS_diag_mediator.o SIS_types.o SIS_ctrl_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/SIS2/src/ice_type.F90
ocean_albedo.o: $(SRCROOT)ice_param/ocean_albedo.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ocean_albedo.F90
ocean_rough.o: $(SRCROOT)ice_param/ocean_rough.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ice_param/ocean_rough.F90
./ice_model.F90: $(SRCROOT)/SIS2/src/ice_model.F90
	cp $(SRCROOT)/SIS2/src/ice_model.F90 .
./SIS_slow_thermo.F90: $(SRCROOT)/SIS2/src/SIS_slow_thermo.F90
	cp $(SRCROOT)/SIS2/src/SIS_slow_thermo.F90 .
./ice_thm.F90: $(SRCROOT)/SIS2/src/ice_thm.F90
	cp $(SRCROOT)/SIS2/src/ice_thm.F90 .
./SIS_ctrl_types.F90: $(SRCROOT)/SIS2/src/SIS_ctrl_types.F90
	cp $(SRCROOT)/SIS2/src/SIS_ctrl_types.F90 .
./SIS2_memory.h: $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h
	cp $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h .
./SIS2_ice_thm.F90: $(SRCROOT)/SIS2/src/SIS2_ice_thm.F90
	cp $(SRCROOT)/SIS2/src/SIS2_ice_thm.F90 .
./SIS_dyn_bgrid.F90: $(SRCROOT)/SIS2/src/SIS_dyn_bgrid.F90
	cp $(SRCROOT)/SIS2/src/SIS_dyn_bgrid.F90 .
./ice_type.F90: $(SRCROOT)/SIS2/src/ice_type.F90
	cp $(SRCROOT)/SIS2/src/ice_type.F90 .
./ice_spec.F90: $(SRCROOT)/SIS2/src/ice_spec.F90
	cp $(SRCROOT)/SIS2/src/ice_spec.F90 .
./SIS_diag_mediator.F90: $(SRCROOT)/SIS2/src/SIS_diag_mediator.F90
	cp $(SRCROOT)/SIS2/src/SIS_diag_mediator.F90 .
./SIS_dyn_cgrid.F90: $(SRCROOT)/SIS2/src/SIS_dyn_cgrid.F90
	cp $(SRCROOT)/SIS2/src/SIS_dyn_cgrid.F90 .
./ocean_rough.F90: $(SRCROOT)ice_param/ocean_rough.F90
	cp $(SRCROOT)ice_param/ocean_rough.F90 .
./combined_ice_ocean_driver.F90: $(SRCROOT)/SIS2/src/combined_ice_ocean_driver.F90
	cp $(SRCROOT)/SIS2/src/combined_ice_ocean_driver.F90 .
./ice_grid.F90: $(SRCROOT)/SIS2/src/ice_grid.F90
	cp $(SRCROOT)/SIS2/src/ice_grid.F90 .
./SIS_get_input.F90: $(SRCROOT)/SIS2/src/SIS_get_input.F90
	cp $(SRCROOT)/SIS2/src/SIS_get_input.F90 .
./ice_albedo.F90: $(SRCROOT)ice_param/ice_albedo.F90
	cp $(SRCROOT)ice_param/ice_albedo.F90 .
./SIS_utils.F90: $(SRCROOT)/SIS2/src/SIS_utils.F90
	cp $(SRCROOT)/SIS2/src/SIS_utils.F90 .
./SIS_fixed_initialization.F90: $(SRCROOT)/SIS2/src/SIS_fixed_initialization.F90
	cp $(SRCROOT)/SIS2/src/SIS_fixed_initialization.F90 .
./SIS_optics.F90: $(SRCROOT)/SIS2/src/SIS_optics.F90
	cp $(SRCROOT)/SIS2/src/SIS_optics.F90 .
./SIS_hor_grid.F90: $(SRCROOT)/SIS2/src/SIS_hor_grid.F90
	cp $(SRCROOT)/SIS2/src/SIS_hor_grid.F90 .
./SIS_sum_output.F90: $(SRCROOT)/SIS2/src/SIS_sum_output.F90
	cp $(SRCROOT)/SIS2/src/SIS_sum_output.F90 .
./ice_age_tracer.F90: $(SRCROOT)/SIS2/src/ice_age_tracer.F90
	cp $(SRCROOT)/SIS2/src/ice_age_tracer.F90 .
./SIS_debugging.F90: $(SRCROOT)/SIS2/src/SIS_debugging.F90
	cp $(SRCROOT)/SIS2/src/SIS_debugging.F90 .
./SIS_tracer_registry.F90: $(SRCROOT)/SIS2/src/SIS_tracer_registry.F90
	cp $(SRCROOT)/SIS2/src/SIS_tracer_registry.F90 .
./SIS_fast_thermo.F90: $(SRCROOT)/SIS2/src/SIS_fast_thermo.F90
	cp $(SRCROOT)/SIS2/src/SIS_fast_thermo.F90 .
./SIS_tracer_advect.F90: $(SRCROOT)/SIS2/src/SIS_tracer_advect.F90
	cp $(SRCROOT)/SIS2/src/SIS_tracer_advect.F90 .
./SIS_transcribe_grid.F90: $(SRCROOT)/SIS2/src/SIS_transcribe_grid.F90
	cp $(SRCROOT)/SIS2/src/SIS_transcribe_grid.F90 .
./ice_boundary_types.F90: $(SRCROOT)/SIS2/src/ice_boundary_types.F90
	cp $(SRCROOT)/SIS2/src/ice_boundary_types.F90 .
./SIS_dyn_trans.F90: $(SRCROOT)/SIS2/src/SIS_dyn_trans.F90
	cp $(SRCROOT)/SIS2/src/SIS_dyn_trans.F90 .
./SIS_tracer_flow_control.F90: $(SRCROOT)/SIS2/src/SIS_tracer_flow_control.F90
	cp $(SRCROOT)/SIS2/src/SIS_tracer_flow_control.F90 .
./SIS_continuity.F90: $(SRCROOT)/SIS2/src/SIS_continuity.F90
	cp $(SRCROOT)/SIS2/src/SIS_continuity.F90 .
./SIS_memory_macros.h: $(SRCROOT)/SIS2/src/SIS_memory_macros.h
	cp $(SRCROOT)/SIS2/src/SIS_memory_macros.h .
./MOM_memory_macros.h: $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h
	cp $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h .
./SIS_types.F90: $(SRCROOT)/SIS2/src/SIS_types.F90
	cp $(SRCROOT)/SIS2/src/SIS_types.F90 .
./ice_ridge.F90: $(SRCROOT)/SIS2/src/ice_ridge.F90
	cp $(SRCROOT)/SIS2/src/ice_ridge.F90 .
./version_variable.h: $(SRCROOT)MOM6/src//framework/version_variable.h
	cp $(SRCROOT)MOM6/src//framework/version_variable.h .
./ice_shortwave_dEdd.F90: $(SRCROOT)/SIS2/src/ice_shortwave_dEdd.F90
	cp $(SRCROOT)/SIS2/src/ice_shortwave_dEdd.F90 .
./SIS_transport.F90: $(SRCROOT)/SIS2/src/SIS_transport.F90
	cp $(SRCROOT)/SIS2/src/SIS_transport.F90 .
./ocean_albedo.F90: $(SRCROOT)ice_param/ocean_albedo.F90
	cp $(SRCROOT)ice_param/ocean_albedo.F90 .
./SIS_sum_output_type.F90: $(SRCROOT)/SIS2/src/SIS_sum_output_type.F90
	cp $(SRCROOT)/SIS2/src/SIS_sum_output_type.F90 .
SRC = $(SRCROOT)/SIS2/src/SIS_ctrl_types.F90 $(SRCROOT)/SIS2/src/ice_model.F90 $(SRCROOT)/SIS2/src/SIS_debugging.F90 $(SRCROOT)/SIS2/src/combined_ice_ocean_driver.F90 $(SRCROOT)ice_param/ocean_rough.F90 $(SRCROOT)/SIS2/src/SIS_types.F90 $(SRCROOT)/SIS2/src/SIS_sum_output_type.F90 $(SRCROOT)/SIS2/src/SIS_hor_grid.F90 $(SRCROOT)/SIS2/src/SIS_transport.F90 $(SRCROOT)/SIS2/src/SIS2_ice_thm.F90 $(SRCROOT)/SIS2/src/SIS_fast_thermo.F90 $(SRCROOT)/SIS2/src/ice_type.F90 $(SRCROOT)ice_param/ice_albedo.F90 $(SRCROOT)/SIS2/src/SIS_optics.F90 $(SRCROOT)/SIS2/src/SIS_utils.F90 $(SRCROOT)/SIS2/src/ice_grid.F90 $(SRCROOT)/SIS2/src/SIS_dyn_cgrid.F90 $(SRCROOT)/SIS2/src/ice_shortwave_dEdd.F90 $(SRCROOT)/SIS2/src/SIS_sum_output.F90 $(SRCROOT)/SIS2/src/SIS_tracer_advect.F90 $(SRCROOT)/SIS2/src/SIS_dyn_trans.F90 $(SRCROOT)/SIS2/src/SIS_tracer_registry.F90 $(SRCROOT)/SIS2/src/SIS_dyn_bgrid.F90 $(SRCROOT)/SIS2/src/ice_ridge.F90 $(SRCROOT)/SIS2/src/SIS_diag_mediator.F90 $(SRCROOT)/SIS2/src/SIS_get_input.F90 $(SRCROOT)/SIS2/src/SIS_transcribe_grid.F90 $(SRCROOT)/SIS2/src/ice_thm.F90 $(SRCROOT)/SIS2/src/ice_boundary_types.F90 $(SRCROOT)/SIS2/src/SIS_fixed_initialization.F90 $(SRCROOT)/SIS2/src/ice_age_tracer.F90 $(SRCROOT)/SIS2/src/SIS_slow_thermo.F90 $(SRCROOT)/SIS2/src/SIS_tracer_flow_control.F90 $(SRCROOT)/SIS2/src/ice_spec.F90 $(SRCROOT)ice_param/ocean_albedo.F90 $(SRCROOT)/SIS2/src/SIS_continuity.F90 $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)MOM6/src//framework/version_variable.h
OBJ = SIS_ctrl_types.o ice_model.o SIS_debugging.o combined_ice_ocean_driver.o ocean_rough.o SIS_types.o SIS_sum_output_type.o SIS_hor_grid.o SIS_transport.o SIS2_ice_thm.o SIS_fast_thermo.o ice_type.o ice_albedo.o SIS_optics.o SIS_utils.o ice_grid.o SIS_dyn_cgrid.o ice_shortwave_dEdd.o SIS_sum_output.o SIS_tracer_advect.o SIS_dyn_trans.o SIS_tracer_registry.o SIS_dyn_bgrid.o ice_ridge.o SIS_diag_mediator.o SIS_get_input.o SIS_transcribe_grid.o ice_thm.o ice_boundary_types.o SIS_fixed_initialization.o ice_age_tracer.o SIS_slow_thermo.o SIS_tracer_flow_control.o ice_spec.o ocean_albedo.o SIS_continuity.o
OFF = $(SRCROOT)/SIS2/src/ice_model.F90 $(SRCROOT)/SIS2/src/SIS_slow_thermo.F90 $(SRCROOT)/SIS2/src/ice_thm.F90 $(SRCROOT)/SIS2/src/SIS_ctrl_types.F90 $(SRCROOT)/SIS2/config_src/dynamic_symmetric/SIS2_memory.h $(SRCROOT)/SIS2/src/SIS2_ice_thm.F90 $(SRCROOT)/SIS2/src/SIS_dyn_bgrid.F90 $(SRCROOT)/SIS2/src/ice_type.F90 $(SRCROOT)/SIS2/src/ice_spec.F90 $(SRCROOT)/SIS2/src/SIS_diag_mediator.F90 $(SRCROOT)/SIS2/src/SIS_dyn_cgrid.F90 $(SRCROOT)ice_param/ocean_rough.F90 $(SRCROOT)/SIS2/src/combined_ice_ocean_driver.F90 $(SRCROOT)/SIS2/src/ice_grid.F90 $(SRCROOT)/SIS2/src/SIS_get_input.F90 $(SRCROOT)ice_param/ice_albedo.F90 $(SRCROOT)/SIS2/src/SIS_utils.F90 $(SRCROOT)/SIS2/src/SIS_fixed_initialization.F90 $(SRCROOT)/SIS2/src/SIS_optics.F90 $(SRCROOT)/SIS2/src/SIS_hor_grid.F90 $(SRCROOT)/SIS2/src/SIS_sum_output.F90 $(SRCROOT)/SIS2/src/ice_age_tracer.F90 $(SRCROOT)/SIS2/src/SIS_debugging.F90 $(SRCROOT)/SIS2/src/SIS_tracer_registry.F90 $(SRCROOT)/SIS2/src/SIS_fast_thermo.F90 $(SRCROOT)/SIS2/src/SIS_tracer_advect.F90 $(SRCROOT)/SIS2/src/SIS_transcribe_grid.F90 $(SRCROOT)/SIS2/src/ice_boundary_types.F90 $(SRCROOT)/SIS2/src/SIS_dyn_trans.F90 $(SRCROOT)/SIS2/src/SIS_tracer_flow_control.F90 $(SRCROOT)/SIS2/src/SIS_continuity.F90 $(SRCROOT)/SIS2/src/SIS_memory_macros.h $(SRCROOT)MOM6/src//framework/MOM_memory_macros.h $(SRCROOT)/SIS2/src/SIS_types.F90 $(SRCROOT)/SIS2/src/ice_ridge.F90 $(SRCROOT)MOM6/src//framework/version_variable.h $(SRCROOT)/SIS2/src/ice_shortwave_dEdd.F90 $(SRCROOT)/SIS2/src/SIS_transport.F90 $(SRCROOT)ice_param/ocean_albedo.F90 $(SRCROOT)/SIS2/src/SIS_sum_output_type.F90
clean: neat
	$(RM) -f .libsis2.a.cppdefs $(OBJ) libsis2.a *.mod

neat:
	$(RM) -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libsis2.a: $(OBJ)
	$(AR) $(ARFLAGS) libsis2.a $(OBJ)

