#
# Use `FC=<compiler> ... make` to select compiler
# Use `DEBUG=1 ... make` for a debug build
# Use `NC=1 ... make` to build with NetCDF using `nf-config`

# Compiler
FC ?= gfortran
ifeq ($(FC), f77)  # override possible Make default
  FC := gfortran
endif
$(info FC setting: '$(FC)')

# Default to non-debug build
DEBUG ?= 0

# Default to NetCDF build
NC ?= 1

# Default to first `nf-config` on PATH
NF_CONFIG ?= nf-config

# Compile flags
$(info DEBUG setting: '$(DEBUG)')
ifeq ($(DEBUG), 0)
  FCFLAGS := -O3
else ifeq ($(DEBUG), 1)
  ifeq ($(findstring gfortran,$(notdir $(FC))),gfortran)
    FCFLAGS := -g -Og -Wall -Wextra -Wconversion -pedantic \
      -frecursive -fcheck=bounds -fall-intrinsics -fmax-errors=5 \
      -std=f2003
  else ifeq ($(FC),ifort)
    FCFLAGS := -g -O0 -warn all -check bounds -implicitnone -error-limit 5
  endif
else ifeq ($(DEBUG), 2)
  ifeq ($(findstring gfortran,$(notdir $(FC))),gfortran)
    FCFLAGS := -g -Og -Wall -Wextra -Wconversion -pedantic \
      -frecursive -fcheck=all -fall-intrinsics -fmax-errors=0 \
      -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999 \
      -std=f2003
  else ifeq ($(FC),ifort)
    FCFLAGS := -g  -O0 -warn all -check all -implicitnone -error-limit 0 \
      -fpe0 -traceback
  endif
else
  $(error invalid setting for DEBUG, should be 0, 1 or 2 but is '$(DEBUG)')
endif
$(info FCFLAGS: '$(FCFLAGS)')

# NETCDF Settings here
LIBS :=
INC :=
$(info NC setting: '$(NC)')
$(info NF_CONFIG setting: '$(NF_CONFIG)')
ifeq ($(NC), 1)
  NETCDF_FLIBS := $(shell $(NF_CONFIG) --flibs)
  NETCDF_INC := -I$(shell $(NF_CONFIG) --includedir)
  #
  LIBS += $(NETCDF_FLIBS)
  INC += $(NETCDF_INC)
  FCFLAGS += -DNETCDF
else ifeq ($(NC), 0)
  #
else
  $(error invalid setting for NC, should be 0 or 1 but is '$(NC)')
endif
$(info LIBS: '$(LIBS)')
$(info INC:  '$(INC)')

# Source objects
OBJS := \
  canopy_const_mod.o \
  canopy_date_mod.o \
  canopy_coord_mod.o \
  canopy_canopts_mod.o \
  canopy_canmet_mod.o \
  canopy_canvars_mod.o \
  canopy_utils_mod.o \
  canopy_files_mod.o \
  canopy_readnml.o \
  canopy_alloc.o \
  canopy_init.o \
  canopy_txt_io_mod.o \
  canopy_ncf_io_mod.o \
  canopy_check_input.o \
  canopy_read_txt.o \
  canopy_dxcalc_mod.o \
  canopy_profile_mod.o \
  canopy_var3din_mod.o \
  canopy_phot_mod.o \
  canopy_rad_mod.o \
  canopy_tleaf_mod.o \
  canopy_wind_mod.o \
  canopy_fire_mod.o \
  canopy_eddy_mod.o \
  canopy_bioparm_mod.o \
  canopy_bioemi_mod.o \
  canopy_drydep_mod.o \
  canopy_calcs.o \
  canopy_write_txt.o \
  canopy_dealloc.o \
  canopy_app.o

ifeq ($(NC), 0)
  _ncf_objs := canopy_check_input.o canopy_ncf_io_mod.o
  OBJS := $(filter-out $(_ncf_objs),$(OBJS))
endif

# Program name
PROGRAM := ../canopy

# Targets
.PHONY: all clean
all: $(PROGRAM)

$(PROGRAM): $(OBJS_APP) $(OBJS)
	$(FC) $(FCFLAGS) $^ -o $@ $(LIBS) $(INC)

%.o: %.F90
	$(FC) $(FCFLAGS) $(INC) -c $<

clean:
	rm -f *.o *.mod $(PROGRAM)
