ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

ifndef ARCH
$(error build.make: ARCH is not set)
endif

ifndef NDIM
$(error build.make: NDIM is not set)
endif

SRC_DIRS := . modules amrvacio physics rho hd mhd rhd particle nonlinear rd mf sir
SRC_DIRS := $(addprefix $(AMRVAC_DIR)/src/, $(SRC_DIRS))
LIB_AMRVAC := libamrvac.a
PPFLAGS :=

.PHONY: libamrvac clean
.PRECIOUS: %.f			# Don't remove intermediate .f files

libamrvac: $(LIB_AMRVAC) amrvac.f

# Include makefiles, which define FOBJECTS and dependencies
include $(addsuffix /makefile, $(SRC_DIRS))

# Include architecture and rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

# Get .t files from SRC_DIRS
vpath %.t $(SRC_DIRS)

OBJECTS := $(FOBJECTS:.t=.o) $(INCLUDES:.t=.o)

# # Include multigrid coupling
# ifneq ($(NDIM), 1)
# LIBOMGDIR := $(AMRVAC_DIR)/external_libs/octree-mg
# vpath %.f90 $(LIBOMGDIR)/src

# ifeq ($(NDIM), 3)
# vpath %.f90 $(LIBOMGDIR)/poisson_3d_fft
# include $(LIBOMGDIR)/poisson_3d_fft/definitions.make
# endif

# include $(LIBOMGDIR)/src/definitions.make
# include $(LIBOMGDIR)/makerules.make
# endif

$(LIB_AMRVAC): $(OBJECTS)
	$(RM) $@
	$(AR) rcs $@ $^

clean:
	$(RM) *.o *.mod *.f $(LIB_AMRVAC)

# INCLUDES are always compiled before FOBJECTS
$(FOBJECTS:.t=.o): $(INCLUDES:.t=.mod)
