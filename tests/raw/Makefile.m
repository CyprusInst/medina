# Part of the MEDINA: MECCA - KPP Fortran to CUDA source-to-source pre-processor
# Minimal Makefile 
# ----------------------------------------------
SRCS0 := $(wildcard *.f90)
SRCS_ACC  := $(SRCS_ACC:.cu=.o) 
SRCS_ACC  := $(wildcard *.cu) 
SRCS  := $(filter-out F%.f90, $(SRCS0))
OBJSCUDA  := $(SRCS_ACC:.cu=.o)
OBJS  := $(SRCS:.f90=.o)
MODS  := $(SRCS:.f90=.mod)
MODS_INST := $(addprefix $(PREFIX)/include/, $(SRCS))

.SUFFIXES: $(SUFFIXES) .f90 .md5 .cu

%.o: %.f90
	$(F90) $(F90NOR8) $(INCLUDES) -c $< -o $@
# ----------------------------------------------------------------------

all: messy_main_compilerinfo_mem.f90 $(LIB)

$(LIB): depend $(OBJS) 
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)  
	$(AR) -dv $(LIB) messy_ncregrid_interface.o

# check files
list:
	@echo "------------------------------------------------"
	@echo "SRCS = $(SRCS)"
	@echo "------------------------------------------------"
	@echo
	@echo "------------------------------------------------"
	@echo "OBJS = $(OBJS)"
	@echo "------------------------------------------------"

depend $(MAKEFILE_INC): $(SRCS) 
	$(F_makedepend)


