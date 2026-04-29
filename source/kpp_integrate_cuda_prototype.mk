# -*- Makefile -*-

# - drop-in include file for MECCA/KPP CUDA integration ------------------------

# - sources --------------------------------------------------------------------
 SRCS_CUDA := $(wildcard *.cu) 
 OBJS_CUDA := $(SRCS_CUDA:.cu=.o)

 SRCS += $(SRCS_CUDA)
 OBJS += $(OBJS_CUDA)

.SUFFIXES: $(SUFFIXES) .cu

# uncomment to clean up object files and ensure re-compilation
# beware: it may cause recursive loop of the dependency checker
# $(shell touch $(SRCS_CUDA))
# $(shell rm $(OBJS_CUDA))

# - parameters & flags that we do MECCA/KPP integration with CUDA --------------
# KPP_CUDA becomes effective in messy_mecca_kpp.f90
 CUDA_DEFS = KPP_CUDA #DEBUG

# - GPU architecture and CUDA options ------------------------------------------
# will be substituted, e.g. CUDA_ARCH = --gpu-architecture=sm_60
 CUDA_ARCH = <CUDA_ARCH>
# will be substituted, e.g. NVCCFLAGS = --ptxas-options=-v --ftz=false --prec-div=true --prec-sqrt=true --fmad=false
# uncomment #-DDEBUG to enable CUDA code debug
 NVCCFLAGS = <NVCCFLAGS> #-DDEBUG

# - nvcc compiler call -------------------------------------------------
%.o: %.cu specific.mk 
	nvcc --verbose $(CUDA_ARCH) $(NVCCFLAGS) -O3 -g -c $<

# - adjusting build options ----------------------------------------------------
# add CUDA_DEFS for compilation of MECCA KPP (excludes other KP4-processed code, e.g. SCAV)
messy_mecca_kpp.o: F90NOR8 += $(addprefix $(DEFOPT),$(CUDA_DEFS))
