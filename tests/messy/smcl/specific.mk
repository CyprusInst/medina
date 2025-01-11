
# Part of the MEDINA: MECCA - KPP Fortran to CUDA source-to-source pre-processor

# Linux64 / Intel Compiler ###################################################
ifeq ($(ARCH), LINUX64)
ifeq ($(COMPILER), INTEL)

messy_mecca.o: messy_mecca.f90
	$(F90) $(F90NOR8) -nocheck -c $<
messy_mecca_kpp.o: messy_mecca_kpp.f90
	$(F90) $(F90NOR8) -nocheck -c $<
endif
endif
##############################################################################




messy_mecca_kpp_acc.o: messy_mecca_kpp_acc.cu specific.mk 
	nvcc  -v  --ptxas-options=-v  --gpu-architecture=compute_70  --ftz=false --prec-div=true --prec-sqrt=true --fmad=false -O3  -g   -c  $<