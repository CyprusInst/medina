
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



