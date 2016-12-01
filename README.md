# MECCA - KPP Fortran to CUDA source-to-source pre-processor
# README.txt

** Disclaimer: This software is in alpha-test mode, 
equivalent to the MESSy red traffic light status.
No unexpected behaviour was observed under testing, and users are 
invited to test with their model setup. However, no express guarantee
is provided for production simulations. 
For assistance or to report problems please contact the maintainers:
christoudias@cyi.ac.cy; m.alvanos@cyi.ac.cy

1. Requirements:

Software: CUDA compiler and python are required for the pre-processor. 
Hardware: CUDA compatible GPU (NVIDIA). 

2. Installation:

There are two files required to enable using the GPUs: 
f2c_alpha.py  and kpp_integrate_cuda_prototype.cu. 

The files have to be available in the messy/util directory. 
No additional changes are required. 

Note: MESSy has to be linked with the -lcudart flag. 
For example, you can append it to the SPEC_NETCDF_LIB variable 
in the configuration file (under config/mh-XXXX).

3. Running the MECCA Fortran to CUDA source-to-source pre-processor:

You have to enter the ./messy/util directory to execute the
preprocessor, by running "python f2c_alpha.py". The preprocessor expects
the following files to be in place:

 messy/smcl/messy_mecca_kpp.f90
 messy/smcl/messy_cmn_photol_mem.f90
 messy/smcl/messy_main_constants_mem.f90
 messy/util/kpp_integrate_cuda_prototype.cu
 messy/smcl/specific.mk
 messy/smcl/Makefile.m

If any of these files is missing or not configured as in the MESSy release,
the preprocessor will stop with an error message.

4. Running EMAC with GPU MECCA and improving performance:

The runtime parameter NPROMA should be set to a value not greater than 128.
This allows for optimal memory allocation and performance on the GPU.

Each CPU process that offloads to GPU requires a chunk of the GPU VRAM memory,
dependent on the number of species and reaction constants in the MECCA mechanism. 
The number of GPUs per node and VRAM memory available in each GPU dictates the
total number of CPU cores that can run simultaneously.

Warning: When running multiple CPU processes per GPU, if
memory is not enough the CUDA runtime will fail silently - without any
error. A solution in that case is to use the Multi-process service provided
by NVIDIA as an alternative.

During experiments with an engineering sample of the next generation 
NVIDIA Pascal architecture, the source application will fail due to 
large local memory requirements. Transforming runtime GPU local access to global 
solves the problem, at a performance cost.

*EOF*
