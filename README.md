# MECCA - KPP Fortran to CUDA source-to-source pre-processor

*Disclaimer: This software is in alpha-test mode, 
equivalent to the MESSy red traffic light status.
No unexpected behaviour was observed under testing, and users are 
invited to test with their model setup. However, no express guarantee
is provided for production simulations. 
For assistance or to report problems please contact:* christoudias@cyi.ac.cy
 
## 1. Requirements:

Software: CUDA compiler and python are required for the pre-processor. 

Hardware: CUDA compatible GPU (Fermi, Kepler, Pascal, Volta, or later). 

## 2. Installation:

To enable using the GPUs the file: 
`xmedina.py`  and folder `source` 
have to be available in the messy/util/medina directory. 
No additional changes are required. 

### Command-line options

The following command line options are available to the user
(and can be used for example to run `xmedina.py` in batch mode):

* `-r / --ros`  An integer value of the Rosenbrock solver produced [1: all (select at runtime), 2: Ros2, 3: Ros3, 4: Rodas3, 5: Rodas4]
* `-g / --gpu`  An integer value of the architecture [1: FERMI, 2: KEPLER, 3: MAXWELL, 4: PASCAL]
* `-m / --mem`  Choose y or n wheter you want to use the global memory version (reduced stack allocation and performance) or not
* `-s / --smcl` MESSy smcl folder location, default: "../../smcl/"'


## 3. Running the MECCA Fortran to CUDA source-to-source pre-processor:

You have to enter the ./messy/util/medina directory to execute the
preprocessor, by running "`python xmedina.py`". The preprocessor expects
the following files to be in place:

* `messy/smcl/messy_mecca_kpp.f90`
* `messy/smcl/messy_cmn_photol_mem.f90`
* `messy/smcl/messy_main_constants_mem.f90`
* `messy/util/medina/source/kpp_integrate_cuda_prototype.cu`
* `messy/smcl/specific.mk`
* `messy/smcl/Makefile.m`
 
If any of these files is missing or not configured as in the MESSy release,
the preprocessor will stop with an error message.

## 4. Running EMAC with GPU MECCA and improving performance:

During testing it was found that the runtime parameter `NPROMA` should be set 
to a value not greater than 128 (preferably 64) for optimal memory allocation 
and performance on the GPU.

Each CPU process that offloads to GPU requires a chunk of the GPU VRAM memory,
dependent on the number of species and reaction constants in the MECCA mechanism. 
The number of GPUs per node and VRAM memory available in each GPU dictates the
total number of CPU cores that can run simultaneously.

### NVIDIA Multi-Process Service
To run multiple CPU processes per GPU, the Multi-process service (MPS) provided 
by NVIDIA should be used.

***Warning: Memory Protection*** 

Volta MPS client processes have fully isolated GPU address spaces. Pre-Volta MPS client 
processes allocate memory from different partitions of the same GPU virtual address space. As a result:

* An out-of-range write in a CUDA Kernel can modify the CUDA-accessible memory state of
another process, and will not trigger an error.
* An out-of-range read in a CUDA Kernel can access CUDA-accessible memory modified by 
another process, and will not trigger an error, leading to undefined behavior.

## 5. Unit testing

A self-contained unit test is included in the ditribution. The test includes 
reference source files implementing a simplified chemistry mechanism and 
compiles, exexutes and compares the FORTRAN (using gfortran) 
and auto-generated CUDA versions.

The test is executed by sourcing `driver.sh` under the `tests` directory. 
A utility script that compares the test solver output is also included in `tests/compare.py`

## 6. References

Alvanos, M. and Christoudias, T.: GPU-accelerated atmospheric chemical kinetics in the ECHAM/MESSy (EMAC) Earth system model (version 2.52), Geosci. Model Dev., 10, 3679-3693, https://doi.org/10.5194/gmd-10-3679-2017, 2017. 

Alvanos, M. and Christoudias, T., 2017. MEDINA: MECCA Development in Accelerators – KPP Fortran to CUDA source-to-source Pre-processor. Journal of Open Research Software, 5(1), p.13. DOI: http://doi.org/10.5334/jors.158

M. Alvanos and T. Christoudias, "Accelerating Atmospheric Chemical Kinetics for Climate Simulations," in IEEE Transactions on Parallel and Distributed Systems. DOI: 10.1109/TPDS.2019.2918798
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8723098&isnumber=4359390

