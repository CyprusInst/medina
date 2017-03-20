
####
#### Testing driver for the f2c script.
####
##########################################


echo "==============================="
echo "==============================="
echo "==============================="
echo "====> STEP 1: Copying files..."

(
set -x
cp raw/* ./messy/smcl/
cp ../f2c_alpha.py  ./messy/util/
cp ../kpp_integrate_cuda_prototype.cu  ./messy/util/
)

echo "====> STEP 2: Running script..."

(
set -x
cd messy/util
python ./f2c_alpha.py > /dev/null
)

status=$?
if [ $status == 1 ]; then
       echo "Python parser - Unsuccessful"
       exit -1
fi

echo "====> STEP 3: Compiling the output files..."

(
set -x
cd messy/smcl
nvcc -O0 -c messy_mecca_kpp_acc.cu  2>&1  | grep error
)

status=$?
if [ $status == 0 ]; then
       echo "NVCC - Unsuccessful"
       exit -1
fi

echo "====> STEP 4: Running the application..."


(
set -x
cat ./raw/main.c >> ./messy/smcl/messy_mecca_kpp_acc.cu
cd messy/smcl
nvcc -O1  messy_mecca_kpp_acc.cu  2>&1  | grep error
./a.out
cuda-memcheck ./a.out
)

status=$?
if [ $status == 0 ]; then
       echo "NVCC - Unsuccessful"
       exit -1
fi



echo "====> STEP 5: Compiling original version in FORTRAN..."


(
set -x
cd messy/fortran
gfortran -c messy_cmn_photol_mem.f90
gfortran -c messy_main_constants_mem.f90
gfortran -c messy_mecca_kpp.f90
gfortran -c main.f90
./a.out
)

echo "====> STEP 6: Comparing the output results..."

echo "====> STEP 7: Cleaning up the directories..."


(
set -x
cd messy/smcl/
rm ./*
cd ../util/
rm ./*
)



echo "====> Testing Finished"


