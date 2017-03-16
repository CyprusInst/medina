
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

echo "====> STEP 5: Comparing the output results..."

