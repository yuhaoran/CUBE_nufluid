module purge
module load gcc/6.3.0 openmpi/gnu/2.0.2 fftw/intel/3.3.6-pl2

export FC='gfortran'
export XFLAG=' -O3 -cpp -fopenmp -fcoarray=single -mcmodel=medium'
#export XFLAG=' -cpp -fopenmp -fcoarray=lib'
export OFLAG=${XFLAG}' -c'
export FFTFLAG='-I${FFTW_ROOT}/include/ -L${FFTW_ROOT}/lib/ -lfftw3f -lm -ldl'
# -fopenmp cause (maybe memory) probelm: Segmentation fault: 11
# in cumsum6

export OMP_STACKSIZE=1000M
export OMP_NUM_THREADS=8
export OMP_THREAD_LIMIT=8
export FOR_COARRAY_NUM_IMAGES=1
