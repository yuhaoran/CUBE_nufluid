module load gcc openmpi/3.1.0rc3 fftw/3.3.7 

#export XFLAG=' -O3 -cpp -fopenmp -fcoarray=single -mcmodel=medium'
export XFLAG=' -cpp -fopenmp -fcoarray=lib'
#export XFLAG='-O3 -cpp -fcoarray=single -fopenmp'
export OFLAG=${XFLAG}' -c'
export FFTFLAG='-I${SCINET_FFTW_ROOT}/include/ -L${SCINET_FFTW_ROOT}/lib/ -lfftw3f -lm -ldl'
# -fopenmp cause (maybe memory) probelm: Segmentation fault: 11
# in cumsum6

export OMP_STACKSIZE=1000M
export OMP_NUM_THREADS=4
export OMP_THREAD_LIMIT=4
export FOR_COARRAY_NUM_IMAGES=1
