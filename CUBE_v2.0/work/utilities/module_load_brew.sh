export FC='gfortran'
export XFLAG=' -O3 -cpp -fopenmp -fcoarray=single -mcmodel=medium'
export XFLAG_NO_OMP='-O3 -cpp -fcoarray=single'
#export XFLAG=' -cpp -fopenmp -fcoarray=single -fcheck=all'
#export XFLAG='-O3 -cpp -fcoarray=single -fopenmp'
export OFLAG=${XFLAG}' -c'
export OFLAG_NO_OMP=${XFLAG_NO_OMP}' -c'

#export FFTFLAG='-I/usr/local/include/ -L/usr/local/lib/ -lfftw3f -lm -ldl'
export FFTFLAG='-I/usr/local/Cellar/fftw/3.3.8_1/include/ -L/usr/local/Cellar/fftw/3.3.8_1/lib/ -lfftw3f -lm -ldl'

export OMP_STACKSIZE=6000M
export OMP_NUM_THREADS=12
#export OMP_THREAD_LIMIT=4
ulimit -s 64000
ulimit
