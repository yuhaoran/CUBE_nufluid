export FC='gfortran'
export XFLAG=' -O3 -cpp -fopenmp -fcoarray=single'
export XFLAG_NO_OMP='-O3 -cpp -fcoarray=single'
#export XFLAG=' -cpp -fopenmp -fcoarray=single -fcheck=all'
#export XFLAG='-O3 -cpp -fcoarray=single -fopenmp'
export OFLAG=${XFLAG}' -c'
export OFLAG_NO_OMP=${XFLAG_NO_OMP}' -c'

#export FFTFLAG='-I/usr/local/include/ -L/usr/local/lib/ -lfftw3f -lm -ldl'
export FFTFLAG='-I/Users/haoran/opt/include/ -L/Users/haoran/opt/lib/ -lfftw3f -lfftw3 -lfftw3_omp -lfftw3f_omp -lm -ldl'

export OMP_STACKSIZE=6000M
export OMP_NUM_THREADS=4
#export OMP_THREAD_LIMIT=4
ulimit -s 61000
ulimit
