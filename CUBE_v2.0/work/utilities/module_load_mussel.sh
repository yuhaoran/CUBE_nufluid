# for SciNet niagara: module load NiaEnv/2018a intel/2018.2 intelmpi/2018.2 fftw/3.3.7

# on CITA, load modules:
# module load intel/intel-18 intelmpi/2018.1.163 fftw/3.3.7-intelmpi-18
# for kingcrab, lobster, homard etc.,
# also
# export I_MPI_FABRICS=shm:tcp
# export I_MPI_TCP_NETMASK=10.5.0.0
# bugfix: if have problems compiling lpt.x, remove the -qopenmp flag

#export XFLAG=' -O3 -cpp -fopenmp -fcoarray=single -mcmodel=medium'
export FC='ifort'
export XFLAG='-O3 -fpp -qopenmp -coarray=single -mcmodel=medium'
export XFLAG_NO_OMP='-O3 -fpp -coarray=single -mcmodel=medium'
#export XFLAG='-O3 -fpp -coarray=single -fopenmp'
export OFLAG=${XFLAG}' -c'
export OFLAG_NO_OMP=${XFLAG_NO_OMP}' -c'
#export FFTFLAG='-I$(SCINET_FFTW_ROOT)/include/ -L$(SCINET_FFTW_ROOT)/lib/ -lfftw3f -lm -ldl'
export FFTFLAG='-lfftw3_omp -lfftw3f_omp -lm -ldl'
# -fopenmp cause (maybe memory) probelm: Segmentation fault: 11
# in cumsum6

export OMP_STACKSIZE=64G
#export KMP_STACKSIZE=64000M
export OMP_NUM_THREADS=64
#export OMP_THREAD_LIMIT=4
export FOR_COARRAY_NUM_IMAGES=1
#ulimit -s 61000
