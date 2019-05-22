module purge
module load gcc/5.2.0 openmpi/gcc/1.8.3 use.experimental caf/gcc/5.2.0-openmpi
module load fftw/3.3.0-gcc-openmpi

#module load gcc/6.2.0  openmpi/gcc6/1.10.6 use.experimental fftw/3.3.0-gcc-openmpi

#module list

export FC='caf'
export XFLAG='-O3 -cpp -march=native -mcmodel=medium'
export OFLAG=${XFLAG}' -c'
export FFTFLAG='-I'${SCINET_FFTW_INC}' ''-L'${SCINET_FFTW_LIB}' -lfftw3f -lm -ldl'
