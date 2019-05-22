module purge
module load intel/16.0.3 intelmpi/5.0.3.048 fftw/3.3.4-intel-impi
module load caf/intel/any
#module list

export FC='caf'
export XFLAG='-O3 -fpp -mcmodel=medium'

#export FC='mpiifort'
#export XFLAG='-O3 -ilp64 -fpp -mcmodel=medium -coarray=distributed'

# GPC nodes do not support intel AVX instructions "-xHost"
#export XFLAG='-O3 -fpp -mcmodel=medium'

export OFLAG=${XFLAG}' -c'
export FFTFLAG='-I'${SCINET_FFTW_INC}' ''-L'${SCINET_FFTW_LIB}' -lfftw3f -lm -ldl'

# use uniformized caf on GPC
# export FOR_COARRAY_NUM_IMAGES=1
