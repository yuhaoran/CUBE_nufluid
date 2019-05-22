source ../utilities/module_load_mac.sh

cp parameters.f90_1 parameters.f90
source run_mac.sh > log1

cp parameters.f90_2 parameters.f90
source run_mac.sh > log2

cp parameters.f90_3 parameters.f90
source run_mac.sh > log3

cp parameters.f90_4 parameters.f90
source run_mac.sh > log4

cp parameters.f90_5 parameters.f90
source run_mac.sh > log5

