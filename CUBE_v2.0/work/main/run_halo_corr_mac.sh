source ../utilities/module_load_mac.sh

cd ../utilities/
make clean
make

cd ../main/
make clean
make EXTRA=-DHALOFIND

cd ../utilities/
./ic.x
./ic_nu.x

cd ../main/
./main.x # with runtime halofinder

cd ../utilities/
./cicpower.x
#./dsp.x
./ang_mom_corr.x # generate halo_init_spin

cd ../main/
#make clean
#make EXTRA=-Dhalo_spin_correlation
#./main.x

cd ../utilities/
#./ang_mom_corr.x

cd ../main/
