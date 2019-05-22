source ../utilities/module_load_mac.sh

#rm -rf ../output/universe1/*
make clean
make

cd ../utilities/
make clean
make

./ic.x
#./ic_nu.x

cd ../main/
./main.x

cd ../utilities/
./cicrsd.x
./cicpower.x
./halofinder.x
./dsp.x
./ang_mom_corr.x
./lpt.x
#source findhalos.sh

cd ../main/
