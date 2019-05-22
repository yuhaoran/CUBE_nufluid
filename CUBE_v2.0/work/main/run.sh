source ../utilities/module_load_mac.sh 

cd ../utilities/
make clean
make ic.x cicpower.x
./ic.x

cd ../main/
make clean
make
./main.x

cd ../utilities/
./cicpower.x

cd ../main/
