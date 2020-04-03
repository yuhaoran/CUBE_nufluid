source module_load_brew.sh
make clean
make
cd ../main
make clean
make

cd ../utilities
./ic.x
./ic_nu.x

cd ../main
./main.x

cd ../utilities
./cicpower.x

