echo "Cleaning..."
./clean.sh
mkdir ppmresults
mkdir logs
echo "Compiling..."
./compile.sh -D POOL_DEBUG
echo "Dumping into debug/* ..."
rm -r debug
mkdir debug
objdump -D pool > debug/dump
echo "Running program..."
mpirun --oversubscribe -np $1 ./pool initialspec.txt ppmresults/finalbrd
