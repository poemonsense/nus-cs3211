echo "Cleaning..."
./clean.sh
mkdir ppmresults
mkdir logs
echo "Compiling..."
./compile.sh
echo "Running program..."
mpirun --oversubscribe -np $1 ./pool initialspec.txt ppmresults/finalbrd
