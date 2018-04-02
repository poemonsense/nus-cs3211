echo "Cleaning..."
./clean.sh
echo "Compiling..."
./compile.sh
echo "Running program..."
mpirun --oversubscribe -np $1 ./pool initialspec.txt finalbrd.ppm
