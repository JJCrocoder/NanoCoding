#!/bin/bash

# Set default input values
Niter=${1:-3}             #Number of iterations
density=${2:-0.7}         #Densisy of particles
temp=${3:-1.9}            #Bath temperature of the system
Nsteps=${4:-50000}        #Number of steps for montecarlo
Nbins=${5:-50}            #Number of bins

# Compile particles.cpp
g++ particles.cpp -o particles_executable

# Create folders for the data
mkdir -p positions
mkdir -p energies

# Loop for save the data
for i in $(seq 1 $Niter); do
    echo "Iteration $i ..............."

    # Run the executable
    ./particles_executable $density $temp $Nsteps

    # Rename and save the files
    mv position.txt ./positions/position$i.txt
    mv energy.txt ./energies/energy$i.txt

    # Create a file of the last 300 energy data
    tail -300 "./energies/energy$i.txt" >> "energies/equilibrium_energies.txt"
done

# Compile histogram.cpp
g++ histogram.cpp -o histogram_executable

# Run the histogram executable
./histogram_executable $Nbins

# Remove the executables
rm particles_executable
rm histogram_executable

# Text
echo "Execution completed successfully."
echo ""
echo "Variables to add:"
echo "  Number of iterations"
echo "  Densisy of particles"
echo "  Bath temperature of the system"
echo "  Number of steps for montecarlo"
echo "  Number of bins"