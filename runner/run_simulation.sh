#!/bin/bash

# The initial sentence allows the file to use the bash shell (Thus, the therminal command window)

# Set default input values: the general structure is var_name=${input_numer:-default_value}
Niter=${1:-3}             #Number of iterations
density=${2:-0.7}         #Density of particles
temp=${3:-1.9}            #Bath temperature of the system
Nsteps=${4:-50000}        #Number of steps for montecarlo
Nbins=${5:-50}            #Number of bins
Nbinsrdf=${6:-50}         #Number of bins for rdf calculation

# Compile particles.cpp
g++ particles.cpp -o particles_executable

# Create folders for the data
mkdir -p positions
mkdir -p energies

# Loop to save the data
# Here "$" performs a command substitution of "seq", 
# so the different values generated by "seq 1 $Niter" are stored in "i" as the loop runs
for i in $(seq 1 $Niter); do
    
    echo "Iteration $i ..............." # Printing command

    # Run the executable
    ./particles_executable $density $temp $Nsteps

    # Rename and save the files with "mv" command
    mv position.txt ./positions/position$i.txt
    mv energy.txt ./energies/energy$i.txt

    # Create a file of the last 300 energy data (Here the append operator ">>" is used)
    tail -300 "./energies/energy$i.txt" >> "energies/equilibrium_energies.txt"
    tail -100 "./positions/position$i.txt" >> "positions/equilibrium_positions.txt"
done

# Compile histogram.cpp
g++ histogram.cpp -o histogram_executable
# Run the histogram executable
./histogram_executable $Nbins

# Compile histogram.cpp
g++ rdf.cpp -o rdf_executable
# Run the histogram executable
./rdf_executable $Nbinsrdf

# Remove the executables
rm particles_executable
rm histogram_executable
rm rdf_executable

# Text output as help
echo "Execution completed successfully."
echo ""
echo "Variables to add:"
echo "  Number of iterations"
echo "  Density of particles"
echo "  Bath temperature of the system"
echo "  Number of steps for montecarlo"
echo "  Number of bins for energy histogram"
echo "  Number of bins for rdf representation"
