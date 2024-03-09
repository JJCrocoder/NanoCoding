#!/bin/bash

# The initial sentence allows the file to use the bash shell (Thus, the therminal command window)

# Set default input values: the general structure is var_name=${input_numer:-default_value}
Niter=${1:-3}             #Number of iterations
density=${2:-0.7}         #Density of particles
temp=${3:-1.9}            #Bath temperature of the system
Nsteps=${4:-50000}        #Number of steps
Lbox=${7:-10}             #Default box length

# Compile brownian.cpp
g++ brownian_sim.cpp -o brownian_executable

# Create folders for the data
mkdir -p positions

# Run the executable
./brownian_executable

# Rename and save the files with "mv" command
mv position.txt ./positions/position.txt

# Remove the executables
rm brownian_executable

# Text output as help
echo "Execution completed successfully."
echo ""
echo "Variables to add:"
echo "  Number of iterations"
echo "  Density of particles"
echo "  Bath temperature of the system"
echo "  Number of steps for simulation"
echo "  Box length"
