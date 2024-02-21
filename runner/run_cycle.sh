#!/bin/bash

while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in 
    --density)
        density="$2"
        shift
        ;;
    --temp)
        temp="$2"
        shift 
        ;;
    --N)
        N="$2"
        shift
        ;;
    --max)
        max="$2"
        shift
        ;;
    *)
    esac
    shift
done


# Compile particles.cpp
g++ particle_metropolis.cpp -o particles_executable

#Execute the code
#Idealmente corremos varias veces este codigo almacenando los txt en un nuevo directorio
mkdir energy_files
i=0
while [[ $i -le $max ]]; do
    ./particles_executable "$density" "$temp"
    #change the name of the energy files and move it to a file so keep everything organize
    mv "./energy.txt" "./energy_{$i}.txt"
    mv "./energy_{$i}.txt" ./energy_files/
    #make a file with the last 300 data
    tail -300 "./energy_files/energy_{$i}.txt" >> "energy_files/equilibrium.txt"
    i=$((i+1))
done 

# Compile histogramN.cpp
g++ histogram.cpp -o histogram_executable -lm

# Run the histogram executable
./histogram_executable "./energy_files/equilibrium.txt" "$N"
