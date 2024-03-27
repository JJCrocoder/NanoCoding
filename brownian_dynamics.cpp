#include <cstdlib>
#include <random>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

typedef enum {VERLETSTEP1, VERLETSTEP2} verletstep;

// Periodic boundary conditions (Minimum image convention)
void mic(float * vec, float Lbox) {
    for (int k = 0; k < 3 ; ++k) vec[k] -= floorf(0.5 + vec[k]/Lbox)*Lbox;
    return;
}

// Function prototypes
void generateInitialVelocities(float Vel[]);
void moveVerlet(verletstep vstep, float Pos[], float Vel[], float For[]);
void forces(float Pos[], float For[]);
float temperature(float Vel[]);

// Global variables
const float sigma = 1.0;
float r_cut = 2.5 * sigma;
const float epsilon = 1.0;
const int dim = 3;
float Lbox = 10.0;
float volume = Lbox*Lbox*Lbox;

int Nstep = 100000;
float deltaT = 0.005;
float density = 0.1;
int Npart=(int)(density*volume);

// Main function
int main(int argc, char *argv[]) {

    // Variables
    float Pos[dim*Npart];
    float Vel[dim*Npart]={0};
    float For[dim*Npart]={0};    

    // LOG of the run
    cout << "Particles = "<< Npart <<endl;    

    // Open data files
    ofstream fich_pos;
    fich_pos.open("brownian_positions.txt");

    // Generate initial positions
    float pos;
    vector<float> x;
    fstream positions("equilibrium_energies.txt");
    while (positions >> pos) x.push_back(pos);
    positions.close();
    for (int i = 0; i < 3*Npart; ++i) Pos[i] = x[i];
    cout << "...Initial positions generated" <<endl;

    // Generate initial velocities
    generateInitialVelocities(Vel);
    cout << "...Initial velocities generated" <<endl;

    // Generate initial forces
    forces(Pos, For);
    cout << "...Initial forces generated" <<endl;

    // Print intial temperature
    cout << "Initial temperature: " << temperature(Vel) <<endl;

    // Molecular Dynamics loop
    for (int istep = 0; istep < Nstep+1; ++istep) {

        // Save positions
        for (int i = 0; i < Npart; ++i) {
            for (int k = 0; k < dim; ++k) fich_pos << Pos[i * dim + k] << " ";       
            fich_pos << endl;
        }

        // Move particles using Verlet 1
        moveVerlet(VERLETSTEP1, Pos, Vel, For);

        // Calculate forces
        forces(Pos, For);

        // Move particles using Verlet 2
        moveVerlet(VERLETSTEP2, Pos, Vel, For);

        // Print temperature
        if (istep % 1000 == 0) cout << "Temperature: " << temperature(Vel) <<endl;
    }

    // Close data files
    fich_pos.close();

    return 0;
}

// Generate the gaussian initial velocities
void generateInitialVelocities(float Vel[]) {
    random_device rand_dev;
    mt19937 generator(rand_dev());
    float temp = 1.9;
    uniform_real_distribution<> gauss(0.0, temp);
    
    // Generate gaussian velocities
    for (int i = 0; i < 3*Npart; ++i) Vel[i] = gauss(generator);
}

// Verlet algorithm for update particles position
void moveVerlet(verletstep vstep, float Pos[], float Vel[], float For[]){
    switch (vstep) {
    case VERLETSTEP1:
        // Loop for update each particle
        for (int i = 0; i < Npart; ++i) { 
            for (int k = 0; k < dim; ++k) {
                // Update position using Verlet method
                Pos[i * dim + k] += Vel[i * dim + k]*deltaT+ For[i * dim + k]*deltaT*deltaT*0.5;
                Vel[i * dim + k] += For[i * dim + k] * deltaT * 0.5;

                // Apply periodic boundary conditions
                Pos[i * dim + k] -= floorf(0.5 + Pos[i * dim + k]/Lbox)*Lbox;
            }
        }
    break;
    case VERLETSTEP2:
        // Loop for update each particle
        for (int i = 0; i < Npart; ++i) { 
            for (int k = 0; k < dim; ++k) {
                // Update position using Verlet method
                Vel[i * dim + k] += For[i * dim + k] * deltaT * 0.5;
                
                // Apply periodic boundary conditions
                Pos[i * dim + k] -= floorf(0.5 + Pos[i * dim + k]/Lbox)*Lbox;
            }
        }
    break;
    }
    return;
}

// Algorithm for update particles forces
void forces(float Pos[], float For[]) {
    float r_2cut = r_cut * r_cut;

    // Set the forces to zero
    for (int i = 0; i < 3*Npart; ++i) For[i] = 0;

    // Loop for all particles
    for (int i = 0; i < Npart-1; ++i) {
        for (int j = i+1; j < Npart; ++j) {
            float r_ij[dim];
            float r2_ij = 0.0;

            // Obtain the distances r_ij between i-j
            for (int k = 0; k < dim; ++k) r_ij[k] = (Pos[j * dim + k] - Pos[i * dim + k]);
            // Apply periodic boundary conditions
            mic(r_ij, Lbox);
            // Obtain the square module
            for (int k = 0; k < dim; ++k) r2_ij += pow(r_ij[k], 2);

            // Lennard-Jones forces (if we are under the cut)
            if (r2_ij < r_2cut) {
                float r_mod = sqrt(r2_ij);
                if (r_mod < 1.2*sigma) continue; // If superposition
                float r6 = pow((sigma / r_mod), 6.0);

                // Calculate the force as the derivative of the Lennard-Jones energy
                float f_ij = -48 * epsilon * r6 * (r6 - 0.5) / r_mod;

                // Calculate the force components for both particles (3rd Newton Law)
                for (int k = 0; k < dim; ++k) {
                    For[j * dim + k] -= (f_ij * r_ij[k]) / r_mod;
                    For[i * dim + k] += (f_ij * r_ij[k]) / r_mod;
                }
            }
        }
    }
}

// Temperature calculation
float temperature(float Vel[]){
    float temp=0.;
    for (int i = 0; i < 3*Npart; ++i) temp += Vel[i]*Vel[i];
    temp /= 3*Npart;
    return temp;
}
