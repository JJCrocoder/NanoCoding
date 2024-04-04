#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>
using namespace std;

// Periodic boundary conditions (Minimum image convention)
void mic(float * vec, float Lbox) {
    for (int i = 0; i < 3 ; ++i) vec[i] -= floorf(0.5 + vec[i]/Lbox)*Lbox;
    return;
}

// Function prototypes
void moveEuler(float Pos[], float For[]);
void forces(float Pos[], float For[]);

// Global variables
const float sigma = 1.0;
const float epsilon = 1.0;
const int dim = 3;
float Lbox = 10.0;
float volume = Lbox*Lbox*Lbox;
float r_cut = 2.5 * sigma;

int nsamp_rdf = 10;
int nsamp_pos = 100;
float DIFF =  1.0;
float deltaT = 0.005;
float density = 0.1;
float Temp = 1.9;
float Kb = 1/Temp;
int Npart=(int)(density*volume);
int Nstep = 50000;
float mu = DIFF / (Kb*Temp);

// Main function
int main(int argc, char *argv[]) {

    // Variables
    float Pos[dim*Npart];
    float For[dim*Npart]={0};

    // LOG of the run
    cout << "Particles = "<< Npart <<endl;    

    // Open data files
    ofstream fich_pos;
    fich_pos.open("motion_brownian.txt");

    // Generate initial positions
    float pos;
    vector<float> x;
    fstream positions("equilibrium_positions.txt");
    while (positions >> pos) x.push_back(pos);
    positions.close();
    for (int i = 0; i < 3*Npart; ++i) Pos[i] = x[i];
    cout << "...Initial positions generated" <<endl;

    // Brownian Dynamics loop
    for (int istep = 0; istep < Nstep+1; ++istep) {

        // Save positions
        for (int i = 0; i < Npart; ++i) {
            for (int k = 0; k < dim; ++k) fich_pos << Pos[i * dim + k] << " ";
            fich_pos << endl;
        }

        // Calculate forces
        forces(Pos, For);

        // Move Particles
        moveEuler(Pos, For);
    }

    // Close data files
    fich_pos.close();

    return 0;
}

// Euler algorithm function for the positional motion of particles
void moveEuler(float Pos[], float For[]){
    // Generate gaussian random numbers
    random_device rand_dev;
    mt19937 generator(rand_dev());
    normal_distribution<double> distribucion_normal(0, 1);
    
    // Loop for update each particle
    for (int i = 0; i < Npart; ++i) {
        float nu = distribucion_normal(generator);

        // Update position using Euler method
        for (int k = 0; k < dim; ++k) {
            Pos[i * dim + k] += mu * For[i * dim + k] * deltaT * sqrt(2.0*Kb*Temp*mu*deltaT) * nu * Lbox;

            // Periodic boundary conditions (Minimum image convention)
            Pos[i * dim + k] -= floorf(0.5 + Pos[i * dim + k]/Lbox)*Lbox;
        }
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

            // Lennard-Jones forces if we are under the cut
            if (r2_ij < r_2cut) {
                float r_mod = sqrt(r2_ij);
                if (r_mod < 1.2*sigma) continue;
                float r6 = pow((sigma / r_mod), 6.0);

                // Calculate the force as the derivative of the Lennard-Jones energy
                float f_ij = -48 * epsilon * r6 * (r6 - 0.5) / r_mod;

                // Calculate the force components for both particles (3rd Newton Law)
                for (int k = 0; k < dim; ++k) {
                    For[j * dim + k] += (f_ij * r_ij[k]) / r_mod;
                    For[i * dim + k] -= (f_ij * r_ij[k]) / r_mod;
                }
            }
        }
    }
}
