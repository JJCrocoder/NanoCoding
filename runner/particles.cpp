// LIBRARIES INCLUDED

// <iostream>: input and output operations
// <fstream>: file manipulation
// <cmath>: Mathematical functions and constants
// <random>: classes and functions for generating random numbers and sampling from different pdf

#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

// This allow to use declarations in "std" namespace without calling it
// Most of the variables or operations tha we declarate are in this namespace, so this line is convenient
using namespace std;

// Minimum image convention (periodic boundary conditions)
// Void functons doesn't return anything but you can change an argument inside them
// In this case, vec is created outside the function, but modified inside of it
void mic(float * vec, float Lbox) {
    for (int i = 0; i < 3 ; ++i) {
        vec[i] -= floorf(0.5 + vec[i]/Lbox)*Lbox;
    }  
    return;
}

// Function definition
float Energy(float Position[], float N_O_Pos[], int itag, int Npart);

// Global variables
float sigma = 1.0;
float epsilon = 1.0;
float Lbox = 10.0;
float volume = Lbox*Lbox*Lbox;
int dim = 3;
float r_cut = 2.5 * sigma;

// rdf variables
float dmax = 0.5 * Lbox;
float dmin = 0.0;

// Main function
int main(int argc, char *argv[]) {
    random_device rand_dev;
    mt19937 generator(rand_dev());

    // Variables
    float densinput = atof(argv[1]);
    float Temp = atof(argv[2]);
    int Npart=(int)(densinput*volume);
    float dens = Npart/volume;
    int Nstep = atof(argv[3]);
    int nsamp_ener = 10;
    int nsamp_pos = 100;
    int naccept = 0;
    float deltaR = 0.1;
    float Position[dim*Npart];

    // LOG of the run
    cout << "   Density = "<< dens <<endl;
    cout << "   Temperature = "<< Temp <<endl;

    // Open data files
    ofstream fich_ener, fich_posi;
    fich_ener.open("energy.txt");
    fich_posi.open("position.txt");

    // Random number generators
    uniform_real_distribution<> dis1(-1.0, 1.0);
    uniform_real_distribution<> dis2(0.0, 1.0);
    uniform_int_distribution<> dist(0, Npart-1);

    // Initial positions
    for (int i=0; i<3*Npart; ++i) Position[i] = dis1(generator)*Lbox/2.;

    // Monte carlo loop
    for (int istep = 0; istep<Nstep; ++istep) 
    {
        int itag = dist(generator);
        float PosNew[dim], PosOld[dim];
        float Enew, Eold, prob;

        // Move particle itag
        for (int k=0; k<dim; ++k) {
            PosOld[k] = Position[dim*itag + k];
            PosNew[k] = Position[dim*itag + k] + deltaR * dis1(generator);
        }

        //Energies of itag
        Eold = Energy(Position, PosOld, itag, Npart);
        Enew = Energy(Position, PosNew, itag, Npart);

        // Probabiliy ratio
        prob = exp(-(Enew-Eold)/Temp);

        float xi = dis2(generator);
        float Esample = Eold;

        // If transition is accept
        if (prob > xi) {
            for (int k=0; k<dim; ++k) Position[itag*dim + k] = PosNew[k];
            naccept = naccept + 1;
            Esample = Enew;
        }

        // Save sample energy and position
	// Here we will also calculate the radial distribution functions (rdf)
        if (istep % nsamp_ener == 0) fich_ener << 0.5*Esample << endl;       
        if (istep % nsamp_pos == 0) 
	{
            for (int i = 0; i < Npart; ++i) 
	    {
                for (int k = 0; k < dim; ++k)
		{
                    float pos_value = Position[i * dim + k];
                    fich_posi << pos_value << " ";
                }
                fich_posi << endl;
            }
        }
    }
    
    // Close data files
    fich_ener.close();
    fich_posi.close();

    // LOG of the run
    cout << "..............Done" <<endl;
    return 0;
}

// Lennard-jones energy function
float Energy(float Position[], float Pos_itag[], int itag, int Npart) {
    float U_tot = 0.0;
    float r_2cut = r_cut * r_cut;

    // Loop for every particle except itag
    for (int i = 0; i < Npart; ++i)  {
        if (i == itag) continue;
        float vec_dist[dim];
        float u_ij = 0.0;
        float r2_ij=0.0;
        
        // Obtain the distances r_ij
        for (int k = 0; k < dim; ++k) vec_dist[k] = (Pos_itag[k] - Position[i * dim + k]);        
        mic(vec_dist, Lbox);
        for (int k = 0; k < dim; ++k) r2_ij += vec_dist[k]*vec_dist[k]; 

        // Lennard-jones potential energy if we are under the cut
        if (r2_ij < r_2cut) {
            float r_mod = sqrt(r2_ij);
	        float r6 = pow((sigma/r_mod), 6.0);
	        float rc6 = pow((sigma/r_cut), 6.0);
	        u_ij = 4 * epsilon * ((r6*r6-r6) - (rc6*rc6-rc6));   
        }
        // Save the total energy
        U_tot += u_ij;
    }
    return U_tot;
}
