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

// the ".h" libraries are refered to files that have been created to be modified outside from this program
#include "parameters.h"
#include "../neighbors.h"
#include "../mic_pbc"
#include "../my_random"

// Various function definitions
float Energy(float Position[], float N_O_Pos[], int itag, int Npart);

// Global variables or constants
const double pi = 3.14159265358979323846;
float volume = pow(Lbox,dim);
int Npart=(int)(densinput*volume);	// We calculate the number of particles with the density
float dens = Npart/volume;		// And fix the real density value

// rdf variables
float dmax = 0.5 * Lbox;
float dmin = 0.0;
int num_bins = 100;
float DeltaX = (dmax - dmin)/num_bins;

// MAIN FUNCTION

// The compiled executable can get some data as input
// "argc" gives the number or arguments that the function has accepted as input, including the exceuting command
// "argv" is an array (string type) that includes each of the arguments given as inputs
// For example, if we run ./main 5 3.0 as an executable: argc = 3, argv = {"./main", "5", "3.0"}
// For this code, the arguments should be the system density, temperature and number of simulation steps.
int main(int argc, char *argv[]) {

    // Variables
    int nsamp_ener = 10;		// Steps for energy sampling
    int nsamp_pos = 100;		// Steps for position sampling
    int naccept = 0;			// counter for accepted particle displacements
    int ncount = 0;			// 
    float deltaR = 0.1;			// Size of the metropolis random displacement
    float Position[dim*Npart];		// Particles positions array
    vector<int> Histo(num_bins, 0);  //Vector which elements are related with each one of the different bins
    float Const = 4.0 / 3.0 * dens * pi; //Constant needed for normalization
	
    // LOG of the run
    cout << " Density = "<< dens <<endl;
    cout << " Temperature = "<< Temp <<endl;

    // Open data files
    ofstream fich_ener, fich_posi, fich_rdf;
    fich_ener.open("energy.txt");
    fich_posi.open("position.txt");
    fich_rdf.open("rdf_hist.txt");

    // Initial positions
    for (int i=0; i<dim*Npart; ++i) Position[i] = uniform(-0.5*Lbox, 0.5*Lbox);

    // Montecarlo loop
    for (int istep = 0; istep<Nstep; ++istep) 
    {
        int itag = rand_num(0,Npart-1);	// A random index for selecting a particle
        float PosNew[dim], PosOld[dim];	// Position arrays for the selected particle
        float Enew, Eold, prob;		// Energies and probabilities variables

        // Move selected particle
        for (int k=0; k<dim; ++k) {
            PosOld[k] = Position[dim*itag + k]; // The current position
            PosNew[k] = Position[dim*itag + k] +  uniform(-deltaR, deltaR); // The modified position
        }
	
        //Energies of the selected particle
        Eold = Energy(Position, PosOld, itag, Npart);	// Current energy
        Enew = Energy(Position, PosNew, itag, Npart);	// Modified Energy

        prob = exp(-(Enew-Eold)/Temp);	// Probabiliy ratio

        float xi = uniform(0.0, 1.0);	// Auxiliar variable for accept-reject
        float Esample = Eold;		// We initialize the sampling with the current energy

        // For acceptance
        if (prob > xi) {
	    // The system configuration is modified
            for (int k=0; k<dim; ++k) Position[itag*dim + k] = PosNew[k];
            naccept = naccept + 1;	// We add an accept
            Esample = Enew;		// The sampled energy is modified
        }
  
        // Save sample energy and position
	// Here we will also calculate the radial distribution function (rdf)
        if (istep % nsamp_ener == 0) fich_ener << 0.5*Esample << endl; // Storing energy      
        if (istep % nsamp_pos == 0) 
	{
            for (int i = 0; i < Npart - 1; ++i)	// For each particle
	    {
		for (int k = 0; k < dim; ++k)		// For each component 
		{
		    fich_posi << Position[i * dim + k] << " "; // we store each component in a single line
		}
	        fich_posi << endl;	// Every particle is a diferent line

		// RDF
		    
		if (i == Npart) continue; // We don't count the last particle in "i" index for RDF calculation
		    
		//We calculate the distances between every pair of particles   
                for (int j = i + 1; j < Npart; ++j) {
                    float Pos_i[dim];
                    float Pos_j[dim];
		    float rij[dim];

                    for (int k = 0; k < dim; ++k) {
                        Pos_i[k] = Position[i * dim + k];
                        Pos_j[k] = Position[j * dim + k];
			rij[k] = Pos_j[k] - Pos_i[k];
                    }
		   
                    //Take into consideration the MIC (Minimum Image Convention)
                    float dist = mic_distance(Pos_i, Pos_j, Lbox);
                    //Analyze if such distance is inside the range of our histogram
                    if (dist > dmax) continue;
		    //Calculate to which bin such distance belongs
		    int IBIN = floor((dist - dmin) / DeltaX);
		    //Remember that the distance is of a pair of particles, so it contributes twice in our histogram.
		    Histo[IBIN] += 2;
		    ncount += 2;
                }
            }
	fich_posi<<endl<<"#"<<endl; // Add time frame separators (minipunto)
        }
    }

    // Normalization of the RDF histogram
    cout<<"Distance\trdf"<<endl;
    for (int k = 0; k < num_bins; ++k) 
    {
        float rrr = dmin + DeltaX * (k + 0.5);
        float Nideal = 4.0 / 3.0 * pi * dens * (pow(rrr + 0.5*DeltaX, 3) - pow(rrr - 0.5*DeltaX, 3));
        float GR = Histo[k] / (Npart * Nideal * ncount);
        cout << rrr << '\t' << GR << endl;
    }

    // Close data files
    fich_rdf.close();
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
