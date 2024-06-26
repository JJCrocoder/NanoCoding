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

// Global variables or constants
#include "parameters.h"
const double pi = 3.14159265358979323846;
float volume = pow(Lbox,dim);

// the ".h" libraries are refered to files that have been created to be modified outside
#include "../neighbors.h"
#include "../mic_pbc.h"
#include "../my_random.h"
#include "../my_potentials.h"

// Various function definitions
float Int_Energy(float * Position, float * N_O_Pos, int itag, int Npart);
float calculate_virial(float * Position, int Npart);

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

    float densinput = atof(argv[1]); 	// We store the input density
    float Temp = atof(argv[2]);		// Same for temperature
    int Npart=(int)(densinput*volume);	// We calculate the number of particles with the density
    float dens = Npart/volume;		// And fix the real density value	
    // Variables
    int Nstep = 200*Npart;
    int nsamp_ener = 4*Npart;		// Steps for energy sampling
    int nsamp_pos = Npart;		// Steps for position sampling
    int naccept = 0;			// counter for accepted particle displacements
    int ncount = 0;			// 
    float deltaR = 0.1;			// Size of the metropolis random displacement
    float Position[dim*Npart];		// Particles positions array
    vector<int> Histo(num_bins, 0);  //Vector which elements are related with each one of the different bins
    float Const = 4.0 / 3.0 * dens * pi; //Constant needed for normalization
    float Int_Energy1 = 0.0;
    
    // LOG of the run
    cout << " Density = "<< dens <<endl;
    cout << " Temperature = "<< Temp <<endl;
    cout << " Number of particles  = "<< Npart <<endl;

    // Open data files
    ofstream fich_ener, fich_posi, fich_rdf;
    fich_ener.open("energy.txt");
    fich_posi.open("position.txt");
    fich_rdf.open("rdf_hist.txt");

    // Initial positions
    for (int i=0; i<dim*Npart; ++i) Position[i] = uniform(-0.5*Lbox, 0.5*Lbox);
    cout << " Initial configuration created"<<endl;
    cout << endl << "Begining of Montecarlo loop:"<<endl;
    cout << "\t" << "Energy per particle"<<endl;
    
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
        Eold = Int_Energy(Position, PosOld, itag, Npart);	// Current energy
        Enew = Int_Energy(Position, PosNew, itag, Npart);	// Modified Int_Energy

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
	Int_Energy1+= 0.5*Esample;
	
        // Save sample energy and position
	// Here we will also calculate the radial distribution function (rdf)
        if (istep % nsamp_ener == 0) {
	  fich_ener << Int_Energy1/nsamp_ener*Npart << endl; // Storing energy
	  cout << "\t" << Int_Energy1/nsamp_ener <<endl;
	  Int_Energy1 = 0;
	}      
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
    // cout<<"Distance\trdf"<<endl;
    for (int k = 0; k < num_bins; ++k) 
    {
        float rrr = dmin + DeltaX * (k + 0.5);
        float Nideal = 4.0 / 3.0 * pi * dens * (pow(rrr + 0.5*DeltaX, 3) - pow(rrr - 0.5*DeltaX, 3));
        float GR = Histo[k] / (Npart * Nideal * ncount);
        // cout << rrr << '\t' << GR << endl;
    }

    // Close data files
    fich_rdf.close();
    fich_ener.close();
    fich_posi.close();

    // LOG of the run
    cout << "..............Done" <<endl;
    return 0;
}

// Interaction energy function
float Int_Energy(float * Position, float * Pos_itag, int itag, int Npart) {
    
    float U_tot = 0.0;

    // Loop for every particle except itag
    for (int i = 0; i < Npart; ++i)  {
        if (i == itag) continue;
        float vec_rij[dim];
        float r2_ij=0.0;
        
        // Obtain the distance r_ij
        for (int k = 0; k < dim; ++k) vec_rij[k] = (Pos_itag[k] - Position[i * dim + k]); // Coordinates       
        mic(vec_rij, Lbox);	// Applying MIC
        for (int k = 0; k < dim; ++k) r2_ij += vec_rij[k]*vec_rij[k]; // Squared modulus

        // pair interaction evaluation and addition
        U_tot += lj_pot(r2_ij);	// Defined in an external file
    }
    return U_tot;
}

// Function for virial calculations
float calculate_virial(float * Position, int Npart) {
    float virial = 0.0;
    for (int i = 0; i < Npart - 1; ++i) {
        for (int j = i + 1; j < Npart; ++j) {
		float rij[dim] = {0.0};
	    	for(int k = 0; k<dim; ++k){
			rij[k] = Position[j * dim + k] -  Position[j * dim + k];
			// Quiero obtener la fuerza sin tener que crear arrays i, j
			// Force[k] = 
			// virial += rij[k]*Force[k];
		}
        }
    }
    return virial/(dim*volume);
}


