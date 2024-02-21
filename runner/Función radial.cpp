#include <iostream>
#include<fstream>
#include<random>
#include<cmath>
#include<stdlib.h>
#include<algorithm>
#include<vector>

using namespace std;
const double pi = 3.14159265358979323846;

//We declare the random generators
double uniform(double min, double max) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(min, max);

    return dis(gen);
}

int rand_num(int min, int max) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(min, max);

    return dis(gen);
}




//We declare the functions that we are going to use
double Energy(double Position[], double N_O_Pos[], int itag, double r_cut);
double LJ_pot(double r, double r_cut);
double MIC_distancia(double Position_part_i[], double Postion_part_j[]);


const double temp = 2.5;         // Temperatura of the system
const double L_box = 10;       // Size of the sides of the box
const double r_cut = 2.5;       // Radius cut-off

const int N_part = 1000; //We define the number of particles.

//In your codes it may be easy (or even better) to define the density and volume and from them the number of particles. I am going to do it this way since my C++ editor complains after if not, when building the position vector.

const int N_dim = 3; // Number of considered dimensions (3D)
double dens = N_part / pow(L_box, 3); //Density of particles of the system


//Constants for building the potential of Leonardo-Cojones
const double epsilon = 1.0;
const double sigma = 1.0;


// Variables for the Monte Carlo simulation
const int N_step = 1000000;     // Number of steps in the simulation
int nsamp_pos = 10000;    // Ratio at which we sample the positions
int naccept = 0;        // Number of accepted jumps. May be useful in future codes, but not in this case.
double deltaR = 0.2;     // Size of the jump
int tsteps = N_step / nsamp_pos;
// You can introduce here like a type of cut for only calculating the positions of the steps above this cut, ensuring with that we are sampling the system in equilibrium

//Variables of our histogram
const int num_bins = 100; //Number of bins
double min_val = 0.0;  //Lower value of r we take into account
double max_val = L_box/2;  //Higher value of r we take into account
double rango = max_val - min_val; //Range over which we are going to build the function (histogram)
double DeltaX = rango / num_bins;  //Width of the bins.
vector<int> Histo(num_bins, 0);  //Vector which elements are related with each one of the different bins
double Const = 4.0 / 3.0 * dens * pi; //Constant needed for normalization


//Main code
int main() {
    
    const int len_array = N_part * N_dim;  //I define this variable since my editor was giving me problems when defining the position vector. It may not be needed in your cases.
    
    
    cout <<"Densidad: " << dens << endl; //Density of particles
    cout <<"Temperatura: "<< temp << endl; //Temperature of the system
    double Position[len_array];  //Array where we save and update the positions of the particles.

    //We open the text file we are going to use
    ofstream fich_rdf;

   
    
    fich_rdf.open("Radial distribution function.txt"); //File for saving the histogram
    
    // We define the initial configuration of our particles
    for (int k = 0; k < N_dim * N_part; ++k) {
        Position[k] = uniform(-0.5 * L_box, 0.5 * L_box); // We create a random coordinate 
    }




    // Main loop of the Monte Carlo simulation

    for (int istep = 0; istep < N_step; ++istep) {
        int itag = rand_num(0, N_part);  //random integer for chosing the moving particle
        double PosNew[N_dim], PosOld[N_dim];
        double Enew, Eold;

        for (int k = 0; k < N_dim; ++k)
        {
            PosOld[k] = Position[N_dim * itag + k]; //Old position of the particle
            PosNew[k] = Position[N_dim * itag + k] + uniform(-deltaR, deltaR); //We calculate the new position of the particle
            //We take into consideration if the particle leaves the box (Periodic Boundary Conditions)
            if (PosNew[k] > L_box / 2) PosNew[k] -= L_box; 
            else if (PosNew[k] < -L_box) PosNew[k] += L_box;
            else continue;
        }
        //Old energy of the system
        Eold = Energy(Position, PosOld, itag, r_cut);
        //New energy of the system after the particle jump
        Enew = Energy(Position, PosNew, itag, r_cut);




        // Probability ratio
        double prob = exp(-(Enew - Eold) / temp);
        double xi = uniform(0, 1);

        double Esample = Eold;
        
        //Acceptance or rejection of the particle jump
        if (prob > xi) {
            for (int k = 0; k < N_dim; ++k) {
                Position[itag * N_dim + k] = PosNew[k];
            }
            naccept++; //We calculate the number of accepted jumps. As said, not useful for this problem.
            Esample = Enew;
        }

        //Position sampling
        if (istep % nsamp_pos == 0) {
            //We calculate the distances between every pair of particles
            for (int i = 0; i < N_part - 1; ++i) {
                for (int j = i + 1; j < N_part; ++j) {
                    double Pos_i[3];
                    double Pos_j[3];

                    for (int k = 0; k < 3; ++k) {
                        Pos_i[k] = Position[i * N_dim + k];
                        Pos_j[k] = Position[j * N_dim + k];
                    }
                    //Take into consideration the MIC (Minimum Image Convention) and calculate the distance between particle i and j
                    double distancia = MIC_distancia(Pos_i, Pos_j);
                    //Analyze if such distance is inside the range of our histogram
                    if (distancia < max_val) {
                        //Calculate to which bin such distance belongs
                        int IBIN = static_cast<int>((distancia - min_val) / DeltaX);
                        //Remember that the distance is of a pair of particles, so it contributes twice in our histogram.
                        Histo[IBIN] += 2;
                    }
                }
            }
        }
    }
        
        // Normalization of the histogram
        for (int k = 0; k < num_bins; ++k) {
            double rlow = min_val + DeltaX * k;
            double rup = min_val + DeltaX * (k + 1);
            double rrr = (rup + rlow) / 2;
            double Nideal = Const * (pow(rup, 3) - pow(rlow, 3));
            double GR = Histo[k] / (N_part * Nideal * tsteps);
            fich_rdf << rrr << ' ' << GR << endl;
        }



    //Cerramos los ficheros que entra el frío
    fich_rdf.close();
   

    return 0;

}



double Energy(double Position[], double N_O_Pos[], int itag, double r_cut)
{
    double U_itag = 0;  // The energy of the selected particle
    // Loop for the particles
    for (int i_sum = 0; i_sum < N_part; ++i_sum) {
        if (i_sum != itag) {

            //Now we construct the distance between 2 atoms
            double rj[N_dim];  // Position array of 2nd atom
            double rij[N_dim];  // Relative position vector
            double mod_rij; // modulus of the vector

            // We perform the substraction and calculate the modulus
            for (int k = 0; k < N_dim; ++k) {
                rj[k] = Position[i_sum * 3 + k];
                rij[k] = N_O_Pos[k] - rj[k];
                rij[k] -= L_box * floor(rij[k] / L_box + 0.5);
            }
            mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
            U_itag += LJ_pot(mod_rij, r_cut);

        }
        else continue;
    }

    return 0.5 * U_itag;
}

double LJ_pot(double r, double r_cut) {
    if (r < r_cut * epsilon) {
        double term1 = pow(epsilon / r, 12);
        double termcut = pow(epsilon / r, 6);
        return 4 * epsilon * (term1 - termcut);
    }
    else {
        return 0.0;
    }
}


double MIC_distancia(double Position_part_i[], double Postion_part_j[]) {
    double rij[3];
    for (int k = 0; k < 3; ++k) {
        rij[k] = Postion_part_j[k] - Position_part_i[k];
        rij[k] -= L_box * floor(rij[k] / L_box + 0.5);    
    }
    double mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
    return mod_rij;
    }


