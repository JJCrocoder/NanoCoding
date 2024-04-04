#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <vector>
using namespace std;
const float pi = 3.14159265358979323846;

// Define functions
float g_interpolated(const vector<float>& x, const vector<float>& y, float x_new);
float func_ener(float r, const vector<float>& r_rdf, const vector<float>& u_rdf);
float func_press(float r, const vector<float>& r_rdf, const vector<float>& u_rdf);
float integral(float (*func)(float, const vector<float>&, const vector<float>&), float a, float b, int n, const vector<float>& r_rdf, const vector<float>& u_rdf);
float Energy(float position[], float N_O_Pos[], int itag, int Npart);
float distance_MIC(float Pos_i[], float Pos_j[], float Lbox, int dim);

// Global variables
const float sigma = 1.0;
const float epsilon = 1.0;
const float Lbox = 15.0;
const float volume = Lbox*Lbox*Lbox;
const int dim = 3;
const float r_cut = 2.5 * sigma;
const int nBins = 75;
const float dr = r_max / nBins;
const float deltaR = 0.1;
const double Nwid = 100;

// Sampling parameters
float density = 1;  // Modifica esto
float Temp = 1;  // Modifica esto
int Npart=(int)(density*volume);
const int Nstep = 1;  vModifica esto
const int nsamp_ener = 1;  // Modifica esto
const int nsamp_wid = 1;  // Modifica esto

// Main function
int main(int argc, char *argv[]) {
    random_device rand_dev;
    mt19937 generator(rand_dev());

    // Inicialization vectors
    float position[dim*Npart];
    float pressure_rdf = 0.0;
    float energy_rdf = 0.0;
    vector<int> g(nBins, 0);
    vector<float> u_rdf(nBins, 0);
    vector<float> r_rdf(nBins, 0);
    vector<float> energySamples;
    vector<float> chempotSamples;
  
    // Random number generators
    uniform_real_distribution<> dis11(-1.0, 1.0);
    uniform_real_distribution<> dis01(0.0, 1.0);
    uniform_int_distribution<> dis0N(0, Npart-1);

    // Initial positions
    for (int i=0; i<3*Npart; ++i) position[i] = dis11(generator)*Lbox/2.;

    // Monte carlo loop
    for (int istep = 0; istep<Nstep+1; ++istep) {
        float PosNew[dim], PosOld[dim];
        float E_new, E_old, prob, E_samp;

        // Generate itag
        int itag = dis0N(generator);

        // Move particle itag
        for (int k=0; k<dim; ++k) {
            PosOld[k] = position[dim*itag + k];
            PosNew[k] = position[dim*itag + k] + deltaR * dis11(generator);
        }

        //Energies of itag
        E_old = Energy(position, PosOld, itag, Npart);
        E_new = Energy(position, PosNew, itag, Npart);

        // Probabiliy ratio
        prob = exp(-(E_new-E_old)/Temp);
        float xi = dis01(generator);

        // If transition is accept
        if (prob > xi) {
            for (int k=0; k<dim; ++k) position[itag*dim + k] = PosNew[k];
            E_samp = E_new;
        }else{
            E_samp = E_old;
        }

        // Save energies
        if (istep % nsamp_ener == 0) energySamples.push_back(0.5 * E_samp);

        // Histogram for the RDF
        if (istep > (Nstep - 5)) { //Modificable cuanto coger de la RDF
            float Pos_i[dim];
            float Pos_j[dim];

            // Loop for i-j interactions
            for (int i = 0; i < Npart - 1; ++i) {
                for (int j = i + 1; j < Npart; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        // Calculate the positions between i-j
                        Pos_i[k] = position[i * dim + k];
                        Pos_j[k] = position[j * dim + k];
                    }

                    // Calculate the distance i-j
                    float r_mod = distance_MIC(Pos_i, Pos_j, Lbox, dim);

                    // If is in our range add to RDF
                    if (r_mod < r_max){
                        int index = static_cast<int>(r_mod / dr);
                        g[index] += 2.0; // We sum 2 beacause we use i->j as j->i
                    }
                }
            }
        }
    }
    fich_eq.close();

    // Normalization of the RDF and save it
    ofstream fich_rdf;
    fich_rdf.open("rdf.txt");
    for (int i = 0; i < nBins; ++i) {
        float r_in = dr * i;
        float r_out = dr * (i + 1);
        float r_mid = (r_in + r_out) / 2;

        // Calculate the ideal gas contribution
        float Nideal = (4.0 / 3.0) * pi * (r_out*r_out*r_out - r_in*r_in*r_in) * density * Npart * nsamp_rdf;

        // Normalize the rdf
        float rdf = g[i] / Nideal;

        // Save for interpolation
        u_rdf.push_back(rdf);
        r_rdf.push_back(r_mid);

        // Save the RDF
        fich_rdf << r_mid << " " << rdf << endl;
    }
    fich_rdf.close();
    
    // Save energy and pressure from RDF
    ofstream fich_pres_rdf("pressure_rdf.txt");
    pressure_rdf = -(2.0 / 3.0) * pi * density * density * integral(func_press, 0.1, r_cut, 1000, r_rdf, u_rdf);
    fich_pres_rdf << pressure_rdf << endl;
    fich_pres_rdf.close();
    
    ofstream fich_ener_rdf("energy_rdf.txt");
    energy_rdf = 2 * pi * density * integral(func_ener, 0.1, r_cut, 1000, r_rdf, u_rdf);
    fich_ener_rdf << energy_rdf << endl;
    fich_ener_rdf.close();

    // Widom insertion method and save it
    ofstream fich_wid("widom.txt");
    float Pos_ins[dim];
    for (int i = 0; i < nsamp_wid; i++) {
        float sum_Ewid = 0.0;

        //Insert particles
        for (int j = 0; j < Nwid; j++) {
            // Add a new particle
            for (int k = 0; k < dim; k++) Pos_ins[k] = dis11(generator)*Lbox/2.;
            // Obtain its energy
            float E_wid = Energy(position, Pos_ins, Npart + 1, Npart + 1);
            sum_Ewid += exp(-E_wid / Temp) / Nwid;
        }
        float chem_pot = -Temp * log(sum_Ewid);
        fich_wid << chem_pot << endl;
    }
    fich_wid.close();
    return 0;
}

// Linear interpolation function of RDF
float g_interpolated(const vector<float>& x, const vector<float>& y, float x_new) {
    // Find the closest points
    int i = 0;
    while (x[i] < x_new) i++;
    // Perform linear interpolation
    float x0 = x[i - 1];
    float x1 = x[i];
    float y0 = y[i - 1];
    float y1 = y[i];
    return y0 + (y1 - y0) * (x_new - x0) / (x1 - x0);
}

// Function for calculate the energy using the RDF
float func_ener(float r, const vector<float>& r_rdf, const vector<float>& u_rdf) {
    float r6 = pow((sigma/r), 6.0);
    float rc6 = pow((sigma/r_cut), 6.0);
    float u = 4 * epsilon * ((r6 * r6 - r6) - (rc6 * rc6 - rc6));
    return r * r * u * g_interpolated(r_rdf, u_rdf, r);
}

// Function for calculate the pressure using the RDF (quizás la ecuacion está mal)
float func_press(float r, const vector<float>& r_rdf, const vector<float>& u_rdf) {
    float r6 = pow((sigma/r), 6.0);
    float rc6 = pow((sigma/r_cut), 6.0);
    float u = 4 * epsilon * ((r6 * r6 - r6) - (rc6 * rc6 - rc6));
    return r * r * r * u * g_interpolated(r_rdf, u_rdf, r);
}

// Integral resolution (Simpson method)
float integral(float (*func)(float, const vector<float>&, const vector<float>&), float a, float b, int n, const vector<float>& r_rdf, const vector<float>& u_rdf) {
    float h = (b - a) / n;
    float sum = 0.0;

    // Initial and final points
    sum += func(a, r_rdf, u_rdf) + func(b, r_rdf, u_rdf);

    // Simpson method
    for (int i = 1; i < n; i += 2) sum += 4 * func(a + i * h, r_rdf, u_rdf);
    for (int i = 2; i < n - 1; i += 2) sum += 2 * func(a + i * h, r_rdf, u_rdf);
    return (h / 3) * sum;
}

// Lennard-Jones energy function
float Energy(float position[], float Pos_itag[], int itag, int Npart) {
    float U_tot = 0.0;
    float r2_cut = r_cut * r_cut;

    // Loop for every particle except itag
    for (int i = 0; i < Npart; ++i)  {
        if (i == itag) continue;
        float vec_dist[dim];
        float u_ij = 0.0;
        float r2_ij=0.0;
        
        // Obtain the distances r_ij between itag and all the other positions
        for (int k = 0; k < dim; ++k) vec_dist[k] = (Pos_itag[k] - position[i * dim + k]);
        for (int k = 0; k < dim; ++k) vec_dist[k] -= floorf(0.5 + vec_dist[k]/Lbox)*Lbox;      
        for (int k = 0; k < dim; ++k) r2_ij += pow(vec_dist[k], 2); 

        // Lennard-jones potential energy if we are under the cut
        if (r2_ij < r2_cut) {
            float r_mod = sqrt(r2_ij);
            float r6 = pow((sigma/r_mod), 6.0);
            float rc6 = pow((sigma/r_cut), 6.0);
            u_ij = 4 * epsilon * ((r6 * r6 - r6) - (rc6 * rc6 - rc6));
        }
        // Save the total energy
        U_tot += u_ij;
    }
    return U_tot;
}

// Distance module for Minimum image convention (periodic boundary conditions)
float distance_MIC(float Pos_i[], float Pos_j[], float Lbox, int dim) {
    float r_ij[dim];
    float r2_ij=0.0;

    // Calculate the distance for MIC
    for (int k = 0; k < dim; ++k) {
        r_ij[k] = Pos_j[k] - Pos_i[k];
        r_ij[k] -= floorf(0.5 + r_ij[k]/Lbox)*Lbox;
        r2_ij += pow(r_ij[k], 2);
    }  
    float mod_rij = sqrt(r2_ij);
    return mod_rij;
}
