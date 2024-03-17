#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

// Prototypes of the functions
void mic(float* vec, float Lbox);
void moveParticles(float* Position);
void rdf(float* Pos);

// Global variables
float sigma = 1.0;
float epsilon = 1.0;
float Lbox = 10.0;
float volume = Lbox * Lbox * Lbox;
int dim = 3;
float Temp = 1.9;
float r_cut = 2.5 * sigma;
float DIFF = 1.0;
int Nstep = 5000;
int nrdf = 10;
int ndraw = 100;
float deltat = 0.001;
float densinput = 0.5;
const int Npart = static_cast<int>(densinput * volume);
float dens = Npart / volume;

int main() {
    // Variables


    // LOG of the run
    cout << "Density " << dens << endl;
    cout << "Temperature " << Temp << endl;

    // System variables
    float Pos[dim * Npart];

    // Data files
    ifstream fich_pos;
    fich_pos.open("position.txt");
    // Initial positions: prepared from monte carlo: read pos
    for (int i = 0; i < 3 * Npart; ++i)
        fich_pos >> Pos[i];
    fich_pos.close();

    // START BD LOOP
    for (int istep = 0; istep < Nstep; ++istep) {
        // Calculate forces : For
      moveParticles(Pos);
      
      if (istep % ndraw == 0) cout << Pos[1] << endl;

      if (istep % nrdf == 0) {
          rdf(Pos);
      }
    }

    return 0;
}

// Minimum image convention
void mic(float* vec, float Lbox) {
    for (int i = 0; i < 3; ++i) {
        vec[i] -= floorf(0.5 + vec[i] / Lbox) * Lbox;
    }
}

// Move particles using Brownian dynamics
void moveParticles(float* Position) {
    // Calculate forces
    float force[dim * Npart];
    float r_cut = 2 * sigma;

    for (int i = 0; i < Npart * dim; ++i)
        force[i] = 0;

    for (int i = 0; i < Npart; ++i){
      float Pos_itag[dim]={0};
      for (int k =0; k<dim;++k){
        Pos_itag[k] = Position[i*dim + k];
      }

        for (int j = i + 1; j < Npart; ++j) {
            float r_2cut = r_cut * r_cut;
            float r2_Ij = 0.0;
            float vec_dist[dim];
            float u_Ij = 0.0;

            for (int k = 0; k < dim; ++k) {
                vec_dist[k] = Pos_itag[k] - Position[j * dim + k];
                mic(vec_dist, Lbox);
                r2_Ij += pow(vec_dist[k], 2);
            }

            if (r2_Ij < r_2cut) {
                float r_mod = sqrt(r2_Ij);
                float r6 = pow((sigma / r_mod), 6.0);
                float rc6 = pow((sigma / r_cut), 6.0);
                float f_ij = -48 * epsilon * r6 * (r6 - 0.5) / r_mod;

                for (int k = 0; k < dim; ++k) {
                    force[i * dim + k] -= (f_ij * vec_dist[k] / r_mod);
                    force[j * dim + k] += (f_ij * vec_dist[k] / r_mod);
                }
            }
        }
    }

    // Move particles using Euler method
    for (int i = 0; i < Npart; ++i) {
        for (int k = 0; k < dim; ++k) {
            Position[i * dim + k] += 0.5 * deltat * force[i * dim + k];
            // Apply periodic boundary conditions
            Position[i * dim + k] -= floorf(0.5 + Position[i * dim + k] / Lbox) * Lbox;
        }
    }
}

// RDF calculation
void rdf(float* Pos) {
    std::ofstream fich_rdf;
    fich_rdf.open("RDF.txt", std::ofstream::app); // Open file in append mode

    const float pi = 3.14159265358979323846;
    const int num_bins = 120; // Number of bins
    float min_val = 0.0;
    float max_val = Lbox / 2;
    float rango = max_val - min_val;
    float DeltaX = rango / num_bins;
    std::vector<int> Histo(num_bins, 0);

    float Const = 4.0 / 3.0 * dens * pi;

    for (int p = 0; p < Npart - 1; ++p) {
        for (int q = p + 1; q < Npart; ++q) {
            float Pos_i[3];
            float Pos_j[3];

            for (int k = 0; k < 3; ++k) {
                Pos_i[k] = Pos[q * 3 + k];
                Pos_j[k] = Pos[p * 3 + k];
            }

            float distancia = 0.0;
            for (int k = 0; k < 3; ++k) {
                float dist = Pos_i[k] - Pos_j[k];
                mic(&dist, Lbox);
                distancia += dist * dist;
            }

            distancia = sqrt(distancia);

            if (distancia < max_val) {
                int IBIN = static_cast<int>((distancia - min_val) / DeltaX);
                Histo[IBIN] += 2;
            }
        }
    }

    for (int k = 0; k < num_bins; ++k) {
        float rlow = min_val + DeltaX * k;
        float rup = min_val + DeltaX * (k + 1);
        float rrr = (rup + rlow) / 2;
        float Nideal = Const * (pow(rup, 3) - pow(rlow, 3));
        float GR = Histo[k] / (Npart * Nideal);
        fich_rdf << rrr << ' ' << GR << std::endl;
    }

    fich_rdf.close();
}
