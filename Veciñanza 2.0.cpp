#include <iostream>
#include<fstream>
#include<random>
#include<cmath>
#include<stdlib.h>
#include<algorithm>
#include<vector>

using namespace std;

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





double Energy(vector <double> Position, double N_O_Pos[], int itag, double r_cut, int N_part, int cell[]);
double LJ_pot(double r, double r_cut);
//function declaration
void getcell(double normpos[], int cell[]);
void make_linked_list(vector <double> Position, vector<int>& list, vector<int>& head, int N_part); //Fills head and list with the current system state
int pbc_cells(int icell);
int cell_index(int cell[]);
void mic(double temppos[], double L_box);
void list_reconstruction(double New_Pos[], vector<int>& list, vector<int>& head, int itag);






const int N_part = 455;      // Número de Partículas
const double temp = 2.5;         // Temperatura del sistema
const double L_box = 10;       // Longitud de los lados de la caja
const double r_cut = 2.5;       // Cutof del radio
const int N_dim = 3;          // Número de dimensiones consideradas

//Constantes del potencial de Leonardo-Cojones
const double epsilon = 1.0;
const double sigma = 1.0;

// We calculate the number of cells in each dimension

double cellsize = 1.2 * r_cut;		// A bit bigger than r_cut
int ncells = static_cast<int>(L_box / cellsize);	// we obtain the number of cells
double mcells = ncells * ncells * ncells;	// Now the total number of cell

// Variables del Montecarlo (el algoritmo no el casino, ojalá fuera el casino)
const int N_step = 1000000;     // Número de saltos en el Monte Carlo

int Ener_cut = 500000;

int nsamp_ener = 1000;    // Ratio de toma de medidas de la energía
int nsamp_pos = 10;    // Ratio de toma de medidas de la posición
int naccept = 0;        // Número de saltos aceptados (masonada)
double deltaR = 0.1;     // Tamaño del salto

vector <int>list;
vector<int>head;
vector<double> Energies;

int main() {
    for (int i = 0; i < N_part; i++) list.push_back(-1);
    for (int j = 0; j < mcells; j++) head.push_back(-1);
    
    double dens = N_part / pow(L_box, 3);
    cout << dens << endl;
    cout << temp << endl;
    vector <double> Position; // Array donde vienen las posiciones de las partículas (array tipo GOD)



    // Definimos la faquin configuración inicial
    for (int k = 0; k < N_dim * N_part; ++k) {
        Position.push_back(uniform(-0.5 * L_box, 0.5 * L_box)); // We create a random coordinate
    }
    

    make_linked_list(Position, list, head, N_part);

    
    // Bucle principal del Monte Carlo (33 negro impar y pasa)


    
    for (int j = 0; j < mcells; j++) cout << head[j] << endl;

    


    
    for (int istep = 0; istep < N_step; ++istep) {
        int itag = rand_num(0, N_part - 1);  //random integer
        double PosNew[N_dim], PosOld[N_dim], ideal_Pos[N_dim];
        double Enew, Eold;



        for (int k = 0; k < N_dim; ++k)
        {
            PosOld[k] = Position[N_dim * itag + k];
            PosNew[k] = Position[N_dim * itag + k] + uniform(-deltaR, deltaR);
            if (PosNew[k] > L_box / 2) PosNew[k] -= L_box / 2;
            else if (PosNew[k] < -L_box / 2) PosNew[k] += L_box / 2;


        }

        int ocell[N_dim];
        int ncell[N_dim];
        

        getcell(PosOld,ocell);
        getcell(PosNew,ncell);


        Eold = Energy(Position, PosOld, itag, r_cut, N_part, ocell);
        Enew = Energy(Position, PosNew, itag, r_cut, N_part, ncell);



        
        // Probability ratio

        double prob = exp(-(Enew - Eold) / temp);
        double xi = uniform(0, 1);

        double Esample = Eold;

        if (prob > xi) {
            for (int k = 0; k < N_dim; ++k) {
                Position[itag * N_dim + k] = PosNew[k];
            }
            naccept++;
            Esample = Enew;
        }

        



        if (istep % nsamp_ener == 0 and istep > Ener_cut) {
            if (isnan(Esample)) continue;
            else {
                // Abrir el archivo y escribir el valor solo si no es NaN

                Energies.push_back(0.5 * Esample);
                
            }

        }

        list_reconstruction(PosNew,  list, head,itag);



    }
    //Media
    double suma = 0.0;
    for (int i = 0; i < Energies.size(); i++) {
        suma += Energies[i] / Energies.size();
    }

    cout << "La energía media es: "<<suma << endl;

    
    

    return 0;

}



double Energy(vector <double> Position, double N_O_Pos[], int itag, double r_cut, int N_part, int cell[]) {
    double U_itag = 0.0;  // The energy of the selected particle
    // Loop for the particles

    for (int jx = cell[0] - 1; jx <= cell[0] + 1; jx++) {
        for (int jy = cell[1] - 1; jy <= cell[1] + 1; jy++) {
            for (int jz = cell[2] - 1; jz <= cell[2] + 1; jz++) {

                int prueba = pbc_cells(jx) + pbc_cells(jy) * ncells + pbc_cells(jz) * ncells * ncells;
                
                if (prueba < mcells && prueba>0) {
                    int j = head[pbc_cells(jx) + pbc_cells(jy) * ncells + pbc_cells(jz) * ncells * ncells];
                    
                    while (j != -1) {
                        if (j == itag) j = list[j];
                        else {

                            //Now we construct the distance between 2 atoms
                            double rj[N_dim];  // Position array of 2nd atom
                            double rij[N_dim];  // Relative position vector
                            double mod_rij; // modulus of the vector

                            // We perform the substraction and calculate the modulus
                            for (int k = 0; k < N_dim; ++k) {
                                rj[k] = Position[j * 3 + k];
                                rij[k] = N_O_Pos[k] - rj[k];
                                rij[k] -= L_box * floor(rij[k] / L_box + 0.5);
                            }
                            mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
                            U_itag += LJ_pot(mod_rij, r_cut);

                            j = list[j];
                        }

                    }
                }
            }
        }
    }

    return U_itag;

}

double LJ_pot(double r, double r_cut) {
    if (r < r_cut * epsilon) {
        double term1 = pow(epsilon / r, 12) - pow(epsilon / r, 6);
        double termcut = pow(epsilon / r_cut, 12) - pow(epsilon / r_cut, 6);
        return 4 * epsilon * (term1 );
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

void mic(double temppos[], double L_box) {
    for (int k = 0; k < 3; ++k) {
        
        temppos[k] -= L_box * floor(temppos[k] / L_box + 0.5);
        
    }

}

void make_linked_list(vector <double> Position, vector<int>& list, vector<int>& head, int N_part) {

    for (int i = 0; i < N_part; ++i) {   // We go trough all the particles
        double temppos[N_dim];  // For each particle, we initialize an auxiliar position vector 
        for (int k = 0; k < N_dim; ++k) temppos[k] = Position[N_dim * i + k]; // We store the coordinates in the auxiñiar vector 
        mic(temppos, L_box);// Fold the position of particle i using mic
        //for (int i = 0; i < N_dim; i++) temppos[i] += L_box * 0.5;
        int cell[N_dim];
        getcell(temppos,cell);
        int icell = cell_index(cell); // Find the cell where temppos belongs to
        
        // build the list and head arrays
        list[i] = head[icell];   // The list element for the current particle points to the last saved particle of each cell
        head[icell] = i;         // Now the last saved particle will be the current one
        
    }
}


void getcell(double normpos[], int cell[]) {
    double pos[N_dim];
    for (int i = 0; i < N_dim; i++) pos[i] = normpos[i] / L_box;
    // Calcular las celdas en cada dirección
    for (int k = 0; k < N_dim; ++k) {
        cell[k] = static_cast<int>((0.5 + pos[k]) * ncells);
    }
}

int cell_index(int cell[]) {
    int index = cell[0] + cell[1] * ncells + cell[2] * ncells * ncells;
    return index;
}



void list_reconstruction(double New_Pos[], vector<int>& list, vector<int>& head, int itag) {
    
    for (int i=0;i<N_dim;i++) New_Pos[i]-=L_box * floor(New_Pos[i] / L_box + 0.5);
    

    for (int i = 0; i < N_part; i++) {
        if (list[i] == itag) list[i] = list[itag];
    }
    
    for (int n = 0; n < mcells; n++) {
        if (head[n] == itag) head[n] = list[itag];
    }



    
    
    
    int new_cell[N_dim];
    getcell(New_Pos, new_cell);
    
    int ncell = cell_index(new_cell);    
    
    
    list[itag] = head[ncell];
    head[ncell] = itag;


}


// Give me the cell index in head taking into account PBC
int pbc_cells(int icell) {
    if (icell == -1) return ncells - 1;
    else if (icell == ncells) return 0;
    else return icell;
}
