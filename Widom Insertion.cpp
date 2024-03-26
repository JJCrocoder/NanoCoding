#include <iostream>
#include<fstream>
#include<random>
#include<cmath>
#include<stdlib.h>
#include<algorithm>

using namespace std;
const double pi = 3.14159265358979323846;

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





double Energy(double Position[], double N_O_Pos[], int itag, double r_cut); 
double LJ_pot(double r, double r_cut);  
double MIC_distancia(double Position_part_i[], double Position_part_j[]); 
double media(double vector[], int medidas);
double varianza(double vector[], int medidas);




const int N_part = 1235;      // Número de Partículas
const double temp = 3.1;         // Temperatura del sistema
const double L_box = 15;       // Longitud de los lados de la caja
const double r_cut = 2.5;       // Cutof del radio
const int N_dim = 3;          // Número de dimensiones consideradas
const double dens = N_part / pow(L_box, 3); //Densidad de partículas del sistema

//Constantes del potencial de Leonardo-Cojones
const double epsilon = 1.0;
const double sigma = 1.0;

// Variables del Montecarlo (el algoritmo no el casino, ojalá fuera el casino)
const int N_step = 1000000;     // Número de saltos en el Monte Carlo

int naccept = 0;        // Número de saltos aceptados (masonada)
double deltaR = 0.2;     // Tamaño del salto

//WIDOM INSERTION
double num_ins = 10000; //Número de inserciones que vamos a hacer para medir el potencial químico
const int num_med = 100; //Número de medidas que hacemos del potencial químico para promediar



int main() {

    cout << "La densidad de partículas del sistema es: " << dens << endl;
    cout << "La temperatura del sistema es: " << temp << endl;

    double Position[N_dim * N_part]; // Array donde vienen las posiciones de las partículas (array tipo GOD)
    double Pot_Quim[num_med];  //Array donde introduciremos las distintas medidas obtenidas del potencial químico
    


    // Definimos la faquin configuración inicial
    for (int k = 0; k < N_dim * N_part; ++k) {
        Position[k] = uniform(-0.5 * L_box, 0.5 * L_box); // We create a random coordinate
    }




    // Bucle principal del Monte Carlo (33 negro impar y pasa)

    for (int istep = 0; istep < N_step; ++istep) {
        int itag = rand_num(0, N_part);  //random integer
        double PosNew[N_dim], PosOld[N_dim];
        double Enew, Eold;

        for (int k = 0; k < N_dim; ++k)
        {
            PosOld[k] = Position[N_dim * itag + k];
            PosNew[k] = Position[N_dim * itag + k] + uniform(-deltaR, deltaR);


        }
        Eold = Energy(Position, PosOld, itag, r_cut);
        Enew = Energy(Position, PosNew, itag, r_cut);




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

    }
    cout << "Monte Carlo corrido" << endl;

    //WIDOM INSERTION

    double Pos_New_particle[N_dim]; //Posición de la partículas que insertamos
    int index = 0;   //Índice para conocer que medida estamos realizando

    //Inserción de la partícula
    for (int j = 0; j < num_med; j++) {
        double sum_ins = 0.0;  //En esta variable se va realizando las operaciones necesarias para la computación del sumatorio que hay definiendo el potencial químico
        for (int i = 0; i < num_ins; i++) {
            for (int k = 0; k < N_dim; k++) Pos_New_particle[k] = 0.0; //Devolvemos a cero la posición de la partícula insertada para evitar problemas
            for (int k = 0; k < N_dim; k++) Pos_New_particle[k] = uniform(-0.5 * L_box, 0.5 * L_box); //Posición aleatorio en la que se inserta la partícula
            double New_Ener = Energy(Position, Pos_New_particle, N_part + 1, r_cut); //Energía de la inserción
            sum_ins += exp(-New_Ener / temp) / num_ins;  //Sumamos a la variable que computa el sumatorio
        }
        double mu = -temp * log(sum_ins); //Con el sumatorio hecho computamos el potencial químico de esta medida
        Pot_Quim[index] = mu; //Lo guardamos en el array en el que recogemos los distintos valores obtenidos de este potencial
        index += 1; //Índice de la siguiente medida
    }

    double mu_pot = media(Pot_Quim, num_med); //Computamos la media de todos los valores obtenidos del potencial químico
    cout << "El potencial químico del sistema es: " << mu_pot << endl;
    double var_mu = varianza(Pot_Quim, num_med); //Y su varianza
    cout << "La desviación estandar del potencil quimico es: " << var_mu << endl;

    //Nótese que si se tiene en cuenta la parte truncada del potencial va a dar un poco más de lo que pone en el Moodle.
    //Si se toma el potencial L-J normal, sin truncar, da lo mismo exacto que el Moodle.

    




       

  
   


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

    return U_itag;

}

double LJ_pot(double r, double r_cut) {
    if (r < r_cut * epsilon) {
        double term1 = pow(epsilon / r, 12) - pow(epsilon / r, 6);
        double termcut = pow(epsilon / r_cut, 12) - pow(epsilon / r_cut, 6);
        return 4 * epsilon * (term1-termcut);
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


double media(double vector[], int medidas) {
    double mean;
    double sum = 0;
    for (int i = 0; i < medidas; i++) sum += vector[i];
    mean = sum / medidas;
    return mean;
}


double varianza(double vector[], int medidas) {
    double mean = media(vector, medidas);
    double sum_var = 0.0;
    for (int i = 0; i < medidas; i++) sum_var += (vector[i] - mean) * (vector[i] - mean);
    double std = sqrt(sum_var / (medidas - 1));
    return std;
}