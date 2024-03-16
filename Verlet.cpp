#include <iostream>
#include<fstream>
#include<random>
#include<cmath>
#include<stdlib.h>
#include<algorithm>
#include<vector>

using namespace std;
const double pi = 3.14159265358979323846;





//Definimos los generadores de números aleatorios.

//Generador con probabilidad uniforme de números double entre un max y un min dados.
double uniform(double min, double max) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(min, max);

    return dis(gen);
}

//Generador con probabilidad uniforme de números int entre un max y un min dados.
int rand_num(int min, int max) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<int> dis(min, max);

    return dis(gen);
}

//Generador con probabilidad Gaussiana con media mean y varianza stddev.
double rand_normal(double mean, double stddev) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dis(mean, stddev);
    return dis(gen);
}




//Declaramos las funciones a emplear en el código.
double Energy(double Position[], double N_O_Pos[], int itag, double r_cut);
double LJ_pot(double r, double r_cut);
vector<double> Force_over_particle(int N_part, double Position[]);
double MIC_distancia(double Position_part_i[], double Postion_part_j[]);
double Energy_Brownian(double PosNew[]);
void Move_Verlet(int N_part, double Position[], double Velocities[], int sw);
double Temperature(double Velocities[]);


const double temp = 1.3;         // Temperatura del sistema.
const double L_box = 10;       // Lado de la caja.


const int N_part = 1000; //Número de partículas del sistema.



const int N_dim = 3; // Número de dimensiones del sistema.

double dens = N_part / pow(L_box, 3); //Densidad de partículas en el sistema.


//Constanttes para construir el potencial y fuerzas de Leonardo-Cojones
const double epsilon = 1.0;
const double sigma = 1.0;
const double r_cut = 2.5 * sigma;       // Cut-off del potencial (Mierdón sideral)

//Constantes para el movimiento Browniano (El movimiento Browniano consiste en que las partículas me la agarran con la mano).
const double D = 1;  //Coeficiente de difusión (Creo que así se llama).
const double mu = D / temp; //Movilidad (Se mueven las partículas lo mismo que Koke con el pedazo de culo de dominicana que lleva años cebando)
const double dt = 0.005; //Intervalo temporal (Para todos iguales menos para los del Depor que para ellos seguimos estando en 2001).

// Variables para el Monte Carlo Abuelotti y el movimiento Browniano.
const int N_step = 100000;     // Número de pasos en la simulación de Monte Carlo.
int Brow_steps = 1000;    // Número de saltos con movimiento Browniano (Tardan lo mismo en darse todos que lo que tarda Koke en bajar a defender).
int nsamp_pos = 100;    // Ratio con el que tomamos las medidas del Browniano.
double deltaR = 0.2;     // Tamaño del salto en el Monte Carlo (Igual al número de puntos por partido que tiene el Almería esta temporada).
int tsteps = Brow_steps / nsamp_pos;  //Esto es para la normalización de la RDF. Es el tiempo que hay entre medida y medida. 



//Variables para el histograma
const int num_bins = 120
; //Número de bins (bin=papelera para aquellos que no entiendan el catalán).
double min_val = 0.0;  //Valor mínimo del histograma (Puntos totales que me ha hecho el hijo de puta de Balliu en 4 jornadas jugando el cabrón todos los minutos).
double max_val = L_box / 2;  //Valor más alto del histograma.
double rango = max_val - min_val; //Rango sobre el que construimos el histrograma.
double DeltaX = rango / num_bins;  //Anchura de los bins (Igual que el tamaño del culo de la vacaburra de Koke).
vector<int> Histo(num_bins, 0);  //Vector cuyos elementos estan asociados con los distintos bins del histograma.
double Const = 4.0 / 3.0 * dens * pi; //Constante necesaria para la normalización (Mierdón).





//Main code
int main() {
    cout << "El dt es: " << dt << endl; //Tamaño del salto temporal

    double dens = N_part / pow(L_box, 3);
    cout << "Densidad: " << dens << endl; //Densidad de partículas
    cout << "Temperatura: " << temp << endl; //Temperatura del sistema
    double Position[N_dim * N_part]; // Array donde vienen las posiciones de las partículas (array tipo GOD)



    // Definimos la faquin configuración inicial
    for (int k = 0; k < N_dim * N_part; ++k) {
        Position[k] = uniform(-0.5 * L_box, 0.5 * L_box); // Coordenadas random en nuestra caja.
    }




    // Bucle principal del Monte Carlo (33 negro impar y pasa). No lo voy a explicar mucho que me da pereza (Grupo de música TOP).

    for (int istep = 0; istep < N_step; ++istep) {
        int itag = rand_num(0, N_part);
        double PosNew[N_dim], PosOld[N_dim];
        double Enew, Eold;

        for (int k = 0; k < N_dim; ++k)
        {
            PosOld[k] = Position[N_dim * itag + k];
            PosNew[k] = Position[N_dim * itag + k] + uniform(-deltaR, deltaR);
        }
        Eold = Energy(Position, PosOld, itag, r_cut);
        Enew = Energy(Position, PosNew, itag, r_cut);




        // Ratio de (masonidad) probabilidad

        double prob = exp(-(Enew - Eold) / temp);
        double xi = uniform(0, 1);
        double Esample = Eold;

        if (prob > xi) {
            for (int k = 0; k < N_dim; ++k) {
                Position[itag * N_dim + k] = PosNew[k];
            }
            Esample = Enew;
        }
    }
    double Pos_New[N_part * N_dim];
    double Pos_Old[N_part * N_dim];
    double Velocities_Old[N_dim * N_part];
    double Velocities_New[N_dim * N_part];

    // Copiar el contenido de Position a Pos_New y Pos_Old para el Brownian.
    for (int i = 0; i < N_part * N_dim; ++i) Pos_New[i] = Position[i];

    //Generamos las velocidades iniciales
    
    for (int k = 0; k < N_dim * N_part; k++) Velocities_New[k] = rand_normal(0,sqrt(temp));




    cout << "Monte Carlo corrido" << endl; //Avisa de que el Monte Carlo a terminado de correrse (Vaya puto guarro, yo de vosotros le hacia un test de embarazo a vuestra terminal).
    //Abrimos el fichero que vamos a usar
    ofstream fich_rdf;



    fich_rdf.open("Brownian_RDF.txt"); //File for saving the histogram


    //Dinámica Browniana, la chicha de este programa (chicha=90% del cuerpo de Koke).
    int sw ;

    //Bucle Browniano.
    for (int bstep = 0; bstep < Brow_steps; ++bstep) {
        for (int i = 0; i < N_part * N_dim; ++i) {
            Pos_Old[i] = Pos_New[i];
            Velocities_Old[i] = Velocities_New[i];

        }

        sw = 1;
        Move_Verlet(N_part, Pos_New, Velocities_New, sw);

        

        
        
        
        
               
        
        sw = 2;
        Move_Verlet(N_part, Pos_New, Velocities_New, sw);
        















        //Sampleo de posiciónes. Lo de la RDF de toda la vida.
        if (bstep % nsamp_pos == 0) {
            //We calculate the distances between every pair of particles
            double Temp = Temperature(Velocities_New);
            cout << "La temperatura es ahora: "<<Temp << endl;
            for (int p = 0; p < N_part - 1; ++p) {
                for (int q = p + 1; q < N_part; ++q) {
                    double Pos_i[3];
                    double Pos_j[3];

                    for (int k = 0; k < 3; ++k) {
                        Pos_i[k] = Pos_New[q * N_dim + k];
                        Pos_j[k] = Pos_New[p * N_dim + k];
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


    
    // Normalización del histograma.
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



double Energy(double Position[], double N_O_Pos[], int itag, double r_cut) {
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

    return  U_itag;
}

vector<double> Force_over_particle(int N_part, double Position[]) {
    vector<double> F_itag(N_part * N_dim, 0.0);

    // Bucle sobre todas las partículas
    for (int i = 0; i < N_part - 1; i++) {
        for (int j = i + 1; j < N_part; j++) {
            // Construcción del vector de posición relativa y su módulo

            double rij[N_dim];   // Vector de posición relativa
            double mod_rij;      // Módulo del vector

            // Realizamos la resta y calculamos el módulo
            for (int k = 0; k < N_dim; ++k) {
                rij[k] = Position[j * N_dim + k] - Position[i * N_dim + k];
                rij[k] -= L_box * floor(rij[k] / L_box + 0.5);

            }



            mod_rij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];



            // Calculamos la fuerza y la añadimos al vector de fuerzas
            double factor = -24.0 * epsilon / sigma * (2 * pow((sigma / mod_rij), 13) - pow((sigma / mod_rij), 7)) / mod_rij;

            // Corrección de índices aquí
            for (int k = 0; k < N_dim; ++k) {
                F_itag[i * N_dim + k] += factor * rij[k];
                F_itag[j * N_dim + k] += -factor * rij[k];
            }


        }
    }
    return F_itag;
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


double MIC_distancia(double Position_part_i[], double Position_part_j[]) {
    double rij[3];
    for (int k = 0; k < 3; ++k) {
        rij[k] = Position_part_j[k] - Position_part_i[k];
        rij[k] -= L_box * floor(rij[k] / L_box + 0.5);
    }
    double mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
    return mod_rij;
}


double Energy_Brownian(double PosNew[]) {
    double U_itag = 0;  // The energy of the selected particle
    // Loop for the particles
    double rij[N_dim];
    double mod_rij;
    for (int i = 0; i < N_part - 1; ++i) {
        for (int j = i + 1; j < N_part; j++) {
            for (int k = 0; k < N_dim; k++) {
                rij[k] = PosNew[j * 3 + k] - PosNew[i * 3 + k];
                rij[k] -= L_box * floor(rij[k] / L_box + 0.5);

            }
        }

        mod_rij = sqrt(rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]);
        U_itag += LJ_pot(mod_rij, r_cut);


    }

    return  U_itag;
}



void Move_Verlet(int N_part, double Position[], double Velocities[], int sw) {
    
    if (sw = 1) {
        vector<double> Forces_N=Force_over_particle(N_part, Position);
        for (int k = 0; k < N_dim * N_part; k++) {
            Velocities[k] = Velocities[k] + Forces_N[k] * dt / 2;
            Position[k] = Position[k] + Velocities[k] * dt;
        }
    }
    else {
        vector<double> Forces_N= Force_over_particle(N_part, Position);
        for (int k = 0; k < N_dim * N_part; k++) {
            Velocities[k] = Velocities[k] + Forces_N[k] * dt / 2;
        }
    }

}


double Temperature(double Velocities[]) {
    double Temp=0.0;
    
    for (int k = 0; k < N_part*N_dim; k++) Temp += Velocities[k] * Velocities[k];   
    
    
    
    Temp /= 3 * N_part;
    return Temp;
}