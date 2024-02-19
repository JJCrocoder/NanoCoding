#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
using namespace std;

// Declaración de la función de la convención de imagen mínima
void mic(float * vec, float Lbox) {
    // Implementación de la convención de imagen mínima
    for (int i = 0; i < 3 ; ++i) {
        vec[i] -= floorf(0.5 + vec[i]/Lbox)*Lbox;
    }  
    return;
}

// Prototipo de la función que calcula la energía
float Energy(float Position[], float Pos_itag[], int itag, int Npart);

// Variables globales
float sigma = 1.0;
float epsilon = 1.0;
float Lbox = 10.0;
float volume = Lbox*Lbox*Lbox;
int dim = 3;
float r_cut = 2.5 * sigma;

int main(int argc, char *argv[]) {
    // Verificar si se proporcionan los argumentos correctos en la línea de comandos
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <density> <temperature>" << endl;
        return 1;
    }

    // Inicialización del generador de números aleatorios
    random_device rand_dev;
    mt19937 generator(rand_dev());

    // Obtener los parámetros de densidad y temperatura de la línea de comandos
    float densinput = atof(argv[1]);
    float Temp = atof(argv[2]);

    // Calcular el número de partículas
    int Npart=(int)(densinput*volume);
    float dens = Npart/volume;

    // Registro de la densidad y la temperatura
    cout <<"Density "<< dens <<endl;
    cout <<"Temperature "<< Temp <<endl;

    // Variables del sistema
    float Position[dim*Npart];

    // Variables del algoritmo de Monte Carlo
    int Nstep = 500000;
    int nsamp_ener = 10;
    int nsamp_pos = 100;
    int naccept = 0;
    float deltaR = 0.1;

    // Archivos de datos
    ofstream fich_ener, fich_posi;
    fich_ener.open("energy.txt");
    fich_posi.open("position.txt");

    // Generadores de números aleatorios
    uniform_real_distribution<> dis1(-1.0, 1.0);
    uniform_real_distribution<> dis2(0.0, 1.0);
    uniform_int_distribution<> dist(0, Npart-1);

    // Posiciones iniciales
    for (int i=0; i<3*Npart; ++i) Position[i] = dis1(generator)*Lbox/2.;

    // Bucle del algoritmo de Monte Carlo
    for (int istep = 0; istep<Nstep; ++istep) {
        int itag = dist(generator);	
        float PosNew[dim], PosOld[dim];
        float Enew, Eold, prob;

        // Mover la partícula itag
        for (int k=0; k<dim; ++k) {
            PosOld[k] = Position[dim*itag + k];
            PosNew[k] = Position[dim*itag + k] + deltaR * dis1(generator);
        }
        Eold = Energy(Position, PosOld, itag, Npart);
        Enew = Energy(Position, PosNew, itag, Npart);
        float Esample = Eold;

        // Ratio de probabilidad
        prob = exp(-(Enew-Eold)/Temp);

        // Generar un número aleatorio entre 0 y 1
        float xi = dis2(generator);
        // Si se acepta
        if (prob > xi) {
            for (int k=0; k<dim; ++k) Position[itag*dim + k] = PosNew[k];
            naccept = naccept + 1;
            Esample = Enew;
        }

        // Guardar la energía muestreada
        if (istep % nsamp_ener == 0) fich_ener << 0.5*Esample << endl;
       
        // Guardar las posiciones muestreadas
        if (istep % nsamp_pos == 0) {
            for (int i = 0; i < Npart; ++i) {
                for (int k = 0; k < dim; ++k) {
                    float pos_value = Position[i * dim + k];
                    fich_posi << pos_value << " ";
                }
                fich_posi << endl;
            }
        }
    }
    
    fich_ener.close(); // Cerrar el archivo de energía
    fich_posi.close(); // Cerrar el archivo de posición

    return 0;
}

// Función que calcula la energía
float Energy(float Position[], float Pos_itag[], int itag, int Npart) {
    // Calcular la energía total
    float U_tot = 0.0;
    float r_2cut = r_cut * r_cut;

    for ( int i = 0; i < Npart; ++i)  {
        if (i == itag) continue;
        float vec_dist[dim];
        float u_ij = 0.0;

        // Calcular la distancia mínima entre partículas
        for (int k = 0; k < dim; ++k) vec_dist[k] = (Pos_itag[k] - Position[i * dim + k]);        
        mic(vec_dist, Lbox);        
        float r2_ij=0.0;
        for (int k = 0; k < dim; ++k) r2_ij += pow(vec_dist[k], 2); 

        // Calcular el potencial de Lennard-Jones
        if (r2_ij < r_2cut) {
            float r_mod = sqrt(r2_ij);
            float r6 = pow((sigma/r_mod), 6.0);
            float rc6 = pow((sigma/r_cut), 6.0);   
            u_ij = 4 * epsilon * ((r6*r6-r6) - (rc6*rc6-rc6));   
        }
        U_tot += u_ij; // Sumar la contribución de esta partícula a la energía total
    }
    return U_tot; // Devolver la energía total  
}
