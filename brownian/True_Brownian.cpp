
#include <iostream>
#include <cmath>
#include <fstream>
#include <random>

using namespace std;

// Prototipos de las funciones

/* Minimum image convention */
void mic(float * vec, float Lbox){

  for (int i = 0; i < 3 ; ++i) {
    vec[i] -= floorf(0.5 + vec[i]/Lbox)*Lbox;
  }
  
  return;
}

float Energy(float Position[], float N_O_Pos[], int itag);

// Global variables

float sigma = 1.0;
float epsilon = 1.0;
float Lbox = 10.0;
float densinput = 0.7;
float volume = Lbox*Lbox*Lbox;
int Npart=(int)(densinput*volume);
int dim = 3;
float Temp=1.9;
float r_cut = 2.5 * sigma;
float DIFF =1.0;



int Nstep = 50000;
int nrdf = 10;
int ndraw = 100;
float deltat = 0.1;


int main (void)
{

  float dens = Npart/volume;
  // LOG of the run
  cout <<"Density "<< dens <<endl;
  cout <<"Temperature "<< Temp <<endl;
  // system variables
  
  float Pos[dim*Npart];
  float For[dim*Npart];
  
  
  // Data files
  
  ifstream fich_pos;
  fich_pos.open("position.txt");
  // Initial positions: prepared from monte carlo: read pos
  fich_pos>>Npart;
  for (int i=0; i<3*Npart; ++i) fich_pos>>Pos[i];
  fich_pos.close();
  
  //+++++++++++++++++++++++++++++++++++++
  
  /// START BD LOOP
  
  for (int istep = 0; istep<Nstep; ++istep)
    {
      
      // calculate forces : For
      forces(Pos,For);
      // Move particles
      moveEuler(Pos,For,deltat);
      
      //sample things:
      if (istep % nrdf == 0) {
	rdf(1,pos);
      }
      if (istep % ndraw == 0) {
	draw(pos);
      }      
    }      
    return 0;
}


// Funcion que calcula la energia

float Energy(float Position[], float Pos_itag[], int itag)
{
  float U_tot = 0.0;
  float r_2cut = r_cut * r_cut;

  for ( int i = 0; i < Npart; ++i)  {
    if (i == itag) continue;
    float vec_dist[dim];
    float u_Ij = 0.0;
    
    for (int k = 0; k < dim; ++k)
      vec_dist[k] = (Pos_itag[k] - Position[i * dim + k]);
    
    mic(vec_dist, Lbox);

    /*    cout << vec_dist[0]<<" "<<vec_dist[1]<<" "<<vec_dist[2]<<endl;*/
      
    float r2_Ij=0.0;
    for (int k = 0; k < dim; ++k)
      { r2_Ij += pow(vec_dist[k], 2);   
	/*	cout << r2_Ij<<endl;*/
      }
    //lennard-jones potential
	
    if (r2_Ij < r_2cut)
      {
	float r_mod = sqrt(r2_Ij);
	float r12 = pow((sigma/r_mod), 12.0);
	float r6 = pow((sigma/r_mod), 6.0);
	float rc12 = pow((sigma/r_cut), 12.0);
	float rc6 = pow((sigma/r_cut), 6.0);
	
	u_Ij = 4 * epsilon * (r6*(r6-1) - rc6*(rc6-1));   
      }
    U_tot += u_Ij;
      }
  return U_tot;  
}

/*La fuerza es vectorial, tiene tres dimensiones. Podemos calcular la fuerza en cada dimensión mediante la energía en 
cada dimensión ?

Para ello multiplicamos la energía por cada interacción pairwise en cada dimensión 
por distanciaX/distancia en cada coordenada

*/ 
float forces (float Position[Npart*dim], int itag){
/*Aquí dentro: Calculamos las energías con una función similar a la de la energía quedandonos con lo de u_Ij
Después multiplicamos cada u_Ij por vec_dist[x,y o z]/distancia y acumulamos con un loop for dentro de una variable que 
se llame f_iJ. Cuando se hayan dado todas las pairwise interactions, sale un f_iJ que lo asociamos dentro de un vector
de fuerzas a la posición del itag*/
    float force[Npart*dim] //Definimos aquí el vector donde vamos a acumular las fuerzas
    float r_cut=2*sigma
    for (i=0; i < Npart; ++i){
        float force_iJ[dim]={0.0}; //Definimos la variable donde vamos a guardar los resultados pairwise
        float Pos_itag[dim]={0.0};
        for (int k=0; k< dim; ++k) float Pos_itag[k]=Position[itag*dim+k];
        for (j=i+1; j< Npart-1; ++j){
            //Definimos algunas variables
            float r_2cut = r_cut * r_cut;
            float r2_Ij=0.0;
            float vec_dist[dim];
            float u_Ij = 0.0; //esta es la energia y encima el vector distancia (useful)
            for (int k = 0; k < dim; ++k){
                vec_dist[k] = (Pos_itag[k] - Position[i * dim + k]);
                mic(vec_dist, Lbox);
                r2_Ij += pow(vec_dist[k], 2); 
            }  
            //lennard-jones potential
            
            if (r2_Ij < r_2cut){
                float r_mod = sqrt(r2_Ij);
                float r6 = pow((sigma/r_mod), 6.0);
                float rc6 = pow((sigma/r_cut), 6.0);
                
                u_Ij = 4 * epsilon * (r6*(r6-1) - rc6*(rc6-1));

                for (int k =0;k<dim;++k) force_Ij[k] += -u_Ij*vect_dist[k]/r_mod;                
            }
        }
        for (int k =0;k<dim;++k) force[i*dim+k] = force_Ij[k];                

    }
}


// random_device rand_dev;
//  mt19937 generator(rand_dev());

// create forces
// create move
