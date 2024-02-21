


#include <iostream> // inputs and outputs
#include <cmath> //For math functions
#include <fstream>
#include <random>

using namespace std;

//prototipos de las funciones 
//minimal image conventio
void mic(float * vec, float Lbox){ //The asterisk provokes permanent changes in the variable we introduce inside
    for (int i = 0; i < 3; ++i){
        vec[i]-=floorf(0.5+vec[i]/Lbox)*Lbox; //This way if a particle is very far away we correct it as it is in other box, so rhis way is closer and we remove teh unreal effect caused by an imaginary box 
    }

}



//Global variables
float sigma = 1.0; 
float epsilon = 1.0;
float Lbox = 10.0; //Lbox of the box 
float volume= Lbox*Lbox*Lbox;
int dim = 3; //Number of dimensions


//In C, the function name must be declare before the main but its not neccesary to define it outside of the main

float energy (float position[], float Pos_itag[0], int itag, int Npart); //In C, since it is a compiled language, it doesn't go line by line but reads everything we compiling and then execute it




int main (int argc, char *argv[]){

    //Local variables

    float density=atof(argv[1]);
    float temp = atof(argv[2]); //Temperature of the system, set to 1 so it is easier. We are in the coexistence regime
    int Npart=(int)(density*volume); //Number of particles (int) finds the closest integer
    float dens= Npart/volume;

    //LOG of the run 
    cout << "Density" <<  " " << dens << endl;
    cout << "Temperature" << " " << temp << endl;


    random_device rand_dev;
    mt19937 generator(rand_dev());

    //system variables 


    float position[dim*Npart];


    //MC variables

    int Nstep=50000;
    int nsamp_ener=10; //The same as below but for the energy 
    int nsamp_pos= 100; //At how many n we are gonna sample the position
    int naccept = 0; //Number of time we are accepting
    float deltaR = 0.1; //Lbox of the jump


    //Data files

    ofstream fich_ener, fich_posi;

    fich_ener.open("energy.txt"); //We open two files for both energy and position
    fich_posi.open("position.txt");

    uniform_real_distribution<> dis1 (-1.0, 1.0);
    uniform_real_distribution<> dis2 (0.0, 1.0);
    uniform_int_distribution<> dist (0, Npart-1);

    //Initial positions 

    for (int i=0; i<3*Npart; ++i) //++i is faster than i++ because of some considerations of when it gives the number and things like that  3*Npart because we are going in the three dimensions
    {
        position[i]=dis1(generator)*Lbox/2; //I think this add the inital position for the particles. It assigns each particle three dimensions (that's why the cycle goes to 3*Npart) that goes between 5 to -5 for all three dimensions
    }

    //MonteCarlo loop

    for (int istep=0; istep<Nstep; ++istep){

        int itag = dist(generator); //Just inside this block so we avoid problems of definition, this the tag of each particle
        float posnew[dim], posold[dim]; //We define the positions and the energies, as well as the probabilities and the random numbers
        float enew, eold;
        float xi2;
        float prob;
        float Esample=eold;

        //Move particle itag)
        for (int k=0; k<dim; ++k){
            //We compute both old and new )position
            posold[k]=position[dim*itag+k]; 
            posnew[k]=position[dim*itag+k]+deltaR*dis1(generator); // dim*itag mark the beggining of the three dimensions of itag inside the position vector and k marks one of each dimension (0,1,2)
        }
        
        eold=energy(position, posold, itag, Npart); //We compute both old and new energy

        enew = energy(position, posnew, itag, Npart);

        prob=exp(-(enew-eold)/temp);

        xi2=dis2(generator);

        if (prob>xi2) //if accept 
        {
            for (int k=0; k<dim; ++k) position[itag*dim+k]=posnew[k]; //We update the position
            naccept++; //We add one to the accept count 
            Esample=enew; //We change Esample to accumulate the energy if required
        }

        if (istep % nsamp_ener==0) fich_ener << 0.5*Esample << endl; // If we put everything in one line the {} are not needed
        if (istep % nsamp_pos==0) {
            for (int k=0; k<dim; k++){
                fich_posi << position[itag*dim+k] << endl;
            }
        }

    }

    fich_ener.close();
    fich_posi.close();
}

//Function that computes the energy
float energy(float position[], float pos_itag[], int itag, int Npart){

    float sigma=1;
    float epsilon=1;
    float cutoff=2.5*sigma;
    float pow_cutoff=cutoff*cutoff;
    float U=0.0;

    for (int j = 0; j<Npart; ++j){
        if (j == itag) continue;
        float distance_vector[dim];
        float U_iJ = 0.0;
        
        for (int k=0;k<dim; ++k) distance_vector[k]= pos_itag[k]-position[dim*j+k];
        
        mic(distance_vector, Lbox);
        float distance_sq= 0.0;
        for (int k=0;k<dim; ++k) distance_sq += pow(distance_vector[k],2);

        if (distance_sq < pow_cutoff){
            float distance = sqrt(distance_sq);
            float dist6 = pow ((sigma/distance), 6.0);
            float cutoff6= pow ((sigma/cutoff), 6.0);

            U_iJ += 4*epsilon*((dist6*dist6-dist6) - (cutoff6*cutoff6-cutoff6));
        }
        U+= U_iJ;
    }
    return U;
}

// cout displays what we give it in the terminal. Is a good way of debugging
// cout << Eold << " " << istep << ... << endl
/* this comments a block */
