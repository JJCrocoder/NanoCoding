

/*
For this implementation, what i think i should do is once i have the positions of all molecules,
check using the histogram logic to which cell they belong. Once I have that, I should implement that
the function energy only is executed on those molecules that are on neighbouring cells */


//LETS START DOING THE LINKED CELLS

#include<vector> //vector class

using namespace std;

/*this is the logic for the generation of the linked cells. I have to create the apply_pbc 
(I think mic does something similar) and the getcell function*/

//function declaration
int getcell(float *pos);

int dim = 3;

//Following C++ logic (what i understand) is better to implement this as part of the code instead of doing a function
void make_linked_list(); //Fills head and list with the current system state
void heal_list(int cella, int cellb); //Recalculates cella and cellb in head and list

//Neighbour list head and list creation
vector<int> list, head;
for(int i=0, i<natoms, ++i) list[i] = -1;
for(int icel=0, icel<ncel, ++icel) head[icel] = -1;

/* build head and list */
for(int i=0, i<natoms, ++i){   // We go trough all the particles
    float temppos[dim] = {0};  // For each particle, we initialize an auxiliar position vector 
      for(int k=0; k<dim; ++k) temppos[k] = pos[dim*i+k]; // We store the coordinates in the auxiñiar vector 
    mic(temppos);    // Fold the position of particle i using mic
    icell = getcell(temppos); // Find the cell where temppos belongs to
    /*build the list and head arrays*/
    list[i] = head[icell];   // The list element for the current particle points to the last saved particle of each cell
    head[icell] = i;         // Now the last saved particle will be the current one
}

//function creation
int getcell(temppos){
    //initialization of variables 
    int icell_dim[dim];
    for (int k = 0 ; k < dim; ++k) icell_dim[k]= (int)(1 + (0.5 + temppos[k])*ncells); //since we are dividing the space in cubes, the ncells in each direction is the same
    int icell = icell_dim[0]+(icell_dim[1]-1)*ncells+(icell_dim[2]-1)*ncells*ncells ;
    return icell
}

//implementation into the code od teh logic above

//** Compute the energy of the itag particle when it is located in pos=postag**//
//Use head and list neigbour search****//
void forceNL(int itag, float postag[3]){
    /*Get the position to the primary box*/
    mic(postag);
    int icel = getcell(postag);
    int jcel =get_N_cells(postag); /*For me this could be the best approach: Define a function that takes
    the postag and gives you a list of dimension 6 with the ID of the neighbouring cells.  
    What about the cells at the extremes?, PBC?*/
    for (int i=0; i < 6 ; ++i){ 
        j = head[jcel[i]];
        /*If there is no particles go to the next cell*/
        while (j!= 0){
            if(itag!=(j-1)){
            float tempposj[dim]={0.0};
            for (int k = 0; k< dim ; ++k) tempposj[k]=position[dim*j+k];
            //apply energy calculations here
            //
            j=list[j]

            }
            // CALCULATE PAIRINTERACTION HERE: THE ENERGY UIJ AND SUM UP FOR UTOT //
            j=list[j];
        /*When j=0 (list[j] = 0) then there is no more particles in cell jcel (see the notes!)*/
        }
    }
}
//**Give me the cell index in head taking into account PBC**//
int pbc_cells(int icell){
if(icell==0) return ncells;
else if(icell==ncells+1) return 1;
else return icell;
}
