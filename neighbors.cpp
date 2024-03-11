#include<vector> //vector class

using namespace std;

/*this is the logic for the generation of the linked cells. I have to create the apply_pbc 
(I think mic does something similar) and the getcell function*/

//function declaration
int getcell(float *pos, int *cell);
//Following C++ logic (what i understand) is better to implement this as part of the code instead of doing a function
void make_linked_list(); //Fills head and list with the current system state
void heal_list(int cella, int cellb); //Recalculates cella and cellb in head and list
void mic(float * vec, float Lbox);
int pbc_cells(int icell);

int dim = 3;
int ncells = 20;
float Lbox = 10.0;

mcells = ncells*ncells*ncells;

vector<int> list, head;

void make_linked_list(float *pos){
//Neighbour list head and list creation
for(int i=0, i<natoms, ++i) list[i] = -1;
for(int icel=0, icel<mcells, ++icel) head[icel] = -1;

/* build head and list */
for(int i=0, i<natoms, ++i){   // We go trough all the particles
    float temppos[dim] = {0};  // For each particle, we initialize an auxiliar position vector 
      for(int k=0; k<dim; ++k) temppos[k] = pos[dim*i+k]; // We store the coordinates in the auxiñiar vector 
    mic(temppos);    // Fold the position of particle i using mic
    icell = getcell(temppos); // Find the cell where temppos belongs to
    // build the list and head arrays
    list[i] = head[icell];   // The list element for the current particle points to the last saved particle of each cell
    head[icell] = i;         // Now the last saved particle will be the current one
}

//function creation
int getcell(float *normpos, int *cell){ // IMPORTANT: We use as an argument the normalized position by the box length
    //since we are dividing the space in cubes, the ncells in each direction is the same
    for (int k = 0 ; k < dim; ++k) cell[k]= (int)(1 + (0.5 + normpos[k])*ncells);
    int icell = cell[0]+(cell[1]-1)*ncells+(cell[2]-1)*ncells*ncells ;
    return icell
}

//implementation in the code of the logic above

//** Compute the energy of the itag particle when it is located in pos=postag**//
//Use head and list neigbour search****//
void forceNL(int itag, float *postag){
    /*Get the position to the primary box*/
    mic(postag, Lbox);

    int celli[dim], cellj[dim];
    // IMPORTANT: We need to include the normalized position array
    int icel = getcell(postag, celli);
    int j; //Index of neighbour particle
    int jcel, jcelx, jcely, jcelz; //cell coordinates and cell index for particle j
      
    /*For every neighbouring cell (27 cells in 3D)*/ // THERE SHOULD BE A MORE EFFICIENT WAY
    for(int jx=cell[0]-1; jx<=cell[0]+1;jx++){
        for(int jy=cell[1]-1; jy<=cell[1]+1;jy++){
            for(int jz=cell[2]-1; jz<=cell[2]+1;jz++){
               /*The neighbour cell must take into account pbc! (see pbc_cells!)*/
               cellj[0] = pbc_cells(jx);
               cellj[1] = pbc_cells(jy);
               cellj[2] = pbc_cells(jz);
               //See getcell!
               jcel = cellj[0]+(cellj[1]-1)*ncells+(cellj[2]-1)*ncells*ncells ;
               /*Get the highest index particle in cell jcel*/
               j = head[jcel];
               /*If there is no particles go to the next cell*/
               if(j==-1) continue;
               /*Use list to travel through all the particles, j, in cell jcel*/
               do{
                   if(itag!=j)
                   // CALCULATE PAIRINTERACTION HERE: THE ENERGY UIJ AND SUM UP FOR UTOT //
                   j=list[j];
               } while (j!=-1);

            }
        }
    }
}

// Give me the cell index in head taking into account PBC
int pbc_cells(int icell){
if(icell==0) return ncells;
else if(icell==ncells+1) return 1;
else return icell;
}

void mic(float * vec, float Lbox) {
    for (int k = 0; k < dim ; ++k) vec[k] -= floorf(0.5 + vec[k]/Lbox)*Lbox;
    return;
}
