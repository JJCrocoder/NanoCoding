#include<vector> //vector class

using namespace std;

//function declaration
int getcell(float *pos, int *cell);
     //Fills head and list with the current system state
void make_linked_list(float *pos, vector<int> list, vector<int> head, int Npart); 
void heal_list(int cella, int cellb); //Recalculates cella and cellb in head and list
void mic(float * vec, float Lbox);
int pbc_cells(int icell);

// We calculate the number of cells in each dimension
float cellsize = 1.15 * r_cut;		// A bit bigger than r_cut
int ncells = (int)(Lbox/cellsize);	// we obtain the number of cells

float mcells = ncells*ncells*ncells;	// Now the total number of cells

void make_linked_list(float *pos, vector<int> list, vector<int> head, int Npart){
  //Neighbour list head and list creation
  for(int i=0; i<Npart; ++i) list[i] = -1;
  for(int icel=0; icel<mcells; ++icel) head[icel] = -1;
  /* build head and list */

  for(int i=0; i<Npart; ++i){   // We go trough all the particles
    float temppos[dim] = {0};  // For each particle, we initialize an auxiliar position vector 
      for(int k=0; k<dim; ++k) temppos[k] = pos[dim*i+k]; // We store the coordinates in the auxiñiar vector 
      mic(temppos, Lbox);    // Fold the position of particle i using mic
    int cell[dim];
    int icell = getcell(temppos, cell); // Find the cell where temppos belongs to
    // build the list and head arrays
    list[i] = head[icell];   // The list element for the current particle points to the last saved particle of each cell
    head[icell] = i;         // Now the last saved particle will be the current one
  }
}

//function creation
int getcell(float *normpos, int *cell){ // IMPORTANT: We use as an argument the normalized position by the box length
    //since we are dividing the space in cubes, the ncells in each direction is the same
    for (int k = 0 ; k < dim; ++k) cell[k]= (int)((0.5 + normpos[k])*ncells);
    int icell = cell[0]+cell[1]*ncells+cell[2]*ncells*ncells ;
    return icell;
}

//implementation in the code of the logic above

/* Compute the energy of the itag particle when it is located in pos=postag */
//Use head and list neigbour search****//
 void forceNL(float *pos, float *Forces, int Npart){
   
   // Clear Forces
   for(int i=0; i<dim*Npart; ++i){
     Forces[i]=0.0;
   }

   vector<int> list, head;
   make_linked_list (pos, list, head);

   // We get the cell for every particle
   for(int i = 0; i<Npart; ++i){
     int cell[dim]={0};
     int icell=getcell(&pos[3*i], cell); // Using "&" gives us the memory direction (pointer), so we work with the location

    /*For every neighbouring cell (27 cells in 3D)*/ // THERE SHOULD BE A MORE EFFICIENT WAY
    for(int jx=cell[0]-1; jx<=cell[0]+1;jx++){
        for(int jy=cell[1]-1; jy<=cell[1]+1;jy++){
            for(int jz=cell[2]-1; jz<=cell[2]+1;jz++){

	      int j = head[pbc_cells(jx) + pbc_cells(jy)*ncells + pbc_cells(jz)*ncells*ncells];
               /*If there is no particles go to the next cell*/
               if(j==-1) continue;
               /*Use list to travel through all the particles, j, in cell jcel*/
               do{
                   if(i<j)
                   // CALCULATE PAIRINTERACTION HERE: THE ENERGY UIJ AND SUM UP FOR UTOT //
                   j=list[j];
		   // ADD FORCE FUNCTION: remind that we must apply pbc (mic)
               } while (j!=-1);

            }
        }
    }
     
   }

}

// Give me the cell index in head taking into account PBC
int pbc_cells(int icell){
if(icell==-1) return ncells-1;
else if(icell==ncells) return 0;
else return icell;
}
