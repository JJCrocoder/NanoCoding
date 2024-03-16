

/*for this implementation, waht i think i should do is once i have the positions of all molecules,
check using the histogrma logic to which cell they belong. Once I have that, I should implement that
the function energy only is executed on those molecules that are on neighbouring cells */

//LETS START DOING THE LINKED CELLS


/*this is the logic for the generation of the linked cells. I have to create the apply_pbc 
(I think mic does something similar) and the getcell function*/
#include<vector> //vector class
//function declaration

int getcell(float *pos)

//Following C++ logic (what i understand) is better to implement this as parto of the code instead of doing a function
void make_linked_list(); //Fills head and list with the current system state
void heal_list(int cella, int cellb); //Recalculates cella and cellb in head and list
vector<int> list, head; //Neighbour list head and list
for(int i=0, i<natoms, i++) list[i]=0
for(int icel=0, icel<ncel, icel++) head[icel]=0
/* build head and list */
for(int i=1, i<=natoms,i++){
    float temppos[3];
    /* the particle index is i-1, but we work with i= 1 to N to build head array*/
    for(int k=0; j<3; k++) temppos[j] = pos[3*(i-1)+k]; /*position of particle i-1*/
    mic(temppos); /* fold the position of particle i into the primary box */
    icell = getcell(temppos); /*find the cell where temppos resides*/
    /*build the list and head arrays*/
    list[i] = head[icell]; 
    head[icell] = i;
}

//function creation
int getcell(temppos){
    //initialization of variables 
    int icell_dim[dim];
    for (int k = 0 ; k < dim; ++k) icell_dim[k]= (int)(1 + (0.5 + temppos[k])*ncells); //since we are dividing the space in cubes, the ncells in each direction is the same
    int icell = icell_dim[0]+(icell_dim[1]-1)*ncells+(icell_dim[2]-1)*ncells*ncells ;
    return icell
}


int get_N_cell(temppos){
    int icell_dim[dim];
    for (int k = 0 ; k < dim; ++k) icell_dim[k]= (int)(1 + (0.5 + temppos[k])*ncells);
    for (int i=0; i < 6 ; ++i ){
        int jcell_dim[dim]={0}
        
    }

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




