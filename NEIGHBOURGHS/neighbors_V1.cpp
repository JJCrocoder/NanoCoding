#include <vector>
#include <cmath>

const int dim = 3; // Dimensionality of the system
const int ncells = 10; // Number of cells in each dimension
const int Npart = 100; // Number of atoms

std::vector<int> list(natoms + 1, 0); // Neighbor list
std::vector<int> head(ncells * ncells * ncells + 1, 0); // Head list for cells

void mic(float pos[dim]) {
    // Apply minimum image convention
    for (int k = 0; k < dim; ++k) {
        pos[k] -= ncells * floor(pos[k] / ncells);
    }
}

int getcell(float pos[dim]) {
    // Get the cell index for a given position
    int icell_dim[dim];
    for (int k = 0; k < dim; ++k) {
        icell_dim[k] = static_cast<int>(1 + (0.5 + pos[k]) * ncells);
    }
    return icell_dim[0] + (icell_dim[1] - 1) * ncells + (icell_dim[2] - 1) * ncells * ncells;
}

void make_linked_list(float* pos) {
    // Build head and list arrays
    for (int i = 1; i <= Npart; ++i) {
        float temppos[dim];
        for (int k = 0; k < dim; ++k) {
            temppos[k] = pos[dim * (i - 1) + k]; // Position of particle i-1
        }
        mic(temppos); // Fold the position into the primary box
        int icell = getcell(temppos); // Find the cell where temppos resides

        // Build the list and head arrays
        list[i] = head[icell];
        head[icell] = i;
    }
}

std::vector<int> get_N_cells(float pos[dim]) {
    // Get the neighboring cells for a given position
    int icell = getcell(pos);
    std::vector<int> jcells;

    // Add neighboring cells considering periodic boundary conditions
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                int jcell = pbc_cells(icell + dx + dy * ncells + dz * ncells * ncells);
                jcells.push_back(jcell);
            }
        }
    }
    return jcells;
}

void forceNL(int itag, float postag[dim], float* position) {
    // Compute the force of the itag particle when it is located at postag
    mic(postag);
    int icell = getcell(postag);
    std::vector<int> jcells = get_N_cells(postag);

    // Loop through neighboring cells
    for (int jcel : jcells) {
        int j = head[jcel];
        // Iterate over particles in the neighboring cell
        while (j != 0) {
            if (itag != (j - 1)) {
                float tempposj[dim];
                for (int k = 0; k < dim; ++k) {
                    tempposj[k] = position[dim * (j - 1) + k];
                }

                float r_2cut = r_cut * r_cut;
                float r2_Ij = 0.0;
                float vec_dist[dim];
                float u_Ij = 0.0;

                for (int k = 0; k < dim; ++k) {
                    vec_dist[k] = Pos_itag[k] - temppos[k];
                    mic(vec_dist, Lbox);
                    r2_Ij += pow(vec_dist[k], 2);
                }

                if (r2_Ij < r_2cut) {
                    float r_mod = sqrt(r2_Ij);
                    float r6 = pow((sigma / r_mod), 6.0);
                    float rc6 = pow((sigma / r_cut), 6.0);
                    float f_ij = -48 * epsilon * r6 * (r6 - 0.5) / r_mod;

                    for (int k = 0; k < dim; ++k) {
                        force[i * dim + k] -= (f_ij * vec_dist[k] / r_mod);
                        force[j * dim + k] += (f_ij * vec_dist[k] / r_mod);
                    }
                }
            }
            j = list[j]; // Move to the next particle in the neighboring cell
        }
    }
}

int pbc_cells(int icell) {
    // Apply periodic boundary conditions to cell indices
    if (icell == 0) return ncells;
    else if (icell == ncells + 1) return 1;
    else return icell;
}

int main() {
    // Data files
    ifstream fich_pos;
    fich_pos.open("position.txt");
    // Initial positions: prepared from monte carlo: read pos
    for (int i = 0; i < 3 * Npart; ++i)
        fich_pos >> Pos[i];
    fich_pos.close();
    make_linked_list(Pos[dim*Npart]); // Build the linked list
    // Perform interactions
    for (int i = 0; i < natoms; ++i) {
        for (int k=0; k<3; ++k) float postag[k]=Pos[i*dim+k];

         forceNL(i, postag[dim], Pos[Npart*dim]);
     }
    delete[] positions; // Free allocated memory
    return 0;
}




