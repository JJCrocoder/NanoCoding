// Minimum image convention (periodic boundary conditions)
// Void functons doesn't return anything but you can change an argument inside them
// In this case, vec is created outside the function, but modified inside of it
void mic(float * vec, float Lbox) {
    for (int i = 0; i < 3 ; ++i) {
        vec[i] -= floorf(0.5 + vec[i]/Lbox)*Lbox; // we apply the mic for each component
    }  
    return;
}

// We use the mic to calculate the image distance between 2 vectors
float mic_distance(float * Position_part_i, float *Position_part_j, float Lbox) {
    float rij[dim];
    float mod2_rij = 0.0; 
    for (int k = 0; k < 3; ++k) rij[k] = Position_part_j[k] - Position_part_i[k];
    mic(rij, Lbox);
    for (int k = 0; k < 3; ++k) mod2_rij += rij[k]*rij[k];
    return sqrt(mod2_rij);
}
