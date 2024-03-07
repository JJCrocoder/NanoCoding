









int main(int argc)







void moveVerlet(int options, float Pos[], float Vel[], float For [])
{
  switch(options)
  {
    case 1:
    for (int i = 0; i < Npart; ++i)
    {
      for (int k = 0; k < dim; ++k)
      {
        Pos[i*dim + k] += Vel[i*dim + k]*deltaT + For[i*dim + k]*deltaT*deltaT*0.5;
        Vel[i*dim + k] += For[i*dim + k]*deltaT;
      }
    }
    break;
    case 2:
    for (int i = 0; i < Npart; ++i)
    {
      for (int k = 0; k < dim; ++k)
      {
        Vel[i*dim + k] += For[i*dim + k]*deltaT;
      }
    }
    break;
  }
}

void forces(float Pos[], float For[])
{
  float r2_cut = r_cut*r_cut;
  for (int i = 0; i < 3*Npart; ++i)
  {
    For[i]=0.0;
  }
  for (int i = 0; i<Npart-1; ++i)
  {
    for (int j = 0; j < Npart; ++j)
    {
      float vec_dist[dim];
      float r2_ij = 0.0;
        
      for (int k = 0; k < dim; ++k) vec_dist[k] = Pos[j*dim+k] - Pos[i * dim + k];
mic(vec_dist, Lbox);
      for (int k = 0; k < dim; ++k) r2_ij += vec_dist[k]*vec_dist[k];

      if (r2_ij < r_2cut)
      {
        float r_mod = sqrt(r2_ij);
        float r6 = pow((sigma/r_mod), 6.0);
        float f_ij = 48*epsilon*(r6*r6-0.5*r6)/r_mod
        for (int k = 0; k < dim; ++k)
        {
          for[j*dim + k] += (f_ij*vec_dist[k])/r_mod;
          for[i*dim + k] -= (f_ij*vec_dist[k])/r_mod;
        }
      }
    }
  }
}
