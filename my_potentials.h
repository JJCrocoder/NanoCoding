

float lj_pot (float r2){
  if (r2 > r_cut*r_cut) return 0;
  else {
    float r6 = pow((sigma*sigma/r2), 3.0);
    float rc6 = pow((sigma/r_cut), 6.0);
    return 4 * epsilon * ((r6*r6-r6) - (rc6*rc6-rc6));  
  }
}
