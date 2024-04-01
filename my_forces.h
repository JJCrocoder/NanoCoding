float lj_force (float r2){
  if (r2 > r_cut*r_cut) return 0.0;
  else {
    float r6 = pow((sigma*sigma/r2), 3.0);
    return 24 * epsilon * ((2*r6*r6-r6))/sqrt(r2);  
  }
}
