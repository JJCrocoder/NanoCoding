#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

/* WARNING
THIS CODE WILL WORK WITH THE EQUILIBRIUM TXT FILE, BEWARE OF THAT. 
IF U USE THE FULL ENERGY TXT FILE THE MEAN IS GONNA BE FUCKED UP AND SO THE HISTOGRAM COMPUTING */
int main(int argc, char* argv[]) { 
  
  // First I compute the mean and the variance 
  
  // Open the file
  std::ifstream file(argv[1]); //opens the file, and will not close it until the end cause it raises some errors
  if (file.is_open()) // check that it is open
  {    
    std::string line;     // Initializes each line 
    double total;         // I define the following variables to compute the mean
    int i = 0;            // A counter for the amount of data
    double value;         // An auxiliar variable that stores the double value in the iterated line 
    while (std::getline(file,line)) // I iterate over each line
    {
      value = std::stof(line); // Assign the variable value to each line
      total = total + value;   // And sum all lines
      i++;                     // Counter advance
    }
    float mean = total/((float)i); // This is to compute the mean from the cumsum

    std::cout << "Mean" << mean <<std::endl; // this prints the result for the mean 
  
    // HISTOGRAM
    
    // Setting the parameters and range of representation
    
    // Nbin is defined by scrren input
    const int Nbin=std::stoi(argv[2]);  // Const keeps the value constant (I don't know why but it is needed for the code to work)
    int histo[Nbin]={0};     // Initialize the histo array, {0} assign 0 to each value of the array
    float pdf[Nbin]={0.0};   // Same for the normalized values of the histogram
    int ncount=0;            // Count how many energies are we looking at
    // Initialize the limits for histogram representation
    float xmin; float xmax;
  
    // Here we compute the values of the xmin and xmax. 
    // I decided this method just because I couldn't think of any other one (Maybe using the variance ??). How xmax/min is computed depends on if mean > 0
    xmin = mean - abs(mean)/2.0;
    xmax = mean + abs(mean)/2.0;

    // We print the values on the screen
    std::cout << "Xmax " << xmax << std::endl;
    std::cout << "Xmin " << xmin << std::endl;
    // Compute the range and the step
    float xrange=xmax-xmin; float deltax = xrange/(float)Nbin;
  
    // Histogram creation
  
    std::ifstream file(argv[1]); //Again, read the file
    // Read each line of the file
    while (std::getline(file, line)) 
    {
      double value = std::stof(line);
      ncount+=1; // We add one to the count even if the value is larger than the range so the pdf doesn't sum up to 1
      if (value < xmin ) continue; // If the value is within our range 
      if (xmax < value ) continue;
      //We compute to which bin does the value belongs
      int ibin = std::floor((value-xmin)/deltax);  
      histo[ibin]+=1; // We add one to the bin the number belongs to
    }
  
    for (int i = 0 ; i< Nbin; ++i)
    {
      //Here we compute the pdf based on the formula Rafa put in the whiteboard  
      pdf[i]=(float)(histo[i])/(deltax*(float)ncount);
    }
    // Close the file
    file.close(); // Close the file

    std::ofstream fich_hist;     // We create an output file called fich_hist
  
    fich_hist.open("hist.txt");  // We open it
    
    for (int i=0; i < Nbin ; ++i)
    {
      fich_hist << pdf[i] << std::endl ;  //We store pdf value inside
    }
  
    fich_hist.close();
  } 

  // This is just to print an error if the file doesn't exists  
  else
  {
    std::cerr << "Error opening file!" << std::endl;
  }
  
  return 0;
}
