#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

/*WARNING
THIS CODE WILL WORK WITH THE TXT FILE OF THE EQUILIBRIUM, BEWARE OF THAT. IF U USE THE FULL ENERGY TXT FILE THE MEAN IS GONNA BE FUCKED UP AND SO THE HISTOGRAM COMPUTING*/
int main(int argc, char* argv[]) {
  
  //First I compute the mean and the variance 
  // Open the file
  std::ifstream file(argv[1]); //opens the file, and will not close it until the end cause it raises some errors
  if (file.is_open()){ // check that it is open
    std::string line; //initializes each line 
    double total; //I define the following variables to compute the mean
    int i = 0;
    double value;
    while (std::getline(file,line)){ //I iterate over each line 
      value = std::stof(line); // assign the variable value to each line
      total = total + value; // And sum all lines
      i++;
    }
    float mean=total/((float)i); //This is to compute the mean, as the name shows

    std::cout << "Mean" << mean <<std::endl; //this prints sout the result for the mean 
  
  //Now with the histogram
  //Setting the variables
  const int Nbin=std::stoi(argv[2]); //This way I define the Nbin by the screen. The const keeps its value constant (I don't know why but it is needed for the code to work)
  int histo[Nbin]={0}; //Initialize the histo array, {0} assign 0 to each value of the array
  float pdf[Nbin]={0.0};  //The same for the normalize values of the histogram
  int ncount=0; // count how many energies are we looking at 
  float xmin; //Initialize the cutoffs for the energy
  float xmax;

  if (mean < 0.0){ //Here we compute the values of the xmin and xmax. I decided this method just because I couldn't think of any other one (Maybe using the variance ??). How xmax/min is computed depends on if mean > 0
   xmin= mean + mean/2.0;
   xmax= mean - mean/2.0;
  }else{
   xmin= mean - mean/2.0;
   xmax= mean + mean/2.0;
  }
  std::cout << "Xmax " << xmax << std::endl;
  std::cout << "Xmin " << xmin << std::endl;
  float xrange=xmax-xmin; //Compute the range and the step
  float deltax = xrange/(float)Nbin;

//Histogram creation

  std::ifstream file(argv[1]); //Again, read the file
    // Read each line of the file
    while (std::getline(file, line)) {
      // Process the line (e.g., print it)
        double value = std::stof(line); //As before, we read each line and store the valeu 
        ncount+=1; //We add one to the count even if the value is larger than the range so the pdf doesn't sum up to 1
        if (value < xmin ) continue; //If the value is within our range 
        if (xmax < value ) continue;
        //We compute to which bin does the value belongs
        int ibin = std::floor((value-xmin)/deltax);  
        histo[ibin]+=1; //We add one to the bin the number belongs
  }

    for (int i = 0 ; i< Nbin; ++i){
        pdf[i]=(float)(histo[i])/(deltax*(float)ncount); //Here we compute the pdf based on the formula Rafa put in the whiteboard
    }

//The PDF has to be equal 1 or be less than one when we multiply it by the deltax
    
    // Close the file
    file.close(); //close the file

    std::ofstream fich_rdf; //We create a new file called fich_rdf

    fich_rdf.open("rdf.txt"); //We open it
    for (int i=0; i < Nbin ; ++i){
        fich_rdf << pdf[i] << std::endl ;  //We store pdf value inside
    }
    fich_rdf.close();
      } else {
    std::cerr << "Error opening file!" << std::endl; //This is just to print an error if the file doesn't exists
  }
    return 0;
}
