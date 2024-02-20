// LIBRARIES INCLUDED

// <iostream>: input and output operations
// <fstream>: file manipulation
// <sstream>: string data manipulation
// <cmath>: Mathematical functions and constants

#include <iostream> 
#include <fstream> 
#include <sstream>
#include <cmath>

// This allow to use declarations in "std" namespace without calling it
// Most of the variables or operations tha we declarate are in this namespace, so this line is convenient
using namespace std;

// MAIN FUNCTION
// The compiled executable can get some data as input
// "argc" gives the number or arguments that the function has accepted as input, including the exceuting command
// "argv" is an array (string type) that includes each of the arguments given as inputs
// For example, if we run ./main 5 3.0 as an executable: argc = 3, argv = {"./main", "5", "3.0"}
int main(int argc, char* argv[]) {

    // Create the positions vector from the file
	
    fstream archivo("positions/equilibrium_positions.txt");
    // Check if the file is open
    if (!archivo.is_open()) 
    {
	cerr << "Error: Unable to open file!" << endl;
        return 1;
    }

    vector<vector<float>> positions;	// A vector array (vector of vectors) with all the positions
	
    /* getline is an "std" function that creates the "line" string by reading "file" until a newline
    character or line ending appears, so all the reading characters are stored in "line" */
    while (getline(file, line)) 
    {
        vector<float> position; // A vector array with a certain position
        stringstream ss(line);  // The "ss" stream object is created
        float value;		// An auxiliar float variable is created
        while (ss >> value) {	// We read the line in sstream format
            // Push back the value into the row vector
            position.push_back(value);	// We store a certain value
            ss.ignore();		// Ignore the delimiter (in this case a space)
        }
        positions.push_back(position);	// Push back the "position" vector into the "positions" vector
    }
    // while (archivo >> energy) x.push_back(energy);
    archivo.close();	// We close the file

	
    // RDF construction
	
    int dim = position.size();		// We define the dimension from the position coordinates
    int Npart = positions.size()/dim;	// And here, the number of particles from the toltal array

    // We are not considering the last particle (it is cosidered in all steps before)
    if (i == Npart-1) continue; 
    for (int j = i+1; j < Npart; ++j) // We count the j particles for pairing
    {
	float ri[dim];	// We initialize the particle i position
	float rj[dim];	// We initialize the particle j position
	for (int k = 0; k<dim; ++k)
	{
	     	// Entering values from the complete list
		ri[k] = positions[i * dim + k]; 
		rj[k] = positions[j * dim + k];
	}
	dist = mic_distance(ri,rj);
	if (dist >= dmax) continue;	
    }
	
    // Variables for the histogram
    const int numBins = atoi(argv[1]);
    float xmax = *max_element(x.begin(), x.end());
    float xmin = *min_element(x.begin(), x.end());
    float xrange = xmax - xmin;
    float binWidth = xrange / numBins;
    vector<int> histogram(numBins, 0);

    // Create the histogram
    for (float value : x) {
        int binIndex = static_cast<int>((value - xmin) / binWidth);
        histogram[binIndex]++;
    }

    // Normalization of the histogram
    cout << "Histogram: " <<endl <<endl;
    vector<float> binNormalizedValues(numBins, 0.0);
    for (int i = 0; i < numBins; ++i) {
        binNormalizedValues[i] = static_cast<float>(histogram[i]) / (x.size() * binWidth);

        // Print the histogram
        cout << "   Bin " << i + 1 << "= " << binNormalizedValues[i] << endl <<endl;
    }

    // Calculate variance
    float mean = accumulate(binNormalizedValues.begin(), binNormalizedValues.end(), 0.0) / numBins;
    float variance = 0.0;
    for (int i = 0; i < numBins; ++i) {
        variance += (binNormalizedValues[i] - mean)*(binNormalizedValues[i] - mean);
    }
    variance /= (numBins - 1);

    // Calculate the standard desviation
    float std = sqrt(variance);
    
    // Print the data
    cout << "Variance= " << variance << endl;
    cout << "Std= " << std << endl <<endl;
    
    return 0;
}
