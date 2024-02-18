// Libraries included
#include <iostream>   // input and output operations
#include <fstream>    // file manipulation
#include <vector>     // 
#include <algorithm>
#include <cmath>
#include <numeric>
using namespace std;

// Main function
int main(int argc, char* argv[]) {

    // Create the energies vector from the file
    fstream archivo("energies/equilibrium_energies.txt");
    vector<float> energies;
    float energy;
    vector<float> x;
    while (archivo >> energy) x.push_back(energy);
    archivo.close();

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
