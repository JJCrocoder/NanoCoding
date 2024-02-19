#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

using namespace std;

int main() {
    // Abrir el archivo energy.txt para lectura
    ifstream infile("energy.txt");
    if (!infile.is_open()) { // Verificar si el archivo se abrió correctamente
        cerr << "Error: No se pudo abrir el archivo energy.txt" << endl;
        return 1;
    }

    // Leer los valores de energía del archivo y almacenar los últimos 300
    vector<float> energies; // Vector para almacenar los valores de energía
    float energy;
    while (infile >> energy) { // Leer cada valor de energía del archivo y agregarlo al vector
        energies.push_back(energy);
    }
    infile.close(); // Cerrar el archivo después de leer

    int numValues = energies.size(); // Obtener el número total de valores de energía
    int numToConsider = min(numValues, 300); // Considerar los últimos 300 valores o menos si hay menos de 300 en total
    vector<float> last300Values(energies.end() - numToConsider, energies.end()); // Obtener los últimos 300 valores

    // Construir el histograma
    const int numBins = 50; // Número de bins del histograma
    float minVal = *min_element(last300Values.begin(), last300Values.end()); // Obtener el valor mínimo de los últimos 300 valores
    float maxVal = *max_element(last300Values.begin(), last300Values.end()); // Obtener el valor máximo de los últimos 300 valores
    float binWidth = (maxVal - minVal) / numBins; // Calcular el ancho de cada bin

    vector<int> histogram(numBins, 0); // Inicializar el histograma con ceros

    // Llenar el histograma contando las ocurrencias en cada bin
   for (float val : last300Values) {
    int binIndex = static_cast<int>((val - minVal) / binWidth); // Calcular el índice del bin correspondiente para el valor
    // Asegurarse de que el índice del bin esté dentro de los límites del histograma
    binIndex = min(max(binIndex, 0), numBins - 1);
    if (val == maxVal) {
        // Si el valor es igual al valor máximo, se incluye en el último bin del histograma
        binIndex = numBins - 1;
    }
    histogram[binIndex]++; // Incrementar el conteo del bin correspondiente
}

    // Normalizar el histograma dividiendo cada conteo por la suma total de conteos
    int totalOccurences = accumulate(histogram.begin(), histogram.end(), 0); // Calcular el número total de ocurrencias
    vector<float> normalizedHistogram(numBins);
    transform(histogram.begin(), histogram.end(), normalizedHistogram.begin(),
              [totalOccurences](int count) { return static_cast<float>(count) / totalOccurences; });

    // Abrir un archivo para escribir el histograma normalizado
    ofstream outfile("normalized_histogram.txt");
    if (!outfile.is_open()) { // Verificar si el archivo se abrió correctamente
        cerr << "Error: No se pudo abrir el archivo normalized_histogram.txt para escritura." << endl;
        return 1;
    }

    // Escribir el histograma normalizado en el archivo
    for (int i = 0; i < numBins; ++i) {
        float binStart = minVal + i * binWidth; // Obtener el límite inferior del bin
        float binEnd = binStart + binWidth; // Obtener el límite superior del bin
        outfile << binStart << "\t" << binEnd << "\t" << normalizedHistogram[i] << endl; // Escribir el rango del bin y su valor normalizado
    }

    outfile.close(); // Cerrar el archivo después de escribir

    cout << "Histograma normalizado guardado en normalized_histogram.txt." << endl; // Mensaje de confirmación

     // Calcular la media de los últimos 300 valores de energía
    float mean = accumulate(last300Values.begin(), last300Values.end(), 0.0) / last300Values.size();

    // Calcular la varianza de los últimos 300 valores de energía
    float variance = 0.0;
    for (float val : last300Values) {
        variance += pow(val - mean, 2); // Calcular la diferencia al cuadrado entre cada valor y la media y sumarlas
    }
    variance /= last300Values.size() - 1; // Usar (N-1) para calcular la varianza

    // Calcular la desviación estándar como la raíz cuadrada de la varianza
    float stdDeviation = sqrt(variance);

    // Imprimir la desviación estándar y la varianza de los últimos 300 valores de energía
    cout << "Standard Deviation of last 300 values: " << stdDeviation << endl;
    cout << "Variance of last 300 values: " << variance << endl;

    return 0;
}
