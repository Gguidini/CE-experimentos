// g++ -std=c++17 -O0 -march=native -I ./ -I ~/include sux.cpp
#include <sux/bits/Rank9Sel.hpp>

#include <random>
#include <chrono>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

#define MiB 1048576
#define SAMPLE_SIZE 5000

#define SIZES {1, 16, 64, 128, 256, 512}
#define DENSITIES {0.05, 0.1, 0.25, 0.5, 0.75, 0.9}

std::random_device rand_dev;
std::mt19937 generator(rand_dev());

int64_t onesInBitVector;

using namespace sux::bits;
using namespace std;

// CSV FIle with following headers:
// Base tructure | Operation | Support Structure | Size of base | Size of structure | Time to create | ...statistics
ofstream outputFile;

// AUX functions, declared at the end of the file
int64_t* loadBitVector(string fileName, int64_t size);
vector<double> calculateStatistics(vector<int64_t>& observations);

void testRank(int64_t* b, int64_t sizeInBytes, std::uniform_int_distribution<int64_t>& distribution) {

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    Rank9Sel<> ranker((const uint64_t*) b, sizeInBytes);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    std::cout << "STRUCTURE - Rank9Sel | OPERATION - RANK\n";
    std::cout << "Time to create struture: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "µs" << endl;
    std::cout << "Size of structure in bytes: " << ranker.size() << endl;

    //Warm up structure with 1000 rank operations
    for(int i = 0; i < 1000; i++) ranker.rank(i);

    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    std::cout << "Statistics. Sample size: " << SAMPLE_SIZE << endl;
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = (idx + distribution(generator)) % sizeInBytes;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        ranker.rank(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }

    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    outputFile << "int64_t*,rank,Rank9Sel," << sizeInBytes << "," << ranker.size() << ",";
    outputFile << chrono::duration_cast<chrono::microseconds>(end - begin).count() << ",";
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

void testSelect(int64_t* b, int64_t sizeInBytes, std::uniform_int_distribution<int64_t>& distribution) {
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    Rank9Sel<> ranker((const uint64_t*) b, sizeInBytes);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    std::cout << "STRUCTURE - Rank9Sel | OPERATION - SELECT\n";
    std::cout << "Time to create struture: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "µs" << endl;
    std::cout << "Size of structure in bytes: " << ranker.size() << endl;

    //Warm up structure with 1000 rank operations
    for(int i = 0; i < 1000; i++) ranker.select((i % onesInBitVector));

    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    std::cout << "Statistics. Sample size: " << SAMPLE_SIZE << endl;
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = (idx + distribution(generator)) % onesInBitVector;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        ranker.select(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }

    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    outputFile << "int64_t*,select,Rank9Sel," << sizeInBytes << "," << ranker.size() << ",";
    outputFile << chrono::duration_cast<chrono::microseconds>(end - begin).count() << ",";
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

int main() {

    std::uniform_int_distribution<int64_t> distribution(0, 1024);

    for (int64_t vectorSize : SIZES) {
        for (double density : DENSITIES) {

            stringstream ss;
            ss << "bitvector_" << vectorSize << "MB_density" << (density * 100);
            string inFile = "data/" + ss.str() + ".txt";
            string outFile = "results/sux/" + ss.str() + ".csv";
            // Bitvector size
            size_t sizeInBytes = vectorSize * MiB;

            std::cout << "\nStarting test with " << vectorSize << "MB and density " << density * 100 << "%\n";
            std::cout << "Input file is " << inFile << endl;
            std::cout << "Output file is " << outFile << endl;
    
            outputFile.open(outFile);

            outputFile << "Base structure,Operation,Support Structure,Size of base (bytes),Size of structure (bytes),Time to create (µs),";
            outputFile << "média (ns),variância amostral,desvio padrão amostral, mediana (ns)\n";
            // Load bitvector from test file
            int64_t* b = loadBitVector(inFile, sizeInBytes);


            testRank(b, sizeInBytes, distribution);
            outputFile << "\n";
            testSelect(b, sizeInBytes, distribution);
            outputFile << "\n";
            
            outputFile.close();
            free(b);
        }
    }

}


// AUX Functions
int64_t* loadBitVector(string fileName, int64_t size) {
    onesInBitVector = 0;
    ifstream myfile (fileName);
    if(!myfile.is_open()) {
        cout << "Não conseguiu abrir o arquivo\n";
        exit(1);
    }
    int64_t* bitvector = (int64_t*) calloc(size, 8);  // Number of bytes
    int64_t addr = 0, offset = 0;
    while(myfile) {
        char next = myfile.get();
        if (next == '1') {
            bitvector[addr / 64] |= (1LL << offset);
            onesInBitVector++;
        } 
        addr++;
        offset = (offset + 1)%64;
    }
    return bitvector;
}

vector<double> calculateStatistics(vector<int64_t>& observations) {
    vector<double> r;
    double average = accumulate(observations.begin(), observations.end(), 0) / observations.size();
    r.push_back(average);
    double variance_sample = 0;
    for (int64_t obs : observations) {
        double diff = (obs - average);
        variance_sample += (diff * diff);
    }
    variance_sample /= (observations.size() - 1);
    double desvio_padrao_amostra = sqrt(variance_sample);
    r.push_back(variance_sample);
    r.push_back(desvio_padrao_amostra);

    sort(observations.begin(), observations.end());
    double mediana = observations.size()%2 == 0 
        ? observations[observations.size() / 2] 
        : (observations[observations.size() / 2] + observations[observations.size() / 2 + 1]) / 2;
    r.push_back(mediana);
    return r;
}