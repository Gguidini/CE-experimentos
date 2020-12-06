// COMPILE WITH g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib sdsl-lite.cpp -lsdsl -ldivsufsort -ldivsufsort64
#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <math.h>


using namespace std;
using namespace sdsl;

#define MiB 1048576
#define SAMPLE_SIZE 5000

#define SIZES {1, 16, 64, 128, 256, 512}
#define DENSITIES {0.05, 0.1, 0.25, 0.5, 0.75, 0.9}

bit_vector b;

std::random_device rand_dev;
std::mt19937 generator(rand_dev());

int64_t onesInBitVector;
// CSV FIle with following headers:
// Base tructure | Operation | Support Structure | Size of base | Size of structure | Time to create | ...statistics
ofstream outputFile;

// AUX functions, declared at the end of the file
bit_vector& loadBitVector(string fileName, int64_t size);
vector<double> calculateStatistics(vector<int64_t>& observations);

void testBitVector(bit_vector& b, int64_t size, std::uniform_int_distribution<int64_t>& distribution) {
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    rank_support_v<> rb(&b);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    std::cout << "STRUCTURE - BITVECTOR (MUTABLE) | OPERATION - RANK\n";
    std::cout << "Time to create struture: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "µs" << endl;
    std::cout << "Size of the base structure in bytes: " << size_in_bytes(b) << endl;
    std::cout << "Size of structure in bytes: " << size_in_bytes(rb) << endl;

    //Warm up structure with 1000 rank operations
    for(int i = 0; i < 1000; i++) rb.rank(i);

    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    std::cout << "Statistics. Sample size: " << SAMPLE_SIZE << endl;
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = (idx + distribution(generator)) % size;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rb.rank(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }

    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    outputFile << "bitvector,rank,rank_support_v," << size_in_bytes(b) << "," << size_in_bytes(rb) << ",";
    outputFile << "0," << chrono::duration_cast<chrono::microseconds>(end - begin).count() << ","; // First time is to create base structure, but bit_vector is already crated
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

void testBitVectorV5(bit_vector& b, int64_t size, std::uniform_int_distribution<int64_t>& distribution) {
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    rank_support_v5<> rb(&b);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    std::cout << "STRUCTURE - BITVECTOR (MUTABLE) | OPERATION - RANK\n";
    std::cout << "Time to create struture: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "µs" << endl;
    std::cout << "Size of the base structure in bytes: " << size_in_bytes(b) << endl;
    std::cout << "Size of structure in bytes: " << size_in_bytes(rb) << endl;

    //Warm up structure with 1000 rank operations
    for(int i = 0; i < 1000; i++) rb.rank(i);

    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    std::cout << "Statistics. Sample size: " << SAMPLE_SIZE << endl;
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = (idx + distribution(generator)) % size;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rb.rank(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }

    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    outputFile << "bitvector,rank,rank_support_v5," << size_in_bytes(b) << "," << size_in_bytes(rb) << ",";
    outputFile << "0," << chrono::duration_cast<chrono::microseconds>(end - begin).count() << ",";
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

void testBitVectorSelect(bit_vector& b, int64_t size, std::uniform_int_distribution<int64_t>& distribution) {
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    select_support_mcl<> selectb(&b);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    std::cout << "STRUCTURE - BITVECTOR (MUTABLE) | OPERATION - SELECT\n";
    std::cout << "Time to create struture: " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "µs" << endl;
    std::cout << "Size of the base structure in bytes: " << size_in_bytes(b) << endl;
    std::cout << "Size of structure in bytes: " << size_in_bytes(selectb) << endl;

    // Warm up structure with 1000 select operations
    // select(0) gives segfault.
    for(int i = 0; i < 1000; i++) selectb.select((i+1) % onesInBitVector);
    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    std::cout << "Statistics. Sample size: " << SAMPLE_SIZE << endl;
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = ((idx + distribution(generator) + 1) % onesInBitVector);
        if(idx == 0) idx++;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        selectb.select(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }

    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    
    outputFile << "bitvector,select,select_support_v," << size_in_bytes(b) << "," << size_in_bytes(selectb) << ",";
    outputFile << "0," << chrono::duration_cast<chrono::microseconds>(end - begin).count() << ",";
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

void testRRRVector(bit_vector& b, int64_t size, std::uniform_int_distribution<int64_t>& distribution) {
    // Immutable, compressed structure
    int64_t timeToCreateBase = 0, timeToCreateStructure = 0;
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    rrr_vector<> rrrb(b);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    timeToCreateBase = chrono::duration_cast<chrono::microseconds>(end - begin).count();
    begin = chrono::steady_clock::now();
    rrr_vector<>::rank_1_type rank_rrrb(&rrrb);
    end = chrono::steady_clock::now();
    timeToCreateStructure = chrono::duration_cast<chrono::microseconds>(end - begin).count();

    std::cout << "STRUCTURE - RRR_BITVECTOR (IMMUTABLE) | OPERATION - RANK\n";
    std::cout << "Time to create rank struture: " << timeToCreateStructure << "µs" << endl;
    std::cout << "Size of the base structure in bytes (after support): " << size_in_bytes(rrrb) << endl;
    std::cout << "Size of rank structure in bytes: " << size_in_bytes(rank_rrrb) << endl;

    for(int i = 0; i < 1000; i++) rank_rrrb.rank(i);

    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = (idx + distribution(generator)) % size;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rank_rrrb.rank(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }
    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    outputFile << "rrr_bitvector,rank,rrr_rank_support," << size_in_bytes(rrrb) << "," << size_in_bytes(rank_rrrb) << ",";
    outputFile << timeToCreateBase << "," << timeToCreateStructure << ",";
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

void testRRRVectorSelect(bit_vector& b, int64_t size, std::uniform_int_distribution<int64_t>& distribution) {
    // Immutable, compressed structure
    int64_t timeToCreateBase = 0, timeToCreateStructure = 0;
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    rrr_vector<> rrrb(b);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    timeToCreateBase = chrono::duration_cast<chrono::microseconds>(end - begin).count();
    begin = chrono::steady_clock::now();
    rrr_vector<>::select_1_type select_rrrb(&rrrb);
    end = chrono::steady_clock::now();
    timeToCreateStructure = chrono::duration_cast<chrono::microseconds>(end - begin).count();

    std::cout << "STRUCTURE - RRR_BITVECTOR (IMMUTABLE) | OPERATION - SELECT\n";
    std::cout << "Time to create select struture: " <<  timeToCreateStructure << "µs" << endl;
    std::cout << "Size of the base structure in bytes (after support): " << size_in_bytes(rrrb) << endl;
    std::cout << "Size of select structure in bytes: " << size_in_bytes(select_rrrb) << endl;

    for(int i = 0; i < 1000; i++) select_rrrb.select((i + 1) % onesInBitVector);

    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = ((idx + distribution(generator) + 1) % onesInBitVector);
        if(idx == 0) idx++;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        select_rrrb.select(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }
    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    outputFile << "rrr_bitvector,select,rrr_select_support," << size_in_bytes(rrrb) << "," << size_in_bytes(select_rrrb) << ",";
    outputFile << timeToCreateBase << "," << timeToCreateStructure << ",";
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

void testSDVector(bit_vector& b, int64_t size, std::uniform_int_distribution<int64_t>& distribution) {
    // Immutable, compressed structure
    int64_t timeToCreateBase = 0, timeToCreateStructure = 0;
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    sd_vector<> sdb(b);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    timeToCreateBase = chrono::duration_cast<chrono::microseconds>(end - begin).count();
    begin = chrono::steady_clock::now();
    sd_vector<>::rank_1_type rank_sdb(&sdb);
    end = chrono::steady_clock::now();
    timeToCreateStructure = chrono::duration_cast<chrono::microseconds>(end - begin).count();

    std::cout << "STRUCTURE - RRR_BITVECTOR (IMMUTABLE) | OPERATION - RANK\n";
    std::cout << "Time to create rank struture: " << timeToCreateStructure << "µs" << endl;
    std::cout << "Size of the base structure in bytes: " << size_in_bytes(sdb) << endl;
    std::cout << "Size of rank structure in bytes: " << size_in_bytes(rank_sdb) << endl;

    for(int i = 0; i < 1000; i++) rank_sdb.rank(i);

    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = (idx + distribution(generator)) % size;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rank_sdb.rank(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }
    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    outputFile << "sd_bitvector,rank,sd_rank_support," << size_in_bytes(sdb) << "," << size_in_bytes(rank_sdb) << ",";
    outputFile << timeToCreateBase << "," << timeToCreateStructure << ",";
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

void testSDVectorSelect(bit_vector& b, int64_t size, std::uniform_int_distribution<int64_t>& distribution) {
    // Immutable, compressed structure
    int64_t timeToCreateBase = 0, timeToCreateStructure = 0;
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    sd_vector<> sdb(b);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    timeToCreateBase = chrono::duration_cast<chrono::microseconds>(end - begin).count();
    begin = chrono::steady_clock::now();
    sd_vector<>::select_1_type select_sdb(&sdb);
    end = chrono::steady_clock::now();
    timeToCreateStructure = chrono::duration_cast<chrono::microseconds>(end - begin).count();

    std::cout << "STRUCTURE - RRR_BITVECTOR (IMMUTABLE) | OPERATION - SELECT\n";
    std::cout << "Time to create select struture: " << timeToCreateStructure << "µs" << endl;
    std::cout << "Size of the base structure in bytes: " << size_in_bytes(sdb) << endl;
    std::cout << "Size of select structure in bytes: " << size_in_bytes(select_sdb) << endl;

    for(int i = 0; i < 1000; i++) select_sdb.select(((i+1) % onesInBitVector));

    int64_t idx = 0;
    vector<int64_t> observations(SAMPLE_SIZE);
    // Collect sample observations
    for(int i = 0; i < SAMPLE_SIZE; i++) {
        idx = ((idx + distribution(generator) + 1) % onesInBitVector);
        if(idx == 0) idx++;
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        select_sdb.select(idx);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        observations.push_back(chrono::duration_cast<chrono::nanoseconds>(end - begin).count());
    }
    vector<double> statistics = calculateStatistics(observations);
    std::cout << "Average: " << statistics[0] << endl;
    std::cout << "Variancia Amostral: " << statistics[1] << endl;
    std::cout << "Desvio Padrão: " << statistics[2] << endl;
    std::cout << "Mediana: " << statistics[3] << endl;
    std::cout << "========================================================\n";

    outputFile << "sd_bitvector,select,sd_select_support," << size_in_bytes(sdb) << "," << size_in_bytes(select_sdb) << ",";
    outputFile << timeToCreateBase << "," << timeToCreateStructure << ",";
    outputFile << statistics[0] << "," << statistics[1] << "," << statistics[2] << "," << statistics[3];
}

int main() {

    // Random number generator
    std::uniform_int_distribution<int64_t> distribution(0, 1024);

    for (int64_t vectorSize : SIZES) {
        for (double density : DENSITIES) {
            stringstream ss;
            ss << "bitvector_" << vectorSize << "MB_density" << (density * 100);
            string inFile = "data/" + ss.str() + ".txt";
            string outFile = "results/sdsl/" + ss.str() + ".csv";
            // Bitvector size
            int64_t size = vectorSize * MiB * 8; // size in bits

            std::cout << "\nStarting test with " << vectorSize << "MB and density " << density * 100 << "%\n";
            std::cout << "Input file is " << inFile << endl;
            std::cout << "Output file is " << outFile << endl;
            // Load bitvector from test file
            bit_vector b = loadBitVector(inFile, size);
            outputFile.open(outFile);

            outputFile << "Base structure,Operation,Support Structure,Size of base (bytes),Size of structure (bytes),Time to create base structure (µs),Time to create (µs),";
            outputFile << "média (ns),variância amostral,desvio padrão amostral, mediana (ns)\n";
            
            testBitVector(b, size, distribution);
            outputFile << "\n";
            testBitVectorV5(b, size, distribution);
            outputFile << "\n";
            testBitVectorSelect(b, size, distribution);
            outputFile << "\n";
            testRRRVector(b, size, distribution);
            outputFile << "\n";
            testRRRVectorSelect(b, size, distribution);
            outputFile << "\n";
            testSDVector(b, size, distribution);
            outputFile << "\n";
            testSDVectorSelect(b, size, distribution);
            outputFile << "\n";

            outputFile.close();
            
        }
    }
}



// AUX Functions
bit_vector& loadBitVector(string fileName, int64_t size) {
    onesInBitVector = 0;
    ifstream myfile (fileName);
    if(!myfile.is_open()) {
        cout << "Não conseguiu abrir o arquivo\n";
        exit(1);
    }
    b.resize(size + 1);            // bit_vector(size, prefil)
    size_t i = 0;
    while(myfile) {
        char next = myfile.get();
        if (next == '1') {
            b[i] = 1;
            onesInBitVector++;
        } else {
            b[i] = 0;
        }
        i++;
    }
    cout << "Loaded bit vector\n";
    return b;
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
