// COMPILE WITH g++ -std=c++17 -O3 -DNDEBUG -I ./ -I ~/include -L ~/lib distributionTest.cpp -lsdsl -ldivsufsort -ldivsufsort64 -march=native
#include <sdsl/bit_vectors.hpp>
#include <sux/bits/Rank9Sel.hpp>

#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>


using namespace std;
using namespace sdsl;
using namespace sux::bits;

#define MiB 1048576
#define SAMPLE_SIZE 500
#define RUN_SIZE 1000

#define SIZE 64
#define DENSITIES {0.05, 0.5, 0.9}

bit_vector b;

int64_t onesInBitVector;
// CSV FIle with following headers:
// Base tructure | Operation | Support Structure | Size of base | Size of structure | Time to create | ...statistics
ofstream outputFile;

// AUX functions, declared at the end of the file
bit_vector& loadBitVector(string fileName, int64_t size);
int64_t* loadBitVectorSux(string fileName, int64_t size);
vector<pair<string,double>> calculateStatistics(vector<double>& observations);
void saveInfoToFile(string baseStructure, string testStructure, bool isRankOp, vector<pair<string,double>>& statistics, vector<double>& sample);

void testBitVector(bit_vector& b, vector<int64_t>& indexes) {
    std::cout << "bit_vector - rank_support_v - rank\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < SAMPLE_SIZE; test++) {
        std::cout << test+1 << "/500\r";
        // Result = (time to create structure + time of SAMPLE_SIZE operations)
        // Create structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rank_support_v<> rb(&b);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 0; i < 1000; i++) rb.rank(i);
        }

        for(int64_t idx: indexes) {
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            rb.rank(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime = totalTime + (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }
        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("bit_vector", "rank_support_v", true, stats, sample);
}

void testBitVectorV5(bit_vector& b, vector<int64_t>& indexes) {
    std::cout << "bit_vector - rank_support_v5 - rank\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < 500; test++) {
        std::cout << test+1 << "/500\r";
        // Result = (time to create structure + time of SAMPLE_SIZE operations)
        // Create structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rank_support_v5<> rb(&b);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 0; i < 1000; i++) rb.rank(i);
        }

        for(int64_t idx: indexes) {
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            rb.rank(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime += (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }

        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("bit_vector", "rank_support_v5", true, stats, sample);
}

void testBitVectorSelect(bit_vector& b, vector<int64_t>& indexes) {
    std::cout << "bit_vector - select_support_mcl - select\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < 500; test++) {
        std::cout << test+1 << "/500\r";
        
        // Result = (time to create structure + time of SAMPLE_SIZE operations)
        // Create structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        select_support_mcl<> selectb(&b);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 1; i < 1000; i++) selectb.select(i);
        }

        for(size_t i = 0; i < RUN_SIZE; i++) {
            int64_t idx = indexes[i];
            if (idx > onesInBitVector) idx = onesInBitVector;
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            selectb.select(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime += (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }

        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("bit_vector", "select_support_mcl", false, stats, sample);
}

void testRRRVector(bit_vector& b, vector<int64_t>& indexes) {
    std::cout << "rrr_bit_vector - rank_support_rrr - rank\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < 500; test++) {
        std::cout << test+1 << "/500\r";
        // Result = (time to create structure + time of SAMPLE_SIZE operations)
        // Create base structure and rank structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rrr_vector<> rrrb(b);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();
        begin = chrono::steady_clock::now();
        rrr_vector<>::rank_1_type rank_rrrb(&rrrb);
        end = chrono::steady_clock::now();
        totalTime += chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 0; i < 1000; i++) rank_rrrb.rank(i);
        }

        for(int64_t idx: indexes) {
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            rank_rrrb.rank(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime += (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }

        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("rrr_bit_vector", "rank_support_rrr", true, stats, sample);
}

void testRRRVectorSelect(bit_vector& b, vector<int64_t>& indexes) {
    std::cout << "rrr_bit_vector - select_support_rrr - select\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < 500; test++) {
        std::cout << test+1 << "/500\r";
        // Create base structure and select structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        rrr_vector<> rrrb(b);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();
        begin = chrono::steady_clock::now();
        rrr_vector<>::select_1_type select_rrrb(&rrrb);
        end = chrono::steady_clock::now();
        totalTime += chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 1; i < 1000; i++) select_rrrb.select(i);
        }

        for(size_t i = 0; i < RUN_SIZE; i++) {
            int64_t idx = indexes[i];
            if (idx > onesInBitVector) idx = onesInBitVector;
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            select_rrrb.select(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime += (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }

        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("rr_bit_vector", "select_support_rrr", false, stats, sample);
}

void testSDVector(bit_vector& b, vector<int64_t>& indexes) {
    std::cout << "sd_bit_vector - rank_support_sd - rank\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < 500; test++) {
        std::cout << test+1 << "/500\r";
        // Result = (time to create structure + time of SAMPLE_SIZE operations)
        // Create base structure and rank structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        sd_vector<> sdb(b);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();
        begin = chrono::steady_clock::now();
        sd_vector<>::rank_1_type rank_sdb(&sdb);
        end = chrono::steady_clock::now();
        totalTime += chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 0; i < 1000; i++) rank_sdb.rank(i);
        }

        for(int64_t idx: indexes) {
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            rank_sdb.rank(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime += (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }

        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("sd_bit_vector", "rank_support_sd", true, stats, sample);
}

void testSDVectorSelect(bit_vector& b, vector<int64_t>& indexes) {
    std::cout << "sd_bit_vector - rank_support_sd - select\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < 500; test++) {
        std::cout << test+1 << "/500\r";
        // Result = (time to create structure + time of SAMPLE_SIZE operations)
        // Create base structure and rank structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        sd_vector<> sdb(b);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();
        begin = chrono::steady_clock::now();
        sd_vector<>::select_1_type select_sdb(&sdb);
        end = chrono::steady_clock::now();
        totalTime += chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 1; i < 1000; i++) select_sdb.select(i);
        }

        for(size_t i = 0; i < RUN_SIZE; i++) {
            int64_t idx = indexes[i];
            if (idx > onesInBitVector) idx = onesInBitVector;
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            select_sdb.select(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime += (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }

        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("sd_bit_vector", "select_support_sd", false, stats, sample);
}

void testRank9(int64_t* b, int64_t sizeInBytes, vector<int64_t>& indexes) {
    std::cout << "int64_t* - Rank9Sel - rank\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < SAMPLE_SIZE; test++) {
        std::cout << test+1 << "/500\r";
        // Result = (time to create structure + time of SAMPLE_SIZE operations)
        // Create structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        sux::bits::Rank9Sel<> ranker((const uint64_t*) b, sizeInBytes);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 0; i < 1000; i++) ranker.rank(i);
        }

        for(int64_t idx: indexes) {
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            ranker.rank(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime = totalTime + (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }
        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("int64_t*", "Rank9Sel", true, stats, sample);
}

void testSelect9(int64_t* b, int64_t sizeInBytes, vector<int64_t>& indexes) {
    std::cout << "int64_t* - Rank9Sel - select\n";
    // Run 500 tests, save results.
    vector<double> sample(SAMPLE_SIZE);
    for (int test = 0; test < SAMPLE_SIZE; test++) {
        std::cout << test+1 << "/500\r";
        // Result = (time to create structure + time of SAMPLE_SIZE operations)
        // Create structure
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        sux::bits::Rank9Sel<> selector((const uint64_t*) b, sizeInBytes);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        double totalTime = chrono::duration_cast<chrono::microseconds>(end - begin).count();

        if(test == 0) {
            // First run gets warm up round, untimed
            for(int i = 0; i < 1000; i++) selector.select(i);
        }

        for(size_t i = 0; i < RUN_SIZE; i++) {
            int64_t idx = indexes[i];
            if (idx > onesInBitVector) idx = onesInBitVector;
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            selector.select(idx);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            totalTime = totalTime + (chrono::duration_cast<chrono::nanoseconds>(end - begin).count() / 1000.0);
        }
        sample[test] = totalTime;
    }
    std::cout << "test ended ✔\n";
    vector<pair<string,double>> stats = calculateStatistics(sample);
    saveInfoToFile("int64_t*", "Rank9Sel", false, stats, sample);
}

void runSdslTests(string inFile, int64_t sizeInBits, vector<int64_t>& index) {
    // Load bitvector from test file
    std::cout << "==> SDSL Structures ▶\n";
    bit_vector b = loadBitVector(inFile, sizeInBits);
    // Run tests per structure
    testBitVector(b, index);
    outputFile << ",\n";
    testBitVectorV5(b, index);
    outputFile << ",\n";
    // testBitVectorSelect(b, index);
    // outputFile << ",\n";
    testRRRVector(b, index);
    outputFile << ",\n";
    testRRRVectorSelect(b, index);
    outputFile << ",\n";
    testSDVector(b, index);
    outputFile << ",\n";
    testSDVectorSelect(b, index);
    std::cout << "==> All SDSL Structures tested ✔\n";
}

void runSuxTests(string inFile, int64_t sizeInBytes, vector<int64_t>& index) {
    // Load bit vector
    std::cout << "==> Sux Structures ▶\n";
    int64_t* b = loadBitVectorSux(inFile, sizeInBytes);
    // Run tests per structure
    testRank9(b, sizeInBytes, index);
    outputFile << ",\n";
    testSelect9(b, sizeInBytes, index);
    std::cout << "==> All Sux Structures tested ✔\n";
}

int main() {

    // Load index from file
    stringstream indexFile;
    indexFile << "utils/indexes_" << RUN_SIZE << ".txt";
    
    std::cout << "Carregando indices do arquivo " << indexFile.str() << " (RUN_SIZE é " << RUN_SIZE << ")\n";

    ifstream isIndex(indexFile.str());
    vector<int64_t> index(RUN_SIZE);
    int64_t nextIndex;
    isIndex >> nextIndex;
    // First we check the size
    if(nextIndex != RUN_SIZE) {
        std::cerr << "Tamanho do index file ("<< nextIndex <<") e RUN_SIZE ("<< RUN_SIZE << ") são diferentes. Abortando\n";
        exit(1);
    }
    for(size_t i = 0; i < RUN_SIZE; i++) {
        isIndex >> index[i];
    }
    std::cout << "Index carregado ✔\n";

    stringstream outName;
    outName << "statistics/" << SIZE << "MB.json";
    outputFile.open(outName.str());
    outputFile << "{\n";

    int64_t sizeInBits = SIZE * MiB * 8;
    for (double density : DENSITIES) {
        stringstream ss;
        ss << "bitvector_" << SIZE << "MB_density" << (density * 100);
        string inFile = "data/" + ss.str() + ".txt";
        
        std::cout << "\nStarting test with " << SIZE << "MB and density " << density * 100 << "%\n";
        outputFile << "\t\"density_" << (density * 100) << "\": [\n";
        // Run sdsl tests
        runSdslTests(inFile, sizeInBits, index);
        outputFile << ",\n";
        runSuxTests(inFile, sizeInBits >> 3, index);
        outputFile << "\n";

        outputFile << "\t],\n";
    }
    outputFile << "}\n";
    outputFile.close();
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
    cout << "Loaded bit vector ✔\n";
    return b;
}

int64_t* loadBitVectorSux(string fileName, int64_t size) {
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
        } 
        addr++;
        offset = (offset + 1)%64;
    }
    cout << "Loaded bit vector ✔\n";
    return bitvector;
}

vector<pair<string,double>> calculateStatistics(vector<double>& observations) {
    vector<pair<string,double>> r;
    double average = accumulate(observations.begin(), observations.end(), 0) / observations.size();
    r.push_back(make_pair("media", average));
    double variance_sample = 0;
    for (int64_t obs : observations) {
        double diff = (obs - average);
        variance_sample += (diff * diff);
    }
    variance_sample /= (observations.size() - 1);
    double desvio_padrao_amostra = sqrt(variance_sample);
    r.push_back(make_pair("variancia_amostral", variance_sample));
    r.push_back(make_pair("desvio_padrao_amostral", desvio_padrao_amostra));

    sort(observations.begin(), observations.end());
    double mediana = observations.size()%2 == 0 
        ? observations[observations.size() / 2] 
        : (observations[observations.size() / 2] + observations[observations.size() / 2 + 1]) / 2;
    r.push_back(make_pair("mediana", mediana));
    return r;
}

string vectorToString(vector<double>& array) {
    if(array.size() == 0) {
        return "[]";
    }
    stringstream ss;
    ss << "[" << array[0];
    for (size_t i = 1; i < array.size(); i++) {
        ss << ", " << array[i];
    }
    ss << "]";
    return ss.str();
}

string vectorToString(vector<pair<string, double>>& array) {
    if(array.size() == 0) {
        return "[]";
    }
    stringstream ss;
    ss << "[[\"" << array[0].first << "\", " << array[0].second << "]";
    for (size_t i = 1; i < array.size(); i++) {
        ss << ", [\"" << array[i].first << "\", " << array[i].second << "]";
    }
    ss << "]";
    return ss.str();
}

void saveInfoToFile(string baseStructure, string testStructure, bool isRankOp, vector<pair<string,double>>& statistics, vector<double>& sample) {
    outputFile << "\t\t{\n";
    outputFile << "\t\t\t" << "\"base_structure\": \"" << baseStructure << "\",\n"; 
    outputFile << "\t\t\t" << "\"test_structure\": \"" << testStructure << "\",\n"; 
    outputFile << "\t\t\t" << "\"run_size\": " << RUN_SIZE << ",\n"; 
    outputFile << "\t\t\t" << "\"operation\": " << (isRankOp ? "\"rank\"" : "\"select\"") << ",\n"; 
    outputFile << "\t\t\t" << "\"statistics\": " << vectorToString(statistics) << ",\n";
    outputFile << "\t\t\t" << "\"sample\": " << vectorToString(sample) << "\n";
    outputFile << "\t\t}";
}