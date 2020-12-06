#include <random>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>


#define RUN_SIZE 1000
#define SIZE 64
#define MiB 1048576

std::random_device rand_dev;
std::mt19937 generator(rand_dev());


int main() {
    // Random number generator
    std::uniform_int_distribution<int64_t> distribution(0, 1024);
    int64_t sizeInBits = SIZE * MiB * 8;

    // Output file
    std::stringstream ss;
    ss << "utils/indexes_" << RUN_SIZE << ".txt";
    std::ofstream outputFile(ss.str());
    outputFile << RUN_SIZE << std::endl;
    // Não pode haver indice 0 por causa dos selects
    int64_t idx = distribution(generator);
    outputFile << idx;
    if(idx == 0) idx++;

    for(size_t i = 1; i < RUN_SIZE; i++) {
        idx = (idx + distribution(generator)) % sizeInBits;
        if(idx == 0) idx++;
        outputFile << " " << idx;
    }
    outputFile << "\n";
}