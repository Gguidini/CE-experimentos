#include <cstdint>
#include <math.h>
#include <random>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <bitset>
#include <sstream>
#include <string>
#include <fstream>

// Fatores controláveis do experimento
const int64_t bitVectorSize[11] = {1, 16, 64, 128, 256, 512, 1024}; // Tamanho em MB do bitVector
const double densityPercentage[6] = {0.05, 0.1, 0.25, 0.5, 0.75, 0.9}; // ( Num of 1s / length)

// Métricas
// tempo por operação
// tamanho da estrutura (em memória)
// tempo da criação da estrutura

// classes em estudo
// SDSL: M bit_vector, I rrr_vector, I sd_vector
// Sux: M StrideDynRankSel, I Rank9Sel  


void* createTestVector(int64_t size, double density) {
    int64_t bits = size * 1048576 * 8;  // size * 2**20 * 8 --> size in Mi bits
    // int64_t bits = size * 8;
    int64_t ones = floor(bits * density);
    // std::cout << "Gonna have " << ones << " ones" << std::endl;
    int8_t* array = (int8_t*) malloc( bits >> 3); // Reserva a quantidade de bytes
    if (density <= 0.5) {
        memset(array, 0, bits >> 3);
    } else {
        memset(array, 255, bits >> 3);
        ones = bits - ones;
    }
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<int64_t> distribution(0,bits - 1);
    while (ones) {
        int64_t idx = distribution(generator); // idx em bits
        // std::cout << "[" << ones << "] idx:" << idx << std::endl;
        int64_t address = idx / 8, offset = idx % 8;
        // std::cout << "address: " << address << " offset: " << offset << std::endl;
        if (density <= 0.5) {
            if (array[address] & (1 << offset)) 
                continue;
            array[address] |= (1 << offset);
        } else {
            if (!array[address] & (1 << offset)) 
                continue;
            array[address] &= ~(1 << offset);
        }
        ones--;
    }
    
    return array;
}

int main() {
    // Generate files for all vectorSizes and densities
    for (int64_t vectorSize : bitVectorSize) {
        for (double density : densityPercentage) {
            std::stringstream ss;
            ss << "data/bitvector_" << vectorSize << "MB_density" << (density * 100) << ".txt";
            std::string fileName = ss.str();
            std::cout << "Creating file " << fileName << std::endl;
            std::ofstream myfile;
            myfile.open (fileName);

            size_t sizeInMB = vectorSize * 1048576;
            int64_t* a = (int64_t*) createTestVector(vectorSize, density);

            for (size_t i = 0; (i*64) < sizeInMB; i++) {
                std::bitset<64> byte(*(a + (i*64)));
                myfile << byte;
            }
            myfile.close();
	    free(a);
        }
    }
}
