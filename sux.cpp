// g++ -g -std=c++17 -Wall -Wextra -O0 -march=native -I ./ -I ~/include sux.cpp -fsanitize=address -fsanitize=undefined
#include <sux/bits/Rank9Sel.hpp>


using namespace sux::bits;
using namespace std;

int main() {

    size_t size = 1024; // bits
	uint64_t bitvect[size / 64 + 1];
    memset(bitvect, 0, (size / 64 + 1) * 8);
    cout << sizeof(bitvect) << endl;

	//for (size_t i = 1; i < size; i++) bitvect[i / 64] |= UINT64_C(1) << i % 64;
    (*bitvect) += 32;
    (*bitvect) += 4;
    // v is a bit vector represented by an array of uint64_t and n is the number of bits represented therein
    // v = bitvect; n = size
	Rank9Sel newName(bitvect, size);

    cout << "Sizeof bitector: " << sizeof(bitvect) << " Sizeof rank: " << sizeof(newName) << endl;
    cout << "Rank(6): " << newName.rank(6) << endl;
    cout << "Rank(4): " << newName.rank(4) << endl;
    cout << "Rank(3): " << newName.rank(3) << endl;
    cout << "Rank(2): " << newName.rank(2) << endl;
    // Sleect is 0 based
    cout << "Select(1): " << newName.select(0) << endl;
    cout << "Select(2): " << newName.select(1) << endl;

}