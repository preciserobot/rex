#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <set>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>

using namespace std;

// swamples given number of ids from a redundant list

class RandomNumber {
public:
    RandomNumber() {
        srand(time(0));
    }
    
    int operator() ( int n) {
        return (int)((double)n * rand()/(RAND_MAX+1.0));
    }
};


/* MAIN */
int main(int argc, char **argv) {
	/* USAGE */
	if (argc != 3) {
		std::cerr << "Usage: " << argv[0] << " idlist number" << std::endl;
		exit(1);
	}
	std::cerr << std::endl << std::endl << "    SWAMPLE 0.1 " << std::endl << std::endl << "Loading data into memory..." << std::endl;
		
	/* DECLARATIONS */
	typedef set<int> intset;
	intset IDs;
	typedef vector<int> intvec;
	intvec shuffleIDs;
	char line[200];
	int  scan[2];
	int  lico = 0;
	int  idSize = 0;
	
	std::ifstream idFile(argv[1]);
	// READ EXON FILE
	if (idFile.fail()) {
		std::cerr << "Error: could not read from ID file " << argv[1] << std::endl;
		exit(1);
	}  
	std::cerr << "Reading IDs from file " << argv[1] << std::endl;
	while (!idFile.eof()) {
		idFile.getline(line, 200);
		if (idFile.eof()) break;
		if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		sscanf(line, "%d", &scan[0]); // ID
		IDs.insert(scan[0]);
		if (!idFile.fail()) ++idSize; // count sucessfully read line
	}
	idFile.close();
	// stats
	std::cerr << "\r" << "idFile: " << idSize << " IDs read" << std::endl;
	
	// make a shuffeled vector
	for (intset::const_iterator itd = IDs.begin(); itd != IDs.end(); ++itd) shuffleIDs.push_back(*itd);
	// random generator function
	//ptrdiff_t myrandom (ptrdiff_t i) { return rand() % i; }
	//ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;
	//std::random_shuffle( shuffleIDs.begin(), shuffleIDs.end(), p_myrandom ); // shuffle elements of v
	
	RandomNumber rnd;
	std::random_shuffle( shuffleIDs.begin(), shuffleIDs.end(), rnd ); // shuffle elements of v

	// printing the results
	int sw = atoi(argv[2]);
	if (sw > shuffleIDs.size()) {
		std::cerr << "WARNING: sample size exceeds mapped read number => writing a unique list" << std::endl;
		sw = shuffleIDs.size();
	}
	for (int i = 0; i < sw; ++i) {
		std::cout << shuffleIDs[i] << std::endl;
	}
	return EXIT_SUCCESS; 
} // end of main

