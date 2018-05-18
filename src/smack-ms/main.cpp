
/*
 *  smack-ms - split mapping check "Multisplice Edition"
 *
 *  Created by David Brawand on 04.05.10.
 *  Copyright 2010 UNIL. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <set>
#include <map>
#include <string>
#include <boost/progress.hpp> // BOOST PROGRESS

#include "smack-ms.hpp"

//namespace
using namespace std;

/* MAIN */
int main(int argc, char* argv[]) {
	/* USAGE */
	if (argc != 11) {
		//                                      1      2         3        4       5       6          7       8       9              10
		std::cerr << "Usage: " << argv[0] << " pairs genomemap splicemap cigars introns readlength minDist maxDist mapoutput_full mapoutput_spliced" << std::endl;
		exit(1);
	}
	std::cerr << std::endl << "  _______  __   __  _______  _______  ___   _  __  ";
	std::cerr << std::endl << " |       ||  |_|  ||   _   ||       ||   | | ||  | ";
	std::cerr << std::endl << " |  _____||       ||  |_|  ||       ||   |_| ||  |   Split";
	std::cerr << std::endl << " | |_____ |       ||       ||       ||      _||  |   Map";
	std::cerr << std::endl << " |_____  ||       ||       ||      _||     |_ |__|   Check";
	std::cerr << std::endl << "  _____| || ||_|| ||   _   ||     |_ |    _  | __    MS"; 
	std::cerr << std::endl << " |_______||_|   |_||__| |__||_______||___| |_||__| " << std::endl << std::endl;
	
	/* DECLARATIONS */
	char      line[200];
	int       scan[7];
	
	// the infolibs
	regiomap regions;
	intmap pairs;
	intmap mc;
	time_t startTime, endTime;

	int lico = 0;
	/* DATA READ STARTS HERE */
	int readl   = atoi(argv[6]);
	int minDist = atoi(argv[7]);
	int maxDist = atoi(argv[8]);
	
	  /**************/
	 /* READ MATES */
	/**************/
	
	// READ MATES
	std::ifstream pairFile(argv[1]);
	if (pairFile.fail()) {
		std::cerr << "Error: could not read from pairFile " << argv[1] << std::endl;
		exit(1);
	}
	lico = 0;
	std::cerr << "Reading matePair info from file " << argv[1] << std::endl;
	while (!pairFile.eof() && !pairFile.fail()) {
		pairFile.getline(line, 200);
		if (pairFile.eof()) break;
		sscanf(line, "%d%d", &scan[0], &scan[1]);
		if (++lico % 10000 == 0) std::cerr << "\r" << lico;
		pairs.insert(pair<int,int>(scan[0], scan[1]));
		pairs.insert(pair<int,int>(scan[1], scan[0]));
	}
	pairFile.close();
	std::cerr << "\r" << "pairFile: " << pairs.size() / 2 << " pairs (in " << lico << " lines)" << std::endl;	

	
	  /*****************/
	 /* READ MAPPINGS */
	/*****************/
	
	// READ FULL MAPPINGS
	std::ifstream mappingFile(argv[2]);
	lico = 0;
	std::cerr << "Reading genomeMappings from file " << argv[2] << std::endl;
	while (!mappingFile.eof() && !mappingFile.fail()) {
		mappingFile.getline(line, 200);
		if (mappingFile.eof()) break;
		sscanf(line, "%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4]); // readid, strand, region, start, mismatches
		if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		if (regions.find(scan[2]) == regions.end()) regions.insert(pair<int,Region*>(scan[2], new Region(scan[2],minDist,maxDist))); // create region
		if (pairs.find(scan[0]) == pairs.end()) pairs.insert(pair<int,int>(scan[0],0));
		regions.find(scan[2])->second->addMapping(scan, readl, pairs.find(scan[0])->second);
	}
	std::cerr << "\r" << "mappingFile: " << lico << " lines" << std::endl;
	mappingFile.close();
	
	
	// READ SPLICE INDEX
	splicemap spm;
	std::ifstream spliceFile(argv[4]);
	lico = 0;
	int dummy;
	std::cerr << "Reading Splice Index " << argv[4] << std::endl;
	while (!spliceFile.eof() && !spliceFile.fail()) {
		spliceFile.getline(line, 200);
		if (spliceFile.eof()) break;
		sscanf(line, "%d%d%d%d%d%d%d", &scan[0], &dummy, &scan[1], &scan[2], &scan[3], &scan[4], &scan[5]); //splice, centersplice, seq_region, strand, start, end, rightside intron (0 if end)
		if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		if (spm.find(scan[0]) == spm.end()) {
			spm[scan[0]] = new Splice(scan[0],scan[1],scan[2]); // splice,seqRegion,strand
		}
		spm.find(scan[0])->second->addSlice(scan[3], scan[4], scan[5]);
	}
	std::cerr << "\r" << "spliceFile: " << lico << " lines" << std::endl;
	spliceFile.close();
	
	
	// READ SPLICE MAPPINGS
	std::ifstream splicemappingFile(argv[3]);
	lico = 0;
	std::cerr << "Reading spliceMappings from file " << argv[3] << std::endl;
	while (!splicemappingFile.eof() && !splicemappingFile.fail()) {
		splicemappingFile.getline(line, 200);
		if (splicemappingFile.eof()) break;
		sscanf(line, "%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4]); // readid, strand, region, start, mismatches
		if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		int sreg = spm.find(scan[2])->second->region();
		if (regions.find(sreg) == regions.end()) regions.insert(pair<int,Region*>(sreg, new Region(sreg,minDist,maxDist))); // create region
		if (pairs.find(scan[0]) == pairs.end()) pairs.insert(pair<int,int>(scan[0],0));
		regions.find(sreg)->second->addMapping(scan, readl, pairs.find(scan[0])->second, spm.find(scan[2])->second);
	}
	std::cerr << "\r" << "splicemappingFile: " << lico << " lines" << std::endl;
	splicemappingFile.close();
	
	
	

	/***************/
	/* unify reads */
	/***************/
	std::cerr << std::endl << "Unifying read mappings..." << std::endl;
	time(&startTime);
	boost::progress_display show_progress( regions.size() );
	int map_c = 0, kill_c = 0;
	for (regiomap::const_iterator rg = regions.begin(); rg != regions.end(); ++rg) {
		++show_progress;
		for(readmap::const_iterator rd = rg->second->rmap.begin(); rd != rg->second->rmap.end(); ++rd) {
			map_c  += rd->second->countMappings();
			kill_c += rd->second->unify(); // will only keep best matching
			mc[rd->second->rid()] += rd->second->countMappings(); // count read mappings AFTER unify
		}
		rg->second->indexMappings();
	}
	
	int acceptedMappings = map_c - kill_c;
	time(&endTime);
	double finaltime = difftime (endTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed  (" << kill_c << " mappings out of " << map_c << " were rejected, " << acceptedMappings << " were kept)" << std::endl;
	
	// READ INTRONS
	std::ifstream intronFile(argv[5]);
	if (intronFile.fail()) {
		std::cerr << "Error: could not read from intronFile " << argv[5] << std::endl;
		exit(1);
	}
	lico = 0;
	std::cerr << "Reading introns from file " << argv[5] << std::endl;
	while (!intronFile.eof() && !intronFile.fail()) {
		intronFile.getline(line, 200);
		if (intronFile.eof()) break;
		sscanf(line, "%d%d%d", &scan[0], &scan[1], &scan[2]);
		if (regions.find(scan[0]) != regions.end()) {
			regions.find(scan[0])->second->addIntron(scan[1], scan[2]); // add mapping check mapping compatibility (overwrites mapping)
			if (++lico % 10000 == 0) std::cerr << "\r" << lico;
		}
	}
	intronFile.close();
	std::cerr << "\r" << "intronFile: " << regions.size() << " regions with " << lico << " introns" << std::endl;	
	
	
	// INIT STATS
	int statsize = abs(maxDist-minDist)+1;
//	cerr << "1 " << statsize << endl;
	int * sumdist = new int[statsize];
//	cerr << "2 " << endl;
	for (int i=0; i < statsize; i++) sumdist[i] = 0;
//	cerr << "3 " << endl;
	float totalMap = 0, totalAcc = 0, totalSingles = 0;
//	cerr << "4 " << endl;
	

	std::ofstream fullOut(argv[9]);
//	cerr << "5 " << endl;
	std::ofstream spliceOut(argv[10]);
	
	
	std::cerr << std::endl << std::endl << "CHR\tMAPS\tACC\tINSERT\tSINGLES" << std::endl;
	for (regiomap::const_iterator rg = regions.begin(); rg != regions.end(); ++rg) {
		// RESOLVE PAIRS
//		cerr << "A" << endl;
		rg->second->markInRangeMappings(readl,&mc); // max and min are stored in region

		// OUTPUT accepted mappings
//		cerr << "B" << endl;
		for (imap::const_iterator it = rg->second->mmm.begin(); it != rg->second->mmm.end(); ++it) {
			if (it->second->isAccepted()) {
				if (it->second->spliced()) spliceOut << it->second->getString() << endl;
 				else                       fullOut   << it->second->getString() << endl;
			}
		}
		// STATS
//		cerr << "C" << endl;
		std::cout << rg->second->rid() << "\t" << rg->second->totalMappings() << "\t" << rg->second->acceptedMappings() << "\t" << minDist + rg->second->getDistanceMode(statsize) << "\t" << rg->second->singles() << std::endl;
//		cerr << "D" << endl;
		intvec rdist = rg->second->getDistanceDistributon();
//		cerr << "E" << endl;
		for (int i=0; i < statsize; i++) sumdist[i] += rdist[i];
//		cerr << "F" << endl;
		totalAcc     += rg->second->acceptedMappings();
//		cerr << "G" << endl;
		totalMap     += rg->second->totalMappings();
//		cerr << "H" << endl;
		totalSingles += rg->second->singles();
//		cerr << "I" << endl;
	}
 
	// print global distribution (normalized by max mode)
	double distMean = 0, distSD = 0, distCount = 0;
	int maxd = 0, maxm = 0;
	for (int i=0; i <= maxDist - minDist; i++) {
		//cerr << ">>" << i << "=" << sumdist[i] << endl;
		distMean  += (i+minDist)*sumdist[i];
		distCount += sumdist[i];
		if (sumdist[i] > maxd) {
			maxd = sumdist[i];
			maxm = i;
		}
	}
	std::cout << setprecision(1);
	std::cout << "TOTAL\t";
	std::cout << totalMap/1000 << "k\t";
	std::cout << totalAcc/1000 << "k\t"; 
	std::cout << minDist+maxm << "\t";
	std::cout << totalSingles << std::endl;
	
	// calc STDEV
	distMean /= distCount;
	for (int i=0; i <= maxDist - minDist; i++) {
		distSD += sumdist[i] * (i+minDist - distMean) * (i+minDist - distMean);
	}	
	distSD /= distCount;
	double stdev = sqrt(distSD);
	
	std::cout << setprecision(2);
	std::cout << "Mate pair distance (relative to end): " << distMean << "±" << stdev << endl; 
	
	//if (4*stdev+distMean > maxDist)                         cerr << endl << "WARNING: Consider increasing insert size threshold to at least " << int(4*stdev+distMean) << "!" << endl << endl;
	//if (distMean-4*stdev < minDist && distMean-4*stdev > 0) cerr << endl << "WARNING: Consider decreasing insert size threshold to at least " << int(distMean-4*stdev) << "!" << endl << endl;
	
	float stars = 0;
	std::cout << std::endl << "Insert size distribution" << std::endl; 
	for (int i=0; i <= maxDist - minDist; i++) {
		 std::cout << i+minDist << "\t";
		 stars = (float)50*sumdist[i]/maxd;
		 for (int j=0; j < (int)stars; j++) std::cout << "*";
		 std::cout << std::endl;
	}
	
	
	exit(0);	// SUCCESS
} // end of main

