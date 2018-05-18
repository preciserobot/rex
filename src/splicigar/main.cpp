
/*
 *  splicigar - returns cigar lines of possible splice juntion neighboring exons using a given half length
 *
 *  Created by David Brawand on 12.02.10.
 *  Added exon stringency filter on 14.04.10.
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
#include "splicigar.hpp"

//namespace
using namespace std;

/* MAIN */
int main(int argc, char* argv[]) {
	/* USAGE */
	if (argc != 6) {
		std::cerr << "Usage: " << argv[0] << " junctions (exonfile/-) overhang readlength minimum_intron_size --> writes to COUT " << std::endl;
		exit(1);
	}
	std::cerr << std::endl << "                           ";
	std::cerr << std::endl << "              spliCigar    ";
	std::cerr << std::endl << "                 sSs       ";
	std::cerr << std::endl << "   .-------------.S        ";
	std::cerr << std::endl << "  '._|0.6|_______.\"       "; 
	std::cerr << std::endl << "                           ";
	std::cerr << std::endl << "                           " << std::endl << std::endl;
	
	/* DECLARATIONS */
	char      line[200];
	int       scan[5];
	
	// the infolibs
	regiomap regions;
	intmap pairs;
	
	int lico = 0;
	/* DATA READ STARTS HERE */
	std::ifstream intronFile(argv[1]);
	std::ifstream exonFile(argv[2]);
	
	int overh = atoi(argv[3]);
	int readl = atoi(argv[4]);
	int mintron = atoi(argv[5]);
	if (overh >= readl) {
		cerr << "ERROR: OVERHANG greater or equal to READLENGTH" << endl;
		abort();
	}
	int distance = readl - overh;
	int idoffset = 0;
	int skippedIntrons = 0, acceptedIntrons = 0;
	// READ INTRONS
	if (intronFile.fail()) {
		std::cerr << "Error: could not read from intronFile " << argv[1] << std::endl;
		exit(1);
	}
	lico = 0;
	std::cerr << "Reading introns from file " << argv[1] << std::endl;
	while (!intronFile.eof() && !intronFile.fail()) {
		intronFile.getline(line, 200);
		if (intronFile.eof()) break;
		sscanf(line, "%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4]); // id, region, strand, start, end
		if (regions.find(scan[1]) == regions.end()) {
			regions.insert(pair<int,Region*>(scan[1], new Region(scan[1], distance))); // create region
		}
		if (scan[4] - scan[3] + 1 >= mintron) { // check size
			regions.find(scan[1])->second->addIntron(scan); // add mapping check mapping compatibility (overwrites mapping)
			acceptedIntrons++;
		}
		else {
			skippedIntrons++;
		}
		std::cerr << "\r" << ++lico;
	}
	intronFile.close();
	std::cerr << "\r" << "intronFile: " << regions.size() << " regions with " << acceptedIntrons << " introns (" << skippedIntrons << " did not pass length filter of " << mintron << " bases)" << std::endl;	

	// stringency?
	bool stringency;
	if (exonFile.fail()) {
		stringency = false;
	}
	else {
		lico = 0;
		std::cerr << "Reading exons from file " << argv[2] << std::endl;
		while (!exonFile.eof() && !exonFile.fail()) {
			exonFile.getline(line, 200);
			if (exonFile.eof()) break;
			sscanf(line, "%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3]);
			if (regions.find(scan[0]) == regions.end()) {
				regions.insert(pair<int,Region*>(scan[0], new Region(scan[0], distance))); // create region
			}
			regions.find(scan[0])->second->addExon(scan); // add mapping check mapping compatibility (overwrites mapping)
			std::cerr << "\r" << ++lico;
		}
		exonFile.close();
		std::cerr << "\r" << "exonFile: " << lico << " exons" << std::endl;
		stringency = true;
	}

	int maxi = 0;
	int cigs = idoffset;
	// run filter
	for(regiomap::const_iterator it = regions.begin(); it != regions.end(); ++it){
		cigs = it->second->makeCigars(cigs, stringency, readl);
		std::cerr << "Region " << it->first << " -> " << cigs - idoffset << std::endl;
		for (intmap::const_iterator aa = it->second->stats.begin(); aa != it->second->stats.end(); ++aa) {
			maxi = (aa->first > maxi) ? aa->first : maxi;
		}
	}
	//HEADER
	std::cerr << "SPLICECLASS\t"; for (unsigned int s = 1; s <= maxi; s++) std::cerr << s << "\t"; std::cerr << std::endl;
	//REGIONSTATS
	for(regiomap::const_iterator it = regions.begin(); it != regions.end(); ++it){
		// stats
		if (stringency) std::cerr << "R" << it->first << " STRINGENT\t";
		else std::cerr << "R" << it->first << "\t\t";
		for (unsigned int s = 1; s <= maxi; s++) {
			if (it->second->stats.find(s) != it->second->stats.end()) std::cerr << it->second->stats.find(s)->second << "\t";
			else std::cerr << "0\t";
		}
		std::cerr << std::endl;
	}
	
	// WRITE OUTPUT
	std::cerr << "Writing Output...";
	for(regiomap::const_iterator it = regions.begin(); it != regions.end(); ++it){
		for(cigarmultimap::const_iterator rg = it->second->activecigars->begin(); rg != it->second->activecigars->end(); ++rg){
			int centralintron = rg->second->centerIntron(distance);
			std::cout << rg->second->cigarid() << "\t" << centralintron << "\t" << it->first << "\t" << rg->second->cigarstrand() << "\t" << rg->second->starts() << "\t";
			for(intronmap::const_iterator sl = rg->second->slices.begin(); sl != rg->second->slices.end(); ++sl){
				std::cout << sl->second->from() - 1 << "\t" << sl->second->intronid() << std::endl << rg->second->cigarid() << "\t" << centralintron << "\t" << it->first << "\t" << rg->second->cigarstrand() << "\t" << sl->second->to() + 1 << "\t";
			}
			std::cout << rg->second->ends() << "\t0" << std::endl;
		}
	}
	std::cerr << "DONE" << std::endl;

	
	
	exit(0);	// SUCCESS
} // end of main
