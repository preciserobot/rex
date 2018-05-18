
/*
 *  geco - the gene coverage evaluator
 *
 *  Created by David Brawand on 03.02.10.
 *  Copyright 2010 UNIL. All rights reserved.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <set>
#include "geco.hpp"
#include <boost/progress.hpp> // BOOST PROGRESS

// implement method to find smallest positive and negative slope
using namespace std;

/* MAIN */
int main(int argc, char* argv[]) {
	
	std::cerr << std::endl << "                          )/_       "    << "                                "; 
	std::cerr << std::endl << "                _.--..---\"-,--c_    "   << "  ______  ______  _____  _____  ";
	std::cerr << std::endl << "           \\L..'           ._O__)_  "   << " |  ____ |______ |      |     | ";
	std::cerr << std::endl << "   ,-.     _.+  _  \\..--( /         "   << " |_____| |______ |_____ |_____| "; 
    std::cerr << std::endl << "     `\\.-''__.-' \\ (     \\_         " << "                                ";
	std::cerr << std::endl << "       `'''       `\\__   /\\         "  << " GeneCoverage Estimator (QC)    ";
	std::cerr << std::endl << "                   ')              "     << "                                " << std::endl << std::endl;
	
	if (argc != 4) {
		std::cerr << "\tUsage: " << argv[0] << " exons coverage output" << std::endl << std::endl;
		exit(1);
	}
	
	// START GLOBAL TIMER
	//	time_t beginTime, endTime, startTime, currentTime;
	//	time(&beginTime);
	//	double finaltime;
	
	// global variables
	typedef map<int, Region*> regmap;
	regmap regions;
	int    scan[8];
	int    dummy;
	char   line[200];
	int lico = 0, success = 0, failed = 0, totalGenes = 0;
	// READ EXON FILE
	std::ifstream exonFile(argv[1]);
	if (exonFile.fail()) {
		std::cerr << "Error: could not read from exon file " << argv[1] << std::endl;
		exit(1);
	}
	std::cerr << "Reading annotation from file (" << argv[1] << ")" << std::endl;
	while (!exonFile.eof()) {
		exonFile.getline(line, 200);
		if (exonFile.eof()) break;
		if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		sscanf(line, "%d%d%d%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4], &scan[5], &scan[6], &scan[7]); // gene, exon, region, start, end, strand, tag/transcript, score
		if (regions.find(scan[2]) == regions.end()) regions.insert(pair<int,Region*>(scan[2], new Region(scan[2])));
		if (regions.find(scan[2])->second->addExon(scan)) success++;
		else failed++;
	}
	exonFile.close();
	if (failed) { std::cerr << endl << "ERROR: Not all annotation lines could be read" << endl; exit(1); }
	for (regmap::const_iterator it = regions.begin(); it != regions.end(); ++it) totalGenes += it->second->geneCount();
	std::cerr << "\r" << regions.size() << " regions with " << totalGenes << " genes (" << success << " annotation lines)" << std::endl << std::endl;
	
	// count total gene number
	for (regmap::const_iterator it = regions.begin(); it != regions.end(); ++it) {
		it->second->setGeneCount();
		it->second->bucketizeGenes();
	}
	
	/*********************
		
		// getBuckets(start, end) // returns a map to gene pointers // several passes should not affect matching count or pocessed genes (unique set)
		// search only matching genes in buckets
		// create better destuctors, garabage collection
		// calc RPKM -> find best transcript and rerun?
		// add auto bucketsize finder (longest gene length)
		// check if in raparm distribution is done for all exons of a gene or non-tagged exons are ignored? is is mathematically optimal?
	 
	// transcript check
	for(regmap::const_iterator rg = regions.begin(); rg != regions.end(); ++rg) {
		for(genemap::const_iterator ge = rg->second->genes.begin(); ge != rg->second->genes.end(); ++ge) {
			for(transcriptmap::const_iterator tr = ge->second->transcripts.begin(); tr != ge->second->transcripts.end(); ++tr) {
				cerr << rg->second->rid() << "\t";
				cerr << ge->second->gid() << "\t";
				cerr << tr->second->tid() << "\t";
				cerr << tr->second->rankedexons.size() << endl;
			}
		}
	}
	exit(1);
	 **********************/ 
		
	// READ COVERAGE FILE
	lico = 0;
	Region * previousRegion = NULL;
	Region * currentRegion;
	long int geneCoverage = 0, fullCoverage = 0, finishedGeneCount = 0, allGeneCount = 0;
	intset currentGenes, previousGenes, coveredRegions, noregion;
	// open file and read
	std::ifstream coverageFile(argv[2]);
	if (coverageFile.fail()) {
		std::cerr << "Error: could not read from coverage file " << argv[2] << std::endl;
		exit(1);
	}
	std::cerr << "Reading coverage from file (" << argv[2] << ")" << std::endl;
	while (!coverageFile.eof()) {
		coverageFile.getline(line, 200);
		if (coverageFile.eof()) break;
		sscanf(line, "%d%d%d%d%d%d", &scan[0], &dummy, &scan[1], &scan[2], &scan[3], &scan[4]); // region, (strand,) start, end, length, coverage
		coveredRegions.insert(scan[0]);
		// fetch region
		if (regions.find(scan[0]) != regions.end()){
			currentRegion = regions.find(scan[0])->second;
			// if new region set remaining genes as finished (and invoke covanalyze?)
			if ((previousRegion != NULL) && (currentRegion != previousRegion)) {
				for (intset::const_iterator it = previousGenes.begin(); it != previousGenes.end(); ++it) {
					previousRegion->addFinishedGene(*it);
				}
				previousGenes.clear();
				finishedGeneCount += previousRegion->finishedGeneCount();
				allGeneCount      += previousRegion->geneCount();
				std::cerr << "\r                                           \r" << lico << "\r\t\t" << scan[0] << "\t" << scan[2];
			}
			previousRegion = currentRegion;

			// STDERR
			if (++lico % 10000 == 0) std::cerr << "\r" << lico << "\r\t\t" << scan[0] << "\t" << scan[2];
			
			// check if line is new
			if (!(currentRegion->isFollowingSlice(scan[1], scan[2]))) continue;
			
			// addCoverageSlice returns affected genes' names
			currentGenes = currentRegion->addCoverageSliceUsingBuckets(scan);
			// sum up coverage
			if (currentGenes.size() > 0) geneCoverage += scan[3] * scan[4];
			fullCoverage += scan[3] * scan[4];
			// finished genes (push to finished and remove from set)
			for (intset::const_iterator it = previousGenes.begin(); it != previousGenes.end(); ++it) {
				if (currentGenes.count(*it) == 0) { // not in current set
					currentRegion->addFinishedGene(*it);
					previousGenes.erase(it); // remove from previous set
				}
			}
			// add new genes to previous
			for (intset::const_iterator it = currentGenes.begin(); it != currentGenes.end(); ++it) {
				previousGenes.insert(*it); // will be ignored if exists as set has unique keys
			}
			currentGenes.clear();
		}
		else { noregion.insert(scan[0]); continue; } // unknown region
	}
	std::cerr << "\r" << lico << " lines with " << fullCoverage << " bases covered (" << geneCoverage << " in genes)" << std::endl;
	if (noregion.size() > 0) {
		std::cerr << "\rNOTICE: Regions ";
		for (intset::const_iterator rr = noregion.begin(); rr != noregion.end(); ++rr) std::cerr << "(" << *rr << ") ";
		std::cerr << "have no annotation" << std::endl << std::endl;
	}
	std::cerr << std::endl;
	coverageFile.close();

	// after all read or region by region..?
	std::ofstream outFile(argv[3]);
	outFile << setprecision(5);
	boost::progress_display show_progress( regions.size() );
	for(regmap::const_iterator it = regions.begin(); it != regions.end(); ++it){
		++show_progress;
		if (coveredRegions.find(it->second->rid()) != coveredRegions.end()) {
			it->second->evaluateAndPrint(fullCoverage, &outFile);
		}
	}
	cerr << "\r                " << endl;
	outFile.close();
	// calc gene efficiency (sequenced bases on genes sum(Region::geneCoverage) vs fullCoverage
/*	double gecoscore;
	double rgeco[2];
	for(regmap::const_iterator it = regions.begin(); it != regions.end(); ++it){
		if (coveredRegions.find(it->second->rid()) != coveredRegions.end()) {
			cerr << "\r(R" << it->second->rid() << ")";
			it->second->getGecoScore();
		}
	}
*/	
	
	exit(0);	// SUCCESS
} // end of main

