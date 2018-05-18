#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <map>
#include <set>
#include "ames-ms.hpp"

using namespace std;

// Ambiguity estimator script
// Maps mockreads to exons
// estimates uniqueness for each exon (unambiguous mapping positions)
// return uniqueness score (exon uniscore)

#define MAXEXONS 2000000 //the maximal size of the exon catalog
#define MAXMAPPINGS 300000000 //maximum number of exons in one transcript
#define MAXREGIONS 30000 //maximum number of region references (chromosomes / contigs)

/* MAIN */
int main(int argc, char **argv) {
	/* USAGE */
	if (argc != 7) {
		std::cerr << "Usage: " << argv[0] << " exons fullmap splicemap splicecatalog readlength amesout" << std::endl;
		exit(1);
	}
	std::cerr << std::endl << "    *              _____          " << std::endl;
	std::cerr <<              "   / \\   _ __ ___ | ____|___     " << std::endl;
	std::cerr <<              "  / _ \\ | \'_ ` _ \\|  _| / __|  Ambiguity" << std::endl;
	std::cerr <<              " / ___ \\| | | | | | |___\\__ \\  Estimator" << std::endl;
	std::cerr <<              "/_/   \\_\\_| |_| |_|_____|___/  MultiSplice" << std::endl << std::endl;
	
	std::cerr << "Loading data into memory..." << std::endl;
	
	
	/* DECLARATIONS */
	const int exonLen   = MAXEXONS;
	const int readLen   = MAXMAPPINGS;
	const int regionLen = MAXREGIONS;
	char      line[200];
	
	Exon   **exons   = new Exon *[exonLen];
	Read   **reads   = new Read *[readLen];
	Region **regions = new Region *[regionLen];
	
	int  exonSize   = 0;
	int	readSize   = 0;
	int	regionSize = 0;
	
	int  exonMax    = 0; // the highest exon id (for iterators)
	
	int readl       = atoi(argv[5]);

	int scan[8];
	int lico = 0;
	int readCount = 0; // counts reads
	
	//start timer
	time_t startTime, currentTime;

	std::ifstream exonFile(argv[1]);
	// READ EXON FILE
	if (exonFile.fail()) {
		std::cerr << "Error: could not read from exon file " << argv[1] << std::endl;
		exit(1);
	}  
	std::cerr << "Reading exons from file " << argv[1] << std::endl;
	while (!exonFile.eof()) {
		exonFile.getline(line, 200);
		if (exonFile.eof()) break;
		if (exonSize < exonLen) {
			if (++lico % 1000 == 0) std::cerr << "\r" << lico;
			sscanf(line, "%d%d%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4], &scan[5], &scan[6]); // gene, exon, region, start, end, strand, core
			exons[scan[1]] = new Exon(scan);
			
			
			if (scan[1] > exonMax) exonMax = scan[1]; // update exonMax
			// ADD Exon to region 
			bool added = false;
			for(int r=0; r < regionSize; r++) {
				if (regions[r]->region() == exons[scan[1]]->region()) {
					regions[r]->addExonToGene(exons[scan[1]]);
					added = true;
				}
			}
			// if no assignment create a region first
			if (added == false) {
				regions[regionSize] = new Region(exons[scan[1]]->region());
				regions[regionSize]->addExonToGene(exons[scan[1]]);
				++regionSize;
				if (!(regionSize < regionLen)) {
					std::cerr << "Out of bounds (regionLen)!" << std::endl;
				}
			}

		
			//if (scan[3] == 132633561) {
			//	std::cerr << std::endl;
			//	std::cerr << line << std::endl;
			//	std::cerr << scan[0] << "\t" << scan[1] << "\t" << scan[2] << "\t" << scan[3] << "\t" << scan[4] << "\t" << scan[5] << "\t" << scan[6] << std::endl;
			//	exons[scan[1]]->show();
			//}
			
		
		} else {
			std::cerr << "Out of bounds (exonLen)!" << std::endl;
		}
		if (!exonFile.fail()) exonSize++; // count sucessfully read line -> limits array size
	}
	exonFile.close();
	// stats
	std::cerr << "\r" << "exonFile: " << exonSize << " records from " << regionSize << " contiguous regions" << std::endl;
	for (int r = 0; r < regionSize; r++) regions[r]->orderGenesByStart(); //order and index genes (NECESSARY FOR DICOSEARCH! in projection method)
	exonFile.close();
	
	

	
	// READ SPLICE INDEX
	splicemap spm;
	std::ifstream spliceFile(argv[4]);
	lico = 0;
	int dummy; // the centersplice is ignored
	std::cerr << "Reading Splice Index " << argv[4] << std::endl;
	while (!spliceFile.eof() && !spliceFile.fail()) {
		spliceFile.getline(line, 200);
		sscanf(line, "%d%d%d%d%d%d%d", &scan[0], &dummy, &scan[1], &scan[2], &scan[3], &scan[4], &scan[5]); //splice, centersplice, seq_region, strand, start, end, rightside intron (0 if end)
		if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		if (spm.find(scan[0]) == spm.end()) {
			spm[scan[0]] = new Splice(scan[0],scan[1],scan[2]); // splice,seqRegion,strand
		}
		spm.find(scan[0])->second->addSlice(scan[3], scan[4]);
	}
	std::cerr << "\r" << "spliceFile: " << lico << " lines" << std::endl;
	spliceFile.close();
	


	
	// SET EXON SPLICE
	std::cerr << "Setting exon anchors..." << std::endl;
	lico = 0;
	int * intron = new int[3];
	for (splicemap::const_iterator sp = spm.begin(); sp != spm.end(); ++sp) {
		if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		if (sp->second->centralintron(intron)) {
			for (int r = 0; r < regionSize; r++) {
				if (intron[0] == regions[r]->region()) {
					regions[r]->introns.insert(pair<int,int>(intron[1],intron[2]));
				}
			}
		}
		else { cerr << "HUH? Couldn't find central intron -> Check Splice library consistency"; abort(); }
	}
	std::cerr << "\r " << lico << " splicigars analyzed" << std::endl;
	int lonelySplice = 0;
	lico = 0;
	for (int r = 0; r < regionSize; r++) {
		for(intmap::const_iterator in = regions[r]->introns.begin(); in != regions[r]->introns.end(); ++in) {
			if (regions[r]->spliceProjection(in->first, in->second)) {
				//cerr << endl << "SUCESS           : projected splice         " << intron[0] << ":" << intron[1] << "-" << intron[2] << " to exon" << endl;
				cerr << "\r" << ++lico << " ";
			}
			else {
				++lonelySplice;
				//cout << regions[r]->region() << "\t" << in->first << "\t" << in->second << endl;
				//cout << "grep " << in->first << " alljunc | grep " << in->second << " >> lonelyjunc" << endl;
				//abort();
			}
		}
	}
	std::cerr << "\r" << "SpliceJunctions: " << lico << " (" << lonelySplice << " junctions were lonely)" << std::endl;
	
	
	
	// READ FULL MAPPINGS
	std::ifstream mappingFile(argv[2]);
	lico = 0;
	std::cerr << "Reading genomeMappings from file " << argv[2] << std::endl;
	while (!mappingFile.eof() && !mappingFile.fail()) {
		mappingFile.getline(line, 200);
		if (readSize < readLen) {
			sscanf(line, "%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3]); // readid, strand, region, start
			if (++lico % 1000 == 0) std::cerr << "\r" << lico;
			if (reads[scan[0]] == NULL) {
				reads[scan[0]] = new Read(scan[0]); // create read (allocate memory)
				if (readSize < scan[0]) readSize = scan[0]; 
				readCount++;
			}
			reads[scan[0]]->addFullMapping(scan[2], scan[1], scan[3], readl);
		} else {
			if (lico == 0) std::cerr << "Error: could not read from mappingFile " << argv[2] << std::endl;
			else std::cerr << "Out of bounds " << readSize << " (readLen)!" << std::endl;
			exit(1);
		}
	}
	std::cerr << "\r" << "mappingFile: " << lico << " lines" << std::endl;
	mappingFile.close();
	
	
	
	// READ SPLICE MAPPINGS
	std::ifstream splicemappingFile(argv[3]);
	lico = 0;
	std::cerr << "Reading spliceMappings from file " << argv[3] << std::endl;
	while (!splicemappingFile.eof() && !splicemappingFile.fail()) {
		splicemappingFile.getline(line, 200);
		if (readSize < readLen) {
			sscanf(line, "%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3]); // readid, strand, region, pos
			if (++lico % 1000 == 0) std::cerr << "\r" << lico;
			if (reads[scan[0]] == NULL) {
				reads[scan[0]] = new Read(scan[0]); // create read (allocate memory)
				if (readSize < scan[0]) readSize = scan[0]; 
				readCount++;
			}
			reads[scan[0]]->addSpliceMapping(scan[1], scan[3], readl, spm.find(scan[2])->second);
		} else {
			if (lico == 0) std::cerr << "Error: could not read from mappingFile " << argv[2] << std::endl;
			else std::cerr << "Out of bounds " << readSize << " (readLen)!" <<readLen << std::endl;
			exit(1);
		}
	}
	std::cerr << "\r" << "splicemappingFile: " << lico << " lines" << std::endl;
	splicemappingFile.close();
	
	
	// resolve all reads for ambiguity
	lico = 0;
	std::cerr << std::endl << "Determining uniqueness of reads and adding to exons..." << std::endl;
	time(&startTime);
	int mismap = 0, successMaps = 0, totalMaps = 0, totalReads = 0;
	for (int j = 0; j < readSize; j++) {
		if (reads[j] == NULL) continue;
		totalReads++;
		//if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		if (reads[j]->isUnique() || reads[j]->unify()) {
			Mapping * tMap = reads[j]->trustedMapping();
			if (j % 10000 == 0) { time(&currentTime); double estim = (difftime(currentTime,startTime) / j) * (readSize - j); std::cerr << "\r ~ " << (int)estim << " seconds remaining"; }
			bool unknownRegion = true;
			for (int r = 0; r < regionSize; r++) {
				if (tMap->regiomap(regions[r])) {
					unknownRegion = false;
					if (regions[r]->projection(tMap, readl)) successMaps++;
					totalMaps++;
				}
			}
			if (unknownRegion) {
				mismap++;
				std::cerr << "\rREGION " << tMap->region() << "(" << tMap->start() << "-" << tMap->end() << ") ?" << std::endl;
			}
		}
		//NO USE AS WILL FINISH QUICKLY ANYWAY// delete reads[j]; // IMPLEMENT EFFICIENT DESTRUCTOR!!!
	}
	time(&currentTime);
	double finaltime = difftime (currentTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl << std::endl;
	
	std::cerr << "                       Total Reads : " << totalReads << std::endl;
	std::cerr << "              with unique mappings : " << totalMaps << std::endl << std::endl;
	std::cerr << "       Successful mappings on exon : " << successMaps << std::endl;
	std::cerr << "         On annotationless Regions : " << mismap << std::endl;
	
	// printing the results
	time(&startTime);
	std::cerr << std::endl << "Analyzing exon uniqueness and writing " << exonSize << " Exons (" << argv[6] << ")..." << std::endl;
	std::ofstream outFile(argv[6]);
	for(int e = 0; e <= exonMax; e++) {
		if (e % 10000 == 0) {
			time(&currentTime);
			double dif = difftime (currentTime,startTime);
			double estim = (dif / e) * (exonSize - e);
			std::cerr << "\r ~ " << (int)estim << " seconds remaining";
		}
		if (exons[e] == NULL) continue;
		outFile << exons[e]->gid() << "\t";
		outFile << exons[e]->eid() << "\t";
		outFile << exons[e]->region() << "\t";
		outFile << exons[e]->start() << "\t";
		outFile << exons[e]->end() << "\t";
		outFile << exons[e]->str() << "\t";
		outFile << exons[e]->tags() << "\t";
		outFile << exons[e]->uniqEstimate(readl) << "\t" << exons[e]->uniqEstimateBrute(readl) << std::endl;
	}
	outFile.close();
	time(&currentTime);
	finaltime = difftime (currentTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl << std::endl;
	
	exit(0); // Success returns zero (should be verified by wrapping script)
} // end of main

