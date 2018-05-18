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
/* #include <boost/regex.hpp> // BOOST REGEX */ 
#include <boost/progress.hpp> // BOOST PROGRESS

using namespace std;


/* DEFINITIONS (TRIAL AND ERROR) */
#define MAXEXONS 2000000    // maximal number of exons
#define MAXREADS 100000000  // maximal number of reads
#define MAXREGIONS 30000    // maximal number of region references (chromosomes / contigs)
#define MAXREADMAPPINGS 200 // maximal number of mappings per read

/* round double to int */
inline int rounddouble(double x) { return static_cast<int>(x + (x > 0.0 ? +0.5 : -0.5)); }

/* strand compare */
inline bool strandcmp(int a, int b) { return (a != b && a + b == 0) ? false : true; }

#include "raparm-ms.hpp"

/* MAIN */
int main(int argc, char* argv[]) {
	/* USAGE */
	if (argc != 16) {
		std::cerr << "Usage: " << argv[0] << " exons genomemap splicemap splicigars readlength tag strands extend bestMatch (maps rpkms islands coverage unique splicefreq)" << std::endl;
		exit(1);
	}
	std::cerr << std::endl << "    _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_| Retro";
	std::cerr << std::endl << "    _|                                                                                          _| Aware,"; 
	std::cerr << std::endl << "    _|    _|_|      _|_|      _|_|      _|_|      _|_|    _|      _|                            _| Probabilistic,";
	std::cerr << std::endl << "    _|  _|    _|  _|    _|  _|    _|  _|    _|  _|    _|  _|_|  _|_|                            _| Ambiguous";
	std::cerr << std::endl << "    _|  _|_|_|    _|_|_|_|  _|_|_|    _|_|_|_|  _|_|_|    _|  _|  _|        _|_| _|_|     _|_|  _| Read";
	std::cerr << std::endl << "    _|  _|   _|   _|    _|  _|        _|    _|  _|   _|   _|      _|  _|_|  _|  _|  _|    _|    _| Mapping over";
	std::cerr << std::endl << "    _|  _|    _|  _|    _|  _|        _|    _|  _|    _|  _|      _|        _|      _|  _|_| 64 _| Multiple";
	std::cerr << std::endl << "    _|                                                                                          _| Splicesites"; 
	std::cerr << std::endl << "    _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_| 64-bit" << std::endl << std::endl;
	
	// START GLOBAL TIMER
	time_t beginTime, endTime, startTime, currentTime;
	time(&beginTime);
	double finaltime;
	
	
	/* DECLARATIONS */
	const int exonLen   = MAXEXONS;
	const int readLen   = MAXREADS;
	const int regionLen = MAXREGIONS;
	char      line[200];
	int       scan[8];
	
	
	int  readl          = atoi(argv[5]);
	int  tag            = atoi(argv[6]);
	bool strandspecific = (atoi(argv[7]) > 0) ? true : false;
	bool extendRetros   = (atoi(argv[8]) > 0) ? true : false;
	bool bestMatch      = (atoi(argv[9]) > 0) ? true : false;   /////-> adapt for new parameter ////////////////////////////////////////////////////

	// the nonredundant exons tag
	const int mt        = MAXTAG;
	int tags[mt];
	for (int z=0; z<mt; z++) tags[z]=0;
	
	Exon   **exons   = new Exon   *[exonLen];
	Read   **reads   = new Read   *[readLen];
	Region **regions = new Region *[regionLen];
	
	int   exonSize   = 0;
	int	  readSize   = 0;
	int	  regionSize = 0;
	int   exonMax    = 0; // the highest exon id (for iterators)
	int   geneMax    = 0; // the highest gene id (for iterators)
	
	//char     dummy[200];
	
	int lico = 0;
	int readCount = 0; // counts reads
	

	if (strandspecific) std::cerr << "         ### WARNING: USING STRANDS -> IN DEVELOPMENT ### " << std::endl << std::endl;
	
	
	  /*******************/
	 /* READ ANNOTATION */
	/*******************/

	map<int, int> exonsPerGene;
	int acceptedExons = 0;
	// READ EXON FILE
	std::ifstream exonFile(argv[1]);
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
			sscanf(line, "%d%d%d%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4], &scan[5], &scan[6], &scan[7]); // gene, exon, region, start, end, strand, core, uniscore
	
			// EXON PER GENE CHECK
			exonsPerGene[scan[0]]++;
			if (maxexon <= exonsPerGene[scan[0]]) {
				std::cerr << "PAGE FAULT: Raise MAXGENEEXONS!" << std::endl;
				std::cerr << scan[0] << "\t" << scan[1] << "\t" << scan[2] << "\t" << scan[3] << "\t" << scan[4] << "\t";
				std::cerr << scan[5] << "\t" << scan[6] << "\t" << scan[7] << std::endl;
				std::cerr << "(" << exonsPerGene[scan[0]] << " > " << maxexon << ")" << std::endl;
				for (map<int,int>::const_iterator it = exonsPerGene.begin(); it !=  exonsPerGene.end(); ++it) {
					std::cerr << it->first << "<>" << it->second << std::endl;
				}
				
				abort();
			}
			// MAXEXON CHECK
			if (scan[1] >= exonLen) {
				std::cerr << "PAGE FAULT: Raise MAXEXONS!" << std::endl;
				std::cerr << "(" << scan[1] << ")" << std::endl;
				abort();
			}
			// check tag
			tags[scan[6]]++;
			if (scan[6] >= tag) acceptedExons++;
			//create exon
			exons[scan[1]] = new Exon(scan);
			if (scan[1] > exonMax) exonMax = scan[1]; // update exonMax
			if (scan[0] > geneMax) geneMax = scan[0]; // update exonMax
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
				if (regionSize >= regionLen) {
					std::cerr << "PAGE FAULT: Raise MAXREGIONS!" << std::endl;
					std::cerr << "(" << exons[scan[1]]->region() << ")" << std::endl;
					exit(1);
				}
				regions[regionSize] = new Region(exons[scan[1]]->region());
				regions[regionSize]->addExonToGene(exons[scan[1]]);
				regionSize++;
			}
		} else {
			std::cerr << "Out of bounds (exonLen)!" << std::endl;
		}
		if (!exonFile.fail()) exonSize++; // count sucessfully read line -> limits array size
	}
	exonFile.close();
	// stats
	int untaggedGenes = 0;
	std::cerr << "\r" << "exonFile: " << exonSize << " records from " << regionSize << " contiguous regions (" << acceptedExons << " tagged)" << std::endl;
	for (int r = 0; r < regionSize; r++) untaggedGenes += regions[r]->orderGenesByEnd(tag); // writes quickaccess (ordered) array and check for tagged_lengths to be non null
	if (untaggedGenes > 0) std::cerr << std::endl << "WARNING: " << untaggedGenes << " genes will not be evaluated due to chosen analysis tag (They will have -1 Expression)" << std::endl << std::endl;

	// Exon tag stats
	for (int i = 0; i < mt; ++i) {
		if (tags[i]) {
			std::cerr << " TAG " << i << "\t" << tags[i] << " Exons";
			if (i == tag) std::cerr << " <=";
			else if (i > tag) std::cerr << " <-";
			std::cerr << std::endl;
		}
	}
	
	// region index
	typedef map<int, int> simpleset;
	simpleset regionids;
	
	delete [] reads;
	reads   = new Read *[readLen];
	for (int r = 0; r < readLen; r++) {
		if (reads[r] != NULL) {
			std::cerr << "WARNING: Reads array not NULL -> memory leak?" << std::endl;
			reads[r] = NULL; // resets reads
			//exit(1);
		}
	}
	
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
		if (readSize < readLen) {
			sscanf(line, "%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4]); // readid, strand, region, start, mismatches
			if (++lico % 1000 == 0) std::cerr << "\r" << lico;
			
			//set strand to 0 if not strand specific
			if (!(strandspecific)) scan[1] = 0;

			if (reads[scan[0]] == NULL) {
				reads[scan[0]] = new Read(scan[0]); // create read (allocate memory)
				if (readSize < scan[0]) readSize = scan[0]; 
				readCount++;
			}
			reads[scan[0]]->addFullMapping(scan[2], scan[1], scan[3], readl, scan[4]);
		} else {
			if (lico == 0) std::cerr << "Error: could not read from mappingFile " << argv[2] << std::endl;
			else std::cerr << "Out of bounds " << readSize << " (readLen)!" << std::endl;
			exit(1);
		}
	}
	std::cerr << "\r" << "mappingFile: " << lico << " lines" << std::endl;
	mappingFile.close();
	
	
	// READ SPLICE INDEX
	splicemap spm;
	std::ifstream spliceFile(argv[4]);
	lico = 0;
	//int dummy;
	std::cerr << "Reading Splice Index " << argv[4] << std::endl;
	while (!spliceFile.eof() && !spliceFile.fail()) {
		spliceFile.getline(line, 200);
		if (spliceFile.eof()) break;
		sscanf(line, "%d%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4], &scan[5]); //splice, (centersplice), seq_region, strand, start, end, rightside intron (0 if end)
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
		if (readSize < readLen) {
			sscanf(line, "%d%d%d%d%d", &scan[0], &scan[1], &scan[2], &scan[3], &scan[4]); // readid, strand, region, start, mismatches
			if (++lico % 1000 == 0) std::cerr << "\r" << lico;

			//set strand to 0 if not strand specific
			if (!(strandspecific)) scan[1] = 0;

			if (reads[scan[0]] == NULL) {
				reads[scan[0]] = new Read(scan[0]); // create read (allocate memory)
				if (readSize < scan[0]) readSize = scan[0]; 
				readCount++;
			}
			//cerr << "\r                    > " << scan[0] << "\t" << scan[1] << "\t" << scan[2] << "\t" << scan[3] << "\t" << scan[4];
			reads[scan[0]]->addSpliceMapping(spm.find(scan[2])->second, scan[1], scan[3], readl, scan[4]);
		} else {
			if (lico == 0) std::cerr << "Error: could not read from mappingFile " << argv[2] << std::endl;
			else std::cerr << "Out of bounds " << readSize << " (readLen)!" << std::endl;
			exit(1);
		}
	}
	std::cerr << "\r" << "splicemappingFile: " << lico << " lines" << std::endl;
	splicemappingFile.close();
	
	  /**********/
	 /* CHECKS */
	/**********/
	
	// check if region to add ( do we really have to add the region?
	int addedRegions = 0;
	for (simpleset::const_iterator it = regionids.begin(); it!=regionids.end(); it++) {
		bool added = false;
		for(int r=0; r < regionSize; r++) {
			if (regions[r]->region() == it->first) added = true;
		}
		if (added == false) {
			regions[regionSize] = new Region(it->first);
			regionSize++;
			addedRegions++;
		}
	}
	if (addedRegions > 0) std::cerr << "NOTE: Added " << addedRegions << " regions" << std::endl;
	
	// check if reads and mappings are same
	bool abortRun = false;
	for (int r = 0; r < readSize; r++) {
		if (reads[r] == NULL) continue;
		for (int m = 0; m < reads[r]->countMappings(); m++) {
			if (reads[r]->readMapping(m)->getRead()->rid() != reads[r]->rid()) {
				std::cerr << "OOPSIE " << reads[r]->readMapping(m)->getRead()->rid();
				std::cerr << " != " << reads[r]->rid() << std::endl;
				abortRun = true;
			}
		}
	}
	if (abortRun) exit(1);
	
	  /***************/
	 /* unify reads */
	/***************/
	std::cerr << std::endl << "Unifying read mappings..." << std::endl;
	time(&startTime);

	int map_c = 0, kill_c = 0;
	for(int r=0;r < readSize; ++r) {
		if (reads[r] == NULL) continue;
		if (r % 1000 == 0) {
			time(&currentTime);
			double dif = difftime (currentTime,startTime);
			double estim = (dif / r) * (readSize - r);
			if (dif > 5) std::cerr << "\r ~ " << (int)estim << " seconds remaining    ";
		} //TIMER
		map_c  += reads[r]->countMappings();
		kill_c += reads[r]->unify(bestMatch); ///////////CONTINUE HERE FOR STRANDSPECIFIC (??? necessary)
	}
	int acceptedMappings = map_c - kill_c;
	time(&currentTime);
	finaltime = difftime (currentTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed  (" << kill_c << " mappings out of " << map_c << " were removed)" << std::endl;

	
	  /************************************/
	 /* map to annotation and distribute */
	/************************************/
	
	// Map unambiguous reads (weight = 1 in any case) 
	std::cerr << std::endl << "Mapping unambiguous reads..." << std::endl;
	time(&startTime);
	int unambiguousReads = 0, ambiguousReads = 0;
	for (int j = 0; j < readSize; j++) {
		if (reads[j] == NULL) continue;
		if (reads[j]->rid() != j) {
			std::cerr << "ABORT TRAP: Memory leak?" << std::endl << j << "!=" << reads[j]->rid() << std::endl;
			abort();
		}
		if (reads[j]->isUnambiguous()) {
			++unambiguousReads;
			if (j % 1000 == 0) {
				time(&currentTime);
				double dif = difftime (currentTime,startTime);
				double estim = (dif / j) * (readSize - j);
				if (dif > 5) std::cerr << "\r ~ " << (int)estim << " seconds remaining    ";
			} //TIMER
			for (int m = 0; m < reads[j]->countMappings(); m++) {
				for (int r = 0; r < regionSize; r++) {
					if (reads[j]->readMapping(m)->regiomap(regions[r])) {
						regions[r]->projection(reads[j]->readMapping(m), readl);
					}
				}
			}
		}
		else ++ambiguousReads;
	}
	time(&currentTime);
	finaltime = difftime (currentTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed  (" << unambiguousReads << " unambiguous reads out of " << ambiguousReads + unambiguousReads <<")" << std::endl;
	
	// Setting of preliminary rpkms (per exons, independent of tag)
	std::cerr << std::endl << "Evaluating..." << std::endl;
	for (int r = 0; r < regionSize; r++) {
		if (regions[r] == NULL) continue;
		std::cerr << "(" << r << "|" << regions[r]->gCount() << ")";
		for (int g = 0; g < regions[r]->gCount(); g++) {
			regions[r]->getGene(g)->calcPrexpressPerExonUsingTag(readl, tag);
		}
	}
	std::cerr << std::endl;

	// Map all reads (ambiguous reads with redundancy) 
	std::cerr << std::endl << "Mapping ambiguous reads..." << std::endl;
	time(&startTime);
	for (int j = 0; j < readSize; j++) {
		if (reads[j] == NULL) continue;
		if (!(reads[j]->isUnambiguous())) { // also sets frac to 1 if unambiguous
			if (j % 1000 == 0) {
				time(&currentTime);
				double dif = difftime (currentTime,startTime);
				double estim = (dif / j) * (readSize - j);
				if (dif > 5) std::cerr << "\r ~ " << (int)estim << " seconds remaining    ";
			} //TIMER
			for (int m = 0; m < reads[j]->countMappings(); m++) {
				for (int r = 0; r < regionSize; r++) {
					if (reads[j]->readMapping(m)->regiomap(regions[r])) {
						regions[r]->projection(reads[j]->readMapping(m), readl);
					}
				}
			}
		}
	}
	time(&currentTime);
	finaltime = difftime (currentTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed  (" << ambiguousReads << " ambiguous reads out of " << ambiguousReads + unambiguousReads <<")" << std::endl;
	
	
	  /*###########******/
	 /* Create Islands */
	/******###########*/
	
	// ADD STRAND SPECIFIC ISLANDS
	// if not strandspecific set all islands to be direction 0 (mappings should already be 0 strand so set strand upon creation)
	
	
	// ISLAND CREATOR (takes data from projection routine [sorting out unassigned reads])
	std::cerr << std::endl << "Building and connecting islands..." << std::endl;
	int newExpressedRetrocopies = 0;
	boost::progress_display show_progress( regionSize );
	for (int r = 0; r < regionSize; r++) {
		++show_progress;
		regions[r]->buildIslands();
		regions[r]->connectIslands();
		regions[r]->verifyIslands(extendRetros); // if zero argument will assign NULLgene to all islands
		if (extendRetros) {
			int newRetrogenes = regions[r]->typeIslands(tag, readl); // sets Rflag if island is a suspected retrocopy
			newExpressedRetrocopies += newRetrogenes;
		}
	}
	std::cerr << std::endl;

	
	  /***************************/
	 /* Create Retro Annotation */
	/***************************/
	
	if (newExpressedRetrocopies > 0 && (extendRetros)) {
		int gMax = geneMax;
		int eMax = exonMax;
		int addedRetrogenes = 0;
		std::cerr << std::endl << "Completing Annotation and Remapping on new retros..." << std::endl;
		time(&startTime);
		for (int r = 0; r < regionSize; r++) {
			for (int i = 0; i < regions[r]->iCount(); i++) {
				if (regions[r]->getIsland(i)->getRflagchar() == 'R') {
					// island is an expressed retro
					Exon * nex = regions[r]->addIslandAnnotationAndResetReads(regions[r]->getIsland(i), ++gMax, ++eMax,tag, readl);
					std::cerr << "\r" << ++addedRetrogenes << "  -> RETROGENE " << nex->gid();
					//std::cerr << "(ISLAND_" << regions[r]->getIsland(i)->getIslandID();
					//std::cerr << "|PARENT_" << regions[r]->getIsland(i)->bestGeneID();
					//std::cerr << "|RETRO_"  << nex->gid() << ")";
				}
			}
			//reorder genes
			regions[r]->orderGenesByEnd(tag);
		}
		std::cerr << "\r                                         \r" << addedRetrogenes << " retrogenes were added" <<  std::endl;

		  /**********************/
		 /* Redistribute reads */
		/**********************/
		
		int remappedReads = 0;
		for (int r=0; r <= readSize; r++) {
			if (r % 1000 == 0) {
				time(&currentTime);
				double dif = difftime (currentTime,startTime);
				double estim = (dif / r) * (readSize - r);
				if (dif > 5) std::cerr << "\r ~ " << (int)estim << " seconds remaining    ";
			} //TIMER
			if (reads[r] == NULL) continue;
			if (reads[r]->flaggedForRemap() && !(reads[r]->isUnambiguous())) {
				++remappedReads;
				for (int m = 0; m < reads[r]->countMappings(); m++) {
					for (int g = 0; g < regionSize; g++) {
						if (reads[r]->readMapping(m)->regiomap(regions[g])) {
							regions[g]->projection(reads[r]->readMapping(m), readl);
						}
					}
				}
			}
		}
		std::cerr << "\r " << (int)finaltime << " seconds elapsed (Remapped " << remappedReads << " Reads)           " << std::endl;
	}
	
	  /*###########********/
	 /* Distribute Reads */
	/********###########*/
	
	//RAPARM ALGO
	std::cerr << std::endl << "Fractionating reads..." << std::endl;
	time(&startTime);
	for (int r = 0; r < readSize; r++) {
		// skips unassigned and unambiguous reads that may not have more than one gene assignment
		if (reads[r] == NULL || reads[r]->isUnambiguous()) continue;
		if (r % 1000 == 0) {
			time(&currentTime);
			double dif = difftime (currentTime,startTime);
			double estim = (dif / r) * (readSize - r);
			if (dif > 5) std::cerr << "\r ~ " << (int)estim << " seconds remaining    ";
		} //TIMER
		if (reads[r]->totalAssignments() > 0) { // unasss
			reads[r]->fractionate(); // update weightfrac for all assignments
		}
	}
	time(&currentTime);
	finaltime = difftime (currentTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;

	  /*###########****#**/
	 /* Calc Expression */
	/******############*/
	
	
	//RPKM CALC
	int readIsMapped        = 0;
	bool hasGeneMapping     = false;
	int readIsTagMapped     = 0;
	bool hasTaggedMapping	= false;
	int intronFlagCount     = 0;
	int exonFlagCount       = 0;
	int intergenicFlagCount = 0;
	std::cerr << std::endl << "Evaluating Gene Expression..." << std::endl;
	time(&startTime);
	for (int r=0; r <= readSize; r++) {
		if (r % 1000 == 0) {
			time(&currentTime);
			double dif = difftime (currentTime,startTime);
			double estim = (dif / r) * (readSize - r);
			if (dif > 5) std::cerr << "\r ~ " << (int)estim << " seconds remaining    ";
		} //TIMER
		if (reads[r] == NULL) continue;
		hasGeneMapping = false;
		hasTaggedMapping = false;
		for (int m=0; m < reads[r]->countMappings(); m++) {
			//flag counting
			if (reads[r]->readMapping(m)->getFlag() == -1) intergenicFlagCount++;
			if (reads[r]->readMapping(m)->getFlag() == 0)  intronFlagCount++;
			if (reads[r]->readMapping(m)->getFlag() == 1)  exonFlagCount++;
			// add readfracs to gene
			if (reads[r]->readMapping(m)->assignments() > 0) {
				hasGeneMapping = true;
				for(int a=0; a < reads[r]->readMapping(m)->assignments(); a++) {
					if (reads[r]->readMapping(m)->getAssignment(a)->commitExpression(tag)) {
						hasTaggedMapping = true;
					} // writes mapscore to gene
				}
			}
		}
		if (hasGeneMapping) readIsMapped++;
		if (hasTaggedMapping) readIsTagMapped++;
	}
	time(&currentTime);
	finaltime = difftime (currentTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;
	
	 
	  /*###########*****/
	 /* Calc Coverage */
	/******##########*/
	
	std::cerr << std::endl << "Evaluating Coverage..." << std::endl;
	time(&startTime);
	
	int debugsum = 0;
	boost::progress_display cov_progress( readSize );
	for (int r=0; r <= readSize; r++) {
		++cov_progress;
		if (reads[r] == NULL) continue;
		if (reads[r]->calculateCoverage(regions, regionSize)) debugsum++;
	}
	time(&currentTime);
	finaltime = difftime (currentTime,startTime);
	std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;
	if (debugsum != readCount) std::cerr << "WARNING: Read " << readCount << " != debugsum " << debugsum << std::endl;
	

	  /**********/
	 /* OUTPUT */
	/**********/

	//OUTPUT RAW
	int geneAssignCount = 0;
	int noAssignCount = 0;
	if (argv[10][0] != '-') {
		std::cerr << std::endl << "Writing mapping protocol..." << std::endl;
		time(&startTime);
		std::ofstream outFile(argv[10]);
		boost::progress_display show_progress( readSize );
		for (int r=0; r <= readSize; r++) {
			++show_progress;
			if (reads[r] == NULL) continue;
			// count mapping sucess for each read
			if (reads[r]->totalAssignments() > 0) geneAssignCount++;
			else noAssignCount++;
			for (int m=0; m < reads[r]->countMappings(); m++) {
				if (reads[r]->readMapping(m)->assignments() > 0) {
					for(int a=0; a < reads[r]->readMapping(m)->assignments(); a++) {
						// print one line for each assignment
						outFile << reads[r]->rid()                                                 << "\t";
						outFile << reads[r]->readMapping(m)->getFlag()                             << "\t";
						outFile << reads[r]->readMapping(m)->getAssignment(a)->getGeneID()         << "\t";
						outFile << reads[r]->readMapping(m)->getAssignment(a)->getExonID()         << "\t";
						outFile << reads[r]->readMapping(m)->getAssignment(a)->getOverlapFrac(tag) << "\t";
						outFile << reads[r]->readMapping(m)->getAssignment(a)->getExpressFrac()    << "\t";
						outFile << reads[r]->readMapping(m)->region()                              << "\t";
						outFile << reads[r]->readMapping(m)->str()                                 << "\t";
						outFile << reads[r]->readMapping(m)->cigarline()                           << std::endl;
					}
				}
				else {
					// print a non-gene map
					outFile << reads[r]->rid()                       << "\t";
					outFile << reads[r]->readMapping(m)->getFlag()   << "\t";
					outFile << 0                                     << "\t";
					outFile << 0                                     << "\t";
					outFile << 0.0                                   << "\t";
					outFile << 0.0                                   << "\t";				
					outFile << reads[r]->readMapping(m)->region()    << "\t";
					outFile << reads[r]->readMapping(m)->str()       << "\t";
					outFile << reads[r]->readMapping(m)->cigarline() << std::endl;
				}
			}
		}
		outFile.close();
		time(&currentTime);
		finaltime = difftime (currentTime,startTime);
		std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;
	}
	else {//just count the assignments
		for (int r=0; r <= readSize; r++) {
			if (reads[r] == NULL) continue;
			if (reads[r]->totalAssignments() > 0) geneAssignCount++;
			else noAssignCount++;
		}
	}
	
	// Writing RPKM values
	if (argv[11][0] != '-') {
		std::cerr << std::endl << "Writing RPKM values..." << std::endl;
		time(&startTime);
		std::ofstream rpkmOut(argv[11]);
		boost::progress_display show_progress( regionSize );
		for (int r = 0; r < regionSize; r++) {
			++show_progress;
			if (regions[r] == NULL || regions[r]->gCount() == 0) continue;
			for (int g = 0; g < regions[r]->gCount(); g++) {
				//printer
				rpkmOut << regions[r]->getGene(g)->gid() << "\t";
				rpkmOut << ((regions[r]->getGene(g)->gid() > geneMax) ? "NOVEL" : "BASE") << "\t";
				rpkmOut << regions[r]->region() << "\t";
				rpkmOut << regions[r]->getGene(g)->location() << "\t";
				rpkmOut << regions[r]->getGene(g)->tagged_length(tag) << "\t";
				rpkmOut << regions[r]->getGene(g)->mappedReads() << "\t";
				rpkmOut << regions[r]->getGene(g)->rpkm(readIsTagMapped, tag) << std::endl;
			}
		}
		rpkmOut.close();
		time(&currentTime);
		finaltime = difftime (currentTime,startTime);
		std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;
	}
	
	//Writing intergenic/intronic islands
	double exonicL     = 0;
	double exonicR     = 0;
	double intronicL   = 0;
	double intronicR   = 0;
	double intergenicL = 0;
	double intergenicR = 0;
	if (argv[12][0] != '-') {
		std::cerr << std::endl << "Writing Islands..." << std::endl;
		time(&startTime);
		std::ofstream islandOut(argv[12]);
		boost::progress_display show_progress( regionSize );
		for (int r = 0; r < regionSize; r++) {
			++show_progress;
			if (regions[r] == NULL || regions[r]->iCount() == 0) continue;
			for (int i = 0; i < regions[r]->iCount(); i++) {
				// coverage calc
				if      (regions[r]->getIsland(i)->getFlag() == -1) {
					intergenicL += regions[r]->getIsland(i)->length();
					intergenicR += regions[r]->getIsland(i)->readCount();
				}
				else if (regions[r]->getIsland(i)->getFlag() ==  0) {
					intronicL += regions[r]->getIsland(i)->length();
					intronicR += regions[r]->getIsland(i)->readCount();
				}
				else if (regions[r]->getIsland(i)->getFlag() ==  1) {
					exonicL += regions[r]->getIsland(i)->length();
					exonicR += regions[r]->getIsland(i)->readCount();
				}
				
				//printing
				islandOut << regions[r]->getIsland(i)->GFFline(readl) << std::endl;
			}
		}
		islandOut.close();
		time(&currentTime);
		finaltime = difftime (currentTime,startTime);
		std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;
		
	}
	else { // just calc coverage
		for (int r = 0; r < regionSize; r++) {
			for (int i = 0; i < regions[r]->iCount(); i++) {
				// coverage calc
				if      (regions[r]->getIsland(i)->getFlag() == -1) {
					intergenicL += regions[r]->getIsland(i)->length();
					intergenicR += regions[r]->getIsland(i)->readCount();
				}
				else if (regions[r]->getIsland(i)->getFlag() ==  0) {
					intronicL += regions[r]->getIsland(i)->length();
					intronicR += regions[r]->getIsland(i)->readCount();
				}
				else if (regions[r]->getIsland(i)->getFlag() ==  1) {
					exonicL += regions[r]->getIsland(i)->length();
					exonicR += regions[r]->getIsland(i)->readCount();
				}
			}
		}
	}
	
	
	//Writing Coverage streams
	if (argv[13][0] != '-') {
		std::cerr << std::endl << "Writing Coverage..." << std::endl;
		time(&startTime);
		std::ofstream coverageOut(argv[13]);
		boost::progress_display show_progress( regionSize );
		for (int r = 0; r < regionSize; r++) {
			++show_progress;
			if (regions[r] == NULL || regions[r]->cov.size() == 0) continue;
			// run through cov map
			int coco = 0;
			int laststart = 1;
			for(covmap::const_iterator it = regions[r]->cov.begin(); it != regions[r]->cov.end(); ++it){
				// calc & print
				if (it->second != 0) { // change can become zero if reads' kiss
					coverageOut << regions[r]->region() << "\t";
					coverageOut << ((strandspecific) ? "1" : "0" ) << "\t";
					coverageOut << laststart << "\t";
					coverageOut << (it->first - 1) << "\t";
					coverageOut << (it->first - laststart) << "\t";
					coverageOut << coco << std::endl;
					// update
					laststart = it->first;
					coco += it->second;
				}
			}
			if (strandspecific) {
				coco = 0;
				laststart = 1;
				for(covmap::const_iterator it = regions[r]->revcov.begin(); it != regions[r]->revcov.end(); ++it){
					// calc & print
					if (it->second != 0) { // change can become zero if reads' kiss
						coverageOut << regions[r]->region() << "\t";
						coverageOut << "-1" << "\t";
						coverageOut << laststart << "\t";
						coverageOut << (it->first - 1) << "\t";
						coverageOut << (it->first - laststart) << "\t";
						coverageOut << coco << std::endl;
						// update
						laststart = it->first;
						coco += it->second;
					}
				}
			}
		}
		coverageOut.close();
		time(&currentTime);
		finaltime = difftime (currentTime,startTime);
		std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;
	}

	if (argv[14][0] != '-') {
		std::cerr << std::endl << "Writing Unique Coverage..." << std::endl;
		time(&startTime);
		std::ofstream uniqueOut(argv[14]);
		boost::progress_display show_progress( regionSize );
		for (int r = 0; r < regionSize; r++) {
			++show_progress;
			if (regions[r] == NULL || regions[r]->unicov.size() == 0) continue;
			// run through cov map
			int coco = 0;
			int laststart = 1;
			for(covmap::const_iterator it = regions[r]->unicov.begin(); it != regions[r]->unicov.end(); ++it){
				// calc & print
				if (it->second != 0) { // change can become zero if reads' kiss
					uniqueOut << regions[r]->region() << "\t";
					uniqueOut << ((strandspecific) ? "1" : "0" ) << "\t";
					uniqueOut << laststart << "\t";
					uniqueOut << (it->first - 1) << "\t";
					uniqueOut << (it->first - laststart) << "\t";
					uniqueOut << coco << std::endl;
					// update
					laststart = it->first;
					coco += it->second;
				}
			}
			if (strandspecific) {
				coco = 0;
				laststart = 1;
				for(covmap::const_iterator it = regions[r]->revunicov.begin(); it != regions[r]->revunicov.end(); ++it){
					// calc & print
					if (it->second != 0) { // change can become zero if reads' kiss
						uniqueOut << regions[r]->region() << "\t";
						uniqueOut << "-1" << "\t";
						uniqueOut << laststart << "\t";
						uniqueOut << (it->first - 1) << "\t";
						uniqueOut << (it->first - laststart) << "\t";
						uniqueOut << coco << std::endl;
						// update
						laststart = it->first;
						coco += it->second;
					}
				}
			}				
		}
		uniqueOut.close();
		time(&currentTime);
		finaltime = difftime (currentTime,startTime);
		std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;
	}
	
	
	// NOT ADAPTED TO STRAND SPECIFIC!!!! (necessary?)
	if (argv[15][0] != '-') {
		std::cerr << std::endl << "Calculating and Writing Splice Frequencies..." << std::endl;
		time(&startTime);
		intdoublemap splicefreq;
		for (int r=0; r <= readSize; r++) {
			if (r % 1000 == 0) {
				time(&currentTime);
				double dif = difftime (currentTime,startTime);
				double estim = (dif / r) * (readSize - r);
				if (dif > 5) std::cerr << "\r ~ " << (int)estim << " seconds remaining    ";
			} //TIMER
			if (reads[r] == NULL) continue;
			for (int m=0; m < reads[r]->countMappings(); m++) {
				if (reads[r]->readMapping(m)->isUnspliced()) continue;
				double express = reads[r]->readMapping(m)->highestExpressFrac();
				if (sameStrand(reads[r]->readMapping(m)->str(), reads[r]->readMapping(m)->getSplice()->str())) {
					intvec * ii = reads[r]->readMapping(m)->getSplice()->spannedIntrons( reads[r]->readMapping(m)->getSlices() );
					// add intron ids
					for (unsigned int t = 0; t != ii->size(); ++t) splicefreq[ ii->at(t) ] += express;
					delete ii;
				}
			}
		}
		time(&currentTime);
		finaltime = difftime (currentTime,startTime);
		std::cerr << "\r"; for(int x=0; x < 150; x++) std::cerr << " ";
		std::cerr << "\r " << (int)finaltime << " seconds elapsed               " << std::endl;
		// write
		std::ofstream spliceOut(argv[15]);
		for (intdoublemap::const_iterator sf = splicefreq.begin(); sf != splicefreq.end(); ++sf) {
			spliceOut << sf->first << "\t" << sf->second << std::endl;
		}
		spliceOut.close();
	}
	
	
	// some basic stats and say goodbye
	time(&endTime);
	finaltime = difftime(endTime,beginTime);

	std::cerr << std::endl << "=====================================================================";
	std::cerr << std::endl << "========================== RAPARM REPORT ============================";
	std::cerr << std::endl << "=====================================================================" << std::endl;
	std::cerr << std::endl << "\tAnalyzed " << acceptedMappings << " mappings in " << finaltime << " seconds";
	std::cerr << " (" << (int) (36*acceptedMappings/(finaltime*10000)) << " Mm/h)" << std::endl;
	std::cerr << std::endl << "\tMapping Stats\t" << intergenicFlagCount << "\t Intergenic";
	std::cerr << std::endl << "\t\t\t" << intronFlagCount << "\t Intronic";
	std::cerr << std::endl << "\t\t\t" << exonFlagCount << "\t Exonic" << std::endl;
	std::cerr << std::endl << "\tRead Stats\t" << geneAssignCount << "\t Assigned";
	std::cerr << std::endl << "\t\t\t" << noAssignCount << "\t Unassigned" << std::endl;
	std::cerr << std::endl << "\t\t\t"; std::cout << readCount << "\t" << " Total Reads" << std::endl;
	std::cerr << std::endl << "\t" << newExpressedRetrocopies << " new expressed Retrocopies" << std::endl;
	std::cerr << std::endl << "\t";
	std::cout << (intergenicR * readl / intergenicL) << "\t"; 
	std::cerr << " fold intergenic coverage";
	std::cerr << std::endl << "\t";
	std::cout << (exonicR * readl / exonicL) << "\t";
	std::cerr << " fold exonic coverage";
	std::cerr << std::endl << "\t";
	std::cout << (intronicR * readl / intronicL);
	std::cerr << "\t fold intronic coverage" << std::endl;
	
	std::cerr << std::endl << "                    * * * HAVE A NICE DAY * * *                      " << std::endl << std::endl;
	std::cout << std::endl;

	exit(0);	// SUCCESS
} // end of main

