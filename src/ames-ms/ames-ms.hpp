/*
 *  ames.hpp
 *  AMbiguity EStimator (for exon annotations)
 *
 *  Created by David Brawand on 14.10.09.
 *  MAJOR REVISION FOR MULTISPLICE 30.03.10
 *  Copyright 2009 DavidBrawand/CIG. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include "universe.hpp"

using namespace std ;

bool Read::unify() {
	// remove identical (may accelerate graph step afterwards)
	for (unsigned int i=0; i < mappings.size(); i++) {
		for (unsigned int j=i; j < mappings.size(); j++) {
			if (i == j) continue; 
			if (mappings[i]->operator==(*mappings[j])) {
				// remove j
				mappings.erase(mappings.begin()+j);
				--j; // step one back as vector gets shortenend by previous operation
			}
		}
	}
	
	// make sure that it will not be evaluated if there is just one mapping ->speedup by using isUnique()
	using namespace boost ;
	typedef adjacency_list<vecS, vecS, undirectedS> Graph;
	
	// declare a graph object
    Graph G;
	
	for (unsigned int i=0; i < mappings.size(); i++) {
		for (unsigned int j=i; j < mappings.size(); j++) {
			if (mappings[i]->overlaps(mappings[j])) {
				add_edge(i, j, G);// add edge
			}
		}
	}
	// find connected components
	std::vector<int> component(num_vertices(G));
    int num = connected_components(G, &component[0]);
    
	// if more than 1 -> ambiguous mappings -> drop
	if (num != 1) return false;
	else {
		//select among highest splice class
		// determine highest splice class
		intmap spliceclasses;
		for (unsigned int mp = 0; mp != mappings.size(); mp++) {
			spliceclasses[mappings[mp]->gaps()]++;
		}
		intmap::iterator sc = spliceclasses.end();
		--sc; // get hightest spliceclass iterator
		int selected = rand() % sc->second;
		int inc = 0; // the spliceclass counter
		for(unsigned int m = 0; m != mappings.size(); ++m) {
			if (mappings[m]->gaps() == sc->first && selected == inc++) {
				mappings[m]->accept();
				return true;
			}
		}
		//////////////////////////// KEEP THIS FOR FUTURE USE IN RAPARM-MS ///////////////////////////////////
		//std::vector<int>::size_type i;
		//cout << "Total number of components: " << num << endl;
		//for (i = 0; i != component.size(); ++i) {
		//	cout << "Vertex " << i <<" is in component " << component[i] << endl;
		//}
		//////////////////////////// KEEP THIS FOR FUTURE USE IN RAPARM-MS ///////////////////////////////////
	}
	abort(); //a mapping should have been selected at this point
}

//OK
bool Exon::assignOverlap(Mapping* m, int readl) {
	if (this->region() != m->region()) {
		std::cerr << "ERROR: RegionError (" << this->region() << " <> " << m->region() << ")" << std::endl;
		this->show();
		m->show();
		abort();
	}
	//m->show();
	bool added = false;
	for(intmap::const_iterator sl = m->slices->begin(); sl != m->slices->end(); ++sl) {
		// both in
		if ( this->start() <= sl->first && sl->second <= this->end() ) {
			int coded = readl - (sl->second - sl->first + 1);
			//std::cerr << "{B:" << coded << "}";
			this->partmap[coded]++;
			added = true;
		}
		// both out
		else if ( sl->first < this->start() && this->end() < sl->second ) {
			int coded = readl - (this->end() - this->start());
			//std::cerr << "{S:" << coded << "}";
			this->partmap[coded]++;
			added = true;
		}
		// left in
		else if ( this->start() <= sl->first && sl->first <= this->end() && sl->second > this->end() ) {
			int coded = readl - (this->end() - sl->first + 1);
			//std::cerr << "{L:" << coded << "}";
			this->partmap[coded]++;
			added = true;
		}
		// right in
		else if ( this->start() <= sl->second && sl->second <= this->end() && sl->first < this->start() ) {
			int coded = readl - (sl->second - this->start() + 1);
			//std::cerr << "{R:" << coded << "}";
			this->partmap[coded]++;
			added = true;
		}
		else {
			// no overlap (both ends on either side)
		}
	}
	return added;
}

//OK
bool Region::projection(Mapping * mp, int readl) { // requires ordered array of genes? uses dichotomy (oneway)
	if (rid != mp->region()) {
		std::cerr << "ERROR: Cannot project mapping with region " << mp->region() << " on region " << rid << std::endl;
		abort();
	}
	int isMappedToExon = 0;
	for(int i = dicoindex(mp); i < gsize; i++) {
		if (i >= gsize) {
			std::cerr << "FATAL: dicosearch failure" << std::endl;
			abort();
		}
		if (garray[i]->region() != rid) {
			std::cerr << "FATAL: gene region association error" << std::endl;
			abort();
		}
		if (mp->end() < garray[i]->end - maxGeneLength) break; // stops if no following genes may contain readmapping
		if (mp->overlaps(garray[i])) {
			for (int j = 0; j < garray[i]->exonCount; j++) {				
				if (garray[i]->getExon(j)->assignOverlap(mp, readl)) isMappedToExon++;
			}
		}
	}
	if (isMappedToExon > 0) return true;
	return false;
}

//OK
bool Region::spliceProjection(int a, int z) { // return true if success
	int isMappedToExon = 0;
	a -= 1;
	z += 1;
	for(int i = dicoindex(a); i < gsize; i++) {
		if (i >= gsize) {
			std::cerr << "FATAL: dicosearch failure" << std::endl;
			abort();
		}
		if (garray[i]->region() != rid) {
			std::cerr << "FATAL: gene region association error" << std::endl;
			abort();
		}
		if (z < garray[i]->end - maxGeneLength) break; // stops if no following genes may contain readmapping
		for (int j = 0; j < garray[i]->exonCount; j++) {
			if (garray[i]->getExon(j)->end() >= a && garray[i]->getExon(j)->start() <= a) {
				isMappedToExon++;
				garray[i]->getExon(j)->addAnchor(z);
				garray[i]->getExon(j)->rightanchor = 1;
			}
			if (garray[i]->getExon(j)->end() >= z && garray[i]->getExon(j)->start() <= z) {
				isMappedToExon++;
				garray[i]->getExon(j)->addAnchor(a);
				garray[i]->getExon(j)->leftanchor = 1;
			}
		}
	}
	if (isMappedToExon > 0) return true;
	return false;
}

//NEW UNIESTIM
int Exon::uniqEstimate(int readl) {
	// use anchors size to weigh spliced uniqueness
	float afrac = (anchors.size() == 0) ? 0.0 : ((leftanchor + rightanchor) / (float) anchors.size());
	float score = 0.0;
	for (int i = 0 ; i < readl ; ++i) {
		float redufrac = (i == 0) ? 1.0 : afrac;
		if (this->partmap[i] == 0) continue;
		float f = (float)(readl - i) / readl;
		float w = (float)(this->partmap[i] * 0.5 * redufrac); // as reverse-complement is also included
		float a = w * f;
		score += a;
	}
	this->setUniscore(score);
	return static_cast<int>(score + 0.5);
}


//OLD (no splice equilibrate)
int Exon::uniqEstimateBrute(int readl) {
	float score = 0.0;
	for (int i = 0 ; i < readl ; ++i) {
		if (this->partmap[i] == 0) continue;
		float f = (float)(readl - i) / readl;
		float w = (float)(this->partmap[i] * 0.5); // as reverse-complement is also included
		float a = w * f;
		score += a;
	}
	this->setUniscore(score);
	return static_cast<int>(score + 0.5);
}
