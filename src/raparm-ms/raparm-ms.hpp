/*
 *  raparm-ms.hpp
 *  raparm-ms (The ultimate Retro-Aware, Probabilistic Ambiguous Read Mapping program, multisplice)
 *
 *  Created by David Brawand on 11.03.10.
 *  Copyright 2010 Université de Lausanne. All rights reserved.
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
#include <set>
#include <time.h>
#include "universe.hpp"


using namespace std ;

/* ------------------------------------------------------------- */
/* TEMPLATES                                                     */
/* ------------------------------------------------------------- */
template <class T> class Matrix {
	T **data;
	unsigned int x, y;
public:
	Matrix(unsigned int w, unsigned int h) {
		x = w;
		y = h;
		data = new T *[w];
		for(unsigned int a=0; a<w; a++) {
			data[a] = new T[h];
		}
	}
	~Matrix(void) {
		for(unsigned int a=0; a<x; a++) {
			delete [] data[a];
		} 
		delete [] data;
	}
	inline T *operator [] (unsigned int a) {
		return data[a];
	}
};
/* ------------------------------------------------------------- */
/* GLOBAL TYPEDEFS AND STATICS                                   */
/* ------------------------------------------------------------- */
const static int maxexon = MAXGENEEXONS;
const static int maxtag  = MAXTAG;
/* ------------------------------------------------------------- */
/* Inline Functions                                              */
/* ------------------------------------------------------------- */
inline double closed_interval_rand(double x0, double x1) {
	return x0 + (x1 - x0) * rand() / ((double) RAND_MAX);
}

/* ------------------------------------------------------------- */
/* Function declarations (to avoid forward declaration)          */
/* ------------------------------------------------------------- */


int Gene::exonAddOverlaps(Exon * x) {
	// returns how may exons in set overlap exon x
	int ovp = 0;
	// count overlaps
	if (exonCount > 0) {
		for (int ec = 0; ec < exonCount; ec++) {
			if (exons[ec]->tags() == x->tags()) {
				if (!((exons[ec]->end() < x->start()) || (exons[ec]->start() > x->end()))) {
					ovp++;
				}
			}
		}
	}
	// add exon
	x->setGene(this);
	exons[exonCount] = x;
	exonCount++;
	if (exonCount > maxexon) std::cerr << "OVERFLOW: " << maxexon << " " << exonCount << std::cerr;
	//return overlap number 
	return ovp;
}

//OK (eliminates redundancy compared to simpleAddExon)
void Gene::addExon (Exon * ex) {
	if (start > end) {
		start = ex->start();
		end   = ex->end();
	}
	if (ex->start() < start) start = ex->start();
	if (ex->end()   > end)   end   = ex->end();
	
	if (int exonOverlaps = this->exonAddOverlaps(ex)) { // merge exons that overlap
		int oldExonCount = exonCount;
		// write map sorted by begin
		typedef multimap<int, Exon*> em; 
		em * allexons = new em[maxtag];
		for (int e = 0; e < exonCount ; e++) {
			allexons[exons[e]->tags()].insert(pair<int,Exon*>(exons[e]->start(),exons[e]));
		}
		
		// delete array and reallocate
		delete exons;
		exonCount = 0;
		Exon ** exons = new Exon *[maxexon];
		
		// overlap for each tag
		for (int t = 0; t < maxtag; t++) {
			if (allexons[t].size() == 0) continue;
			Exon * previous = NULL;
			// flatadd (keep existing, avoids calling setGene)
			for(em::const_iterator it = allexons[t].begin(); it != allexons[t].end(); ++it){
				// returns ordered by start
				if (previous != NULL) {
					if ((previous->tags() == it->second->tags()) && !(((previous->end() < it->second->start()) || (previous->start() > it->second->end())))) {
						//merge the two exons into previous
						if (previous->end() < it->second->end()) previous->end(it->second->end());
					}
					else {
						// no overlap write previous
						exons[exonCount] = previous;
						exonCount++;
						// set new previous
						previous = it->second;
					}
					
				}
				else {
					previous = it->second; // initial definition
				}
			}
			// write previous to array
			exons[exonCount] = previous;
			exonCount++;
			
		}
		// verify if newExonCount is corecct (according to overlap count)
		if (oldExonCount - exonOverlaps != exonCount) {
			std::cerr << std::endl << "FATAL: Gene::addExon hit abort trap" << std::endl;
			std::cerr << "oldExonCount " << oldExonCount << std::endl;
			std::cerr << "newExonCount " << exonCount << std::endl;
			std::cerr << "exonOverlaps " << exonOverlaps << std::endl;
			abort();
		}
		delete [] allexons;
	}
	return;
}

//CHECK
Gene * Mapping::assignedGene() {
	Gene * gr;
	int gid = 0;
	for (unsigned int a = 0; a < assign.size(); a++) {
		if (gid == 0) {
			gid = assign[a]->getGeneID();
			gr = assign[a]->getGene();
		}
		if (gid != assign[a]->getGeneID()) std::cerr << "\t\tconceptual error (multiple gene assignments in Mapping::assigedGene call" << std::endl;
	}
	return gr;
}

/***********************/
/* CHECKED SUBROUTINES */
/***********************/

//OK//STSP
int Read::unify(bool leastMismatch) {
	int removedMappings = 0;
	
	// kill identicial
	for (unsigned int i=0; i < mappings.size(); i++) {
		for (unsigned int j=i; j < mappings.size(); j++) {
			if (i == j) continue; // no self comparison
			if (mappings[i]->operator==(mappings[j])) {
				mappings.erase(mappings.begin() + j);
				removedMappings++;
				--j; // step one back as vector gets shortenend by previous operation
			}
		}
	}
	if (mappings.size() == 1) return removedMappings;
	
	// filter for least mismatch number
	if (leastMismatch) {
		// FIND LOWEST MISMATCH NUMBER
		int lowest = 1000; // so high it must be lowered
		for (unsigned int i=0; i < mappings.size(); i++) if (lowest > mappings[i]->mism) lowest = mappings[i]->mism;  
		// FILTER ALL MAPPINGS WITH MORE MISMATCHES
		for (unsigned int i=0; i < mappings.size(); i++) {
			if (mappings[i]->mism > lowest) {
				mappings.erase(mappings.begin() + i);
				removedMappings++;
				--i; // step one back as vector gets shortenend by previous operation
			}
		}
		if (mappings.size() == 1) return removedMappings;
	}
	
	// make sure that it will not be evaluated if there is just one mapping ->speedup by using isUnique()
	using namespace boost;
	typedef adjacency_list<vecS, vecS, undirectedS> Graph;
	
	// declare a graph object
    Graph G;
	for (unsigned int i=0; i < mappings.size(); i++) {
		for (unsigned int j=i; j < mappings.size(); j++) {
			if (sameStrand(mappings[i]->strand, mappings[j]->strand) && mappings[i]->overlaps(mappings[j])) {
				add_edge(i, j, G);// add edge
			}
		}
	}
	
	// find connected components
	std::vector<int> component(num_vertices(G));
	unsigned int num = connected_components(G, &component[0]); 
	if (num == mappings.size()) {
		// all vertices have no connection -> no cleanup needed -> accept all
		for (unsigned int i=0; i < mappings.size(); i++) mappings[i]->accept();
		return removedMappings;
	}
	else {
		imap compo;
		std::vector<int>::size_type i;
		for (i = 0; i != component.size(); ++i) compo.insert( pair<int, Mapping*>(component[i], mappings[i]) ); //<componentnumber,mappings>
		
		// foreach component decide which ones to kill (one will be kept)
		//mapset killer;
		for (unsigned int c = 0; c < num; ++c) { // foreach component
			pair<imap::iterator,imap::iterator> ret = compo.equal_range(c);
			// find highest splice class and its size
			intmap scc;
			for (imap::const_iterator cc = ret.first; cc != ret.second; ++cc) scc[cc->second->gaps()]++;
			intmap::iterator hc = scc.end();
			--hc;
			int keeper  = rand() % hc->second;
			int counter = 0;
			for (imap::const_iterator cc = ret.first; cc != ret.second; ++cc) {
				if (cc->second->gaps() == hc->first && keeper == counter++) cc->second->accept();
			}
		}

		
		// kill everything which has not been accepted
		for (unsigned int i = 0; i < mappings.size(); ++i) {
			if ( mappings[i]->notAccepted() ) {
				mappings.erase(mappings.begin() + i);
				--i;
				removedMappings++;
			}
		}
		// check if as many components and mappings
		Graph H;
		for (unsigned int i=0; i < mappings.size(); i++) {
			for (unsigned int j=i; j < mappings.size(); j++) {
				if (sameStrand(mappings[i]->strand, mappings[j]->strand) && mappings[i]->overlaps(mappings[j])) add_edge(i, j, H);// add edge
			}
		}
		std::vector<int> component(num_vertices(H));
		unsigned int connex = connected_components(H, &component[0]); 
		if (connex != mappings.size() || num != connex) {
			cout << "OOPSIE:" << endl << "Total number of components: " << num << endl << "Total number of mappings:   " << mappings.size() << endl;
		}
	}
	return removedMappings;
}

//OK
void Region::projection(Mapping * mp, int readl) { // requires ordered array of genes? uses dichotomy (oneway)
	if (rid != mp->region()) {
		std::cerr << "ERROR: Cannot project mapping with region " << mp->region() << " on region " << rid << std::endl;
		abort();
	}
	mp->setFlag(-1);
	for(int i = dicoindex(mp); i < gsize; i++) {
		if (mp->end() < garray[i]->end - maxGeneLength) break; // stops if no following genes may contain readmapping
		if (mp->overlaps(garray[i])) {
			mp->setFlag(0); // won't downgrade
			for (int j = 0; j < garray[i]->exonCount; j++) {
				if (garray[i]->getExon(j)->assignOverlap(mp, readl)) mp->setFlag(1);
			}
		}
	}
	
	// add slices for island build
	for(intmap::const_iterator it = mp->slices->begin(); it != mp->slices->end(); ++it) {
		allMaps.insert( pair<int, Mapping*>(it->first, mp) );
	}
	return;
}

//OK//STSP
bool Exon::assignOverlap(Mapping* m, int readl) {
	if (this->region() != m->region()) {
		std::cerr << "ERROR: RegionError (" << this->region() << " <> " << m->region() << ")" << std::endl;
		abort();
	}
	// simple optimization (shortens overlap evaluation in most cases and avoid separate overlap check function call as in previous versions of RAPARM)
	if (this->end() < m->start() || m->end() < this->start()) return false;
	
	if (!(sameStrand(strand, m->str()))) return false;
	
	bool added = false;
	for(intmap::const_iterator sl = m->slices->begin(); sl != m->slices->end(); ++sl) {
		int coded = -1;
		// both in (whole slice is within exon)
		if ( this->start() <= sl->first && sl->second <= this->end() ) {
			coded = readl - (sl->second - sl->first + 1);
			//std::cerr << "{B:" << coded << "}";
			this->partmap[coded]++;
			added = true;
		}
		// both out
		else if ( sl->first < this->start() && this->end() < sl->second ) {
			coded = readl - (this->end() - this->start());
			//std::cerr << "{S:" << coded << "}";
			this->partmap[coded]++;
			added = true;
		}
		// left in
		else if ( this->start() <= sl->first && sl->first <= this->end() && sl->second > this->end() ) {
			coded = readl - (this->end() - sl->first + 1);
			//std::cerr << "{L:" << coded << "}";
			this->partmap[coded]++;
			added = true;
		}
		// right in
		else if ( this->start() <= sl->second && sl->second <= this->end() && sl->first < this->start() ) {
			coded = readl - (sl->second - this->start() + 1);
			//std::cerr << "{R:" << coded << "}";
			this->partmap[coded]++;
			added = true;
		}
		else {
			// no overlap (both ends on either side)
		}
		
		// create assignment if coded changed
		if (coded >= 0) { // sucessfully mapped (calculate overlapFrac)
			double ovpf = (double)(readl - coded) / (double)readl;
			double expf = ((m->getRead()->isUnambiguous()) ? 1.0 : -1.0 );
			m->addAssignment(this->getGene(), this, ovpf, expf);
		}
	}
	return added;
}

//OK
void Mapping::addAssignment(Gene * g, Exon * e, double ovp, double ef) {
	assign.push_back(new Assignment(g,e,ovp,ef));
	return;
}


//OK
void Gene::calcPrexpressPerExonUsingTag(int readl, int activetag) {
	// sum the uniscore and score for tagged exons
	double scoresum = 0, uniqsumsum = 0, taggedexons = 0;
	for (int i=0; i < this->exonNumber(); i++) {
		if (this->getExon(i)->tags() >= activetag) {
			scoresum   += this->getExon(i)->overlappingReads(readl);
			uniqsumsum += (double) this->getExon(i)->uniscore();
			taggedexons++;
		}
	}
	double genePrexpress = 0.0;
	if (uniqsumsum > 0) genePrexpress = scoresum/uniqsumsum;
	
	
	// update prexpress value per exon
	for (int i=0; i < this->exonNumber(); i++) {
		double score   = this->getExon(i)->overlappingReads(readl);
		double uniqsum = (double) this->getExon(i)->uniscore();

		if (this->getExon(i)->tags() >= activetag && taggedexons > 0) {
			this->getExon(i)->setPrexpress(genePrexpress);
		}
		else {
			double exonPrexpress = 0.0;
			if (uniqsum > 0) exonPrexpress = score/uniqsum;
			this->getExon(i)->setPrexpress(exonPrexpress);
		}
	}
	return;
}



//LEGACY (DO NOT USE!)
void Gene::calcPrexpressPerExon(int readl) {
	// update prexpress value per exon
	for (int i=0; i < this->exonNumber(); i++) {
		double score   = this->getExon(i)->overlappingReads(readl);
		double uniqsum = (double) this->getExon(i)->uniscore();
		if (uniqsum == 0) this->getExon(i)->setPrexpress(0.0); // avoids div by zero and is analytically correct (conservative)
		else              this->getExon(i)->setPrexpress(score/uniqsum);
	}
	return;
}

//OK
double Exon::overlappingReads(int readl) { // RAPARM
	double score = 0.0;
	for (int i = 0 ; i < readl ; ++i) {
		if (this->partmap[i] == 0) continue;
		score += (double) this->partmap[i] * (readl - i) / readl;
	}
	return score;
}

//OK (checks if annotation overlaps ->overlapping redundant annotations are ignored in assignedGeneSum counting to avoid downweighting of duplicates)
bool Gene::compatibleWith (Gene * g) {
	// similar prexpress
	double bigger = (this->getPrexpressSum() < g->getPrexpressSum()) ? g->getPrexpressSum() : this->getPrexpressSum();
	if ((fabs(this->getPrexpressSum() - g->getPrexpressSum()) / bigger) < 0.1) { //less than 10% difference
		int ovpc = 0;
		if (this->start >= g->start) {
			if ( this->end <= g->end) ovpc = (this->end - this->start + 1);
			else ovpc = (g->end - this->start + 1);
		}
		else {
			if (this->end <= g->end) ovpc = (this->end - g->start + 1);
			else ovpc = (g->end - this->end + 1);
		}
		double sharedLengthFraction = (double) ( ovpc / ( (this->length() > g->length()) ? this->length() : g->length()));
		// overlapping more than 80%
		if (sharedLengthFraction > 0.8) return true;
	}
	return false;
}

//OK
int Region::buildIslands() {
	istack reverse; // the reverse islands (0<= are in initial vector)
	islands.clear(); //rebuild if rerun
	
	for(imap::const_iterator it = allMaps.begin(); it != allMaps.end(); ++it) {
		if (it->second->strand < 0) {
			if ( (reverse.size() == 0) || !(reverse.back()->canBeAddedthenAdd(it->first, it->second)) ) { // adding function raises flag automatically
				reverse.push_back(new Island(this, it->first, it->second)); //postion of start, mapping
			}
		}
		else {
			if ( (islands.size() == 0) || !(islands.back()->canBeAddedthenAdd(it->first, it->second)) ) { // adding function raises flag automatically
				islands.push_back(new Island(this, it->first, it->second)); //postion of start, mapping
			}
		}
	}
	
	// merge
	for (unsigned int r = 0; r != reverse.size(); ++r) islands.push_back(reverse[r]);
	reverse.clear();
	
	//order islands by leftend (to accelerate connect islands afterwards)
	jmap imm;
	for (unsigned int s = 0; s != islands.size(); ++s) {
		imm.insert(pair<int, Island*>(islands[s]->leftend(), islands[s]));
	}
	islands.clear();
	for (jmap::const_iterator it = imm.begin(); it != imm.end(); ++it) {
		islands.push_back(it->second);
	}
	return islands.size();
}

//OK (Evaluates if a mapping [slice] can be added to island and adds if yes)
bool Island::canBeAddedthenAdd(int first, Mapping * m) {
	bool added = false;
	if (m->region() != region->region()) return added;
	if (!(sameStrand(m->strand, strand))) return added;
	if (first <= to + 1 + idist && from <= m->sliceEndWithStart(first)) {
		if (first < from)                     from      = first;
		if (m->start() < leftmost)            leftmost  = m->start();
		if (m->sliceEndWithStart(first) > to) to        = m->sliceEndWithStart(first);
		if (m->end() > rightmost)             rightmost = m->end();
		this->updateLinkage(m);
		totalReads += m->sliceFracWithStart(first);
		allReads++;
		posmap.insert(imap::value_type(first, m));
		this->raiseFlag(m->getFlag()); // raises flag if applicable (-1 intergenic, 0 intronic, 1 exonic)
		added = true;
	}
	return added;
}

//OK (Find connected islands eg multiexonic genes)
int Region::connectIslands() {
	// use graph to connect
	using namespace boost ;
	typedef adjacency_list<vecS, vecS, undirectedS> Graph;
    Graph G;
	for (unsigned int i=0; i < islands.size(); i++) {
		if (islands[i]->linkable()) {
			for (unsigned int j=i; j < islands.size(); j++) {
				if (islands[i]->rightend() < islands[j]->leftend()) break;
				if (islands[i]->leftend() <= islands[j]->rightend() && islands[i]->rightend() >= islands[j]->leftend()) add_edge(i, j, G);// add edge
			}
		}
		else add_edge(i, i, G);// add edge
	}
	
	// find connected components
	std::vector<int> component(num_vertices(G));
	csize = connected_components(G, &component[0]); 
	jmap compo;
	std::vector<int>::size_type i;
	for (i = 0; i != component.size(); ++i) compo.insert( pair<int, Island*>(component[i], islands[i]) ); //<componentnumber,islands>
	
	// foreach component add csize index
	for (int c = 0; c < csize; ++c) { // foreach component
		pair<jmap::iterator,jmap::iterator> ret = compo.equal_range(c);
		// find highest splice class and its size
		int highFlag = -1;
		for (jmap::const_iterator cc = ret.first; cc != ret.second; ++cc) {
			cc->second->connexID(c+1);
			if (cc->second->getFlag() > highFlag) highFlag = cc->second->getFlag();
		}
		for (jmap::const_iterator cc = ret.first; cc != ret.second; ++cc) cc->second->raiseFlag(highFlag); // connected islands get a flag raise to highest flag
	}
	return csize;
}

//OK Evaluate if Gene Coverage
int Region::verifyIslands(bool extend) {
	// get best gene
	int bgIslands = 0;
	if (extend) for(unsigned int i = 0; i < islands.size(); i++) {
		if (islands[i]->evaluateGeneConverage(nullgene)) ++bgIslands;
		else islands[i]->bestGene(nullgene);
	}
	else        for(unsigned int i = 0; i < islands.size(); i++) islands[i]->bestGene(nullgene);
	return bgIslands;
}

//OK Classifies Island Type
int Region::typeIslands(int tag, int readlength) {
	int expressedRetro = 0;
	
	for (unsigned int i=0; i < islands.size(); i++) {
		if (islands[i]->getFlag() > 0) { //EXONIC
			// if island flagged as exonic set Rflag to E-wie-exon
			islands[i]->setRflag('E');
		}
		else if (islands[i]->bestGeneID() != 0) {

			// minimal Size
			bool minSize = (islands[i]->readCount() > 10 && islands[i]->length() > readlength) ? true : false;
		 		
			// ambiguous read dominance
			islands[i]->dominanceTest = (double)islands[i]->bestGeneReads() / (double)islands[i]->ambiguousReads();
			bool dominance = (islands[i]->dominanceTest > 0.5 && islands[i]->bestGeneReads() > 1) ? true : false;
			
			if (dominance && minSize) {
				bool expressed = (islands[i]->uniqueReads() > 0) ? true : false;
				if (expressed) {
					islands[i]->runTests();
					// processing evidence -> spliced mappings in parent / lengthening / multiexonic mappings
					bool isRetro = (islands[i]->splicedGeneMaps > 0 && islands[i]->lengthTest >= 1.5 && islands[i]->multiexonTest > 1) ? true : false;
					if (isRetro) {
						islands[i]->setRflag('R');
						expressedRetro++;
					}
					else {
						islands[i]->setRflag('F');
					}
				}
				else {
					islands[i]->setRflag('P');
				}
			}
			else {
				islands[i]->setRflag('-');
			}
		}
		else {
			islands[i]->setRflag('-');
		}
	}
	return expressedRetro;
}

//OK
bool Island::evaluateGeneConverage(Gene * ng) {
	// add to genemap
	if (flag > 0) return false; // if flagged as exon -> no need to find bestGene
	// reset bGene
	bGene = NULL;
	for(imap::const_iterator it = posmap.begin(); it != posmap.end(); ++it){ //(position mapping)
		for (unsigned int mp = 0; mp < it->second->getRead()->mappings.size(); mp++) {
			if (it->second->getRead()->readMapping(mp)->assignments() > 0) {
				gm.clear(); // reset islands gene counter
				for(int a=0; a < it->second->getRead()->readMapping(mp)->assignments(); a++) gm[it->second->getRead()->readMapping(mp)->getAssignment(a)->getGene()]++;
				for(gmap::const_iterator it2 = gm.begin(); it2 != gm.end(); ++it2) genemap[it2->first]++;
			}
		}
	}
	// evaluate genemap
	bGeneReads = 0;
	if (genemap.size() > 0) {
		for(gmap::iterator it = genemap.begin(); it != genemap.end(); ++it) {
			if (it->second > bGeneReads) {
				bGene      = it->first;
				bGeneReads = it->second;
			}
		}
		return true;
	}
	bGene = ng;
	return false;
}

void Island::runTests() {
	exonvec matchedExons;
	bool assignedToGene;
	splicedGeneMaps = 0;
	int islStart = -1;
	int islEnd   = -1;
	int bgeStart = -1;
	int bgeEnd   = -1;
	int bgreads = 0;

	for(imap::const_iterator it = posmap.begin(); it != posmap.end(); ++it){ //(position mapping)
		for (unsigned int mp = 0; mp < it->second->getRead()->mappings.size(); mp++) {
			//multiexon tester && lengthFactor
			if (it->second->getRead()->readMapping(mp)->assignments() > 0) {
				for(int a=0; a < it->second->getRead()->readMapping(mp)->assignments(); a++) {
					if (it->second->getRead()->readMapping(mp)->getAssignment(a)->getGeneID() == bGene->gid()) {
						bgreads++;
						if (islStart < 0)                                             islStart = it->first;
						if (islEnd   < it->second->sliceEndWithStart(it->first))      islEnd   = it->second->sliceEndWithStart(it->first);
						if (bgeStart < 0)                                             bgeStart = it->second->getRead()->readMapping(mp)->start();
						if (bgeEnd   < it->second->getRead()->readMapping(mp)->end()) bgeEnd   = it->second->getRead()->readMapping(mp)->end();
						matchedExons.push_back(it->second->getRead()->readMapping(mp)->getAssignment(a)->getExon());
					}
				}
			}
			//bestGeneSplicedMaps
			if (it->second == it->second->getRead()->readMapping(mp)) continue; // not itself
			if (it->second->getRead()->readMapping(mp)->spliced()) {
				assignedToGene = false;
				for (int a = 0; a < it->second->getRead()->readMapping(mp)->assignments(); a++) {
					if (it->second->getRead()->readMapping(mp)->getAssignment(a)->getGeneID() == bGene->gid()) {
						assignedToGene = true;
						break;
					}
				}
				if (assignedToGene && !(it->second->spliced())) splicedGeneMaps++;
			}
		}
	}

	// LengthFactor Eval
	if (bgreads < 2) {
		lengthTest = 1;
	}
	else if (islStart > 0 && islEnd > 0 && bgeStart > 0 && bgeEnd > 0) {
		lengthTest = ( double(bgeEnd - bgeStart + 1) / double(islEnd - islStart + 1) );
	}
	else abort();
		
	// Multiexon Eval
	using namespace boost ;
	typedef adjacency_list<vecS, vecS, undirectedS> Graph;
    Graph G;
	for(unsigned int xx = 0; xx < matchedExons.size(); ++xx) {
		for(unsigned int yy = xx; yy < matchedExons.size(); ++yy) {
			if ( xx == yy || matchedExons[xx]->joint(matchedExons[yy])) {
				add_edge(xx, yy, G);// add edge
			}
		}
	}
	std::vector<int> component(num_vertices(G));
	multiexonTest = connected_components(G, &component[0]); 
	
	// set provisional uniscore
	exonset exonList;
	for (unsigned int w=0; w < matchedExons.size(); ++w) exonList.insert(matchedExons[w]);
	int exonLengthSum   = 0;
	int exonUniscoreSum = 0;
	for(exonset::const_iterator ex = exonList.begin(); ex != exonList.end(); ++ex) {
		exonLengthSum += (*ex)->length();
		exonUniscoreSum += (*ex)->length();
	}
	double usc = double(exonUniscoreSum) / double(exonLengthSum);
	double doubleuniscore = usc * this->length();
	uniscore = int(doubleuniscore);

	return;
}

//OK (create retro and prepare remap)
Exon * Region::addIslandAnnotationAndResetReads(Island * isl, int newGeneId, int newExonId, int analysisTag, int rl) {
	int val[8];
	val[0] = newGeneId;
	val[1] = newExonId;
	val[2] = rid;
	val[3] = isl->start();
	val[4] = isl->end();
	val[5] = (isl->strand != 0) ? isl->strand : isl->estimateStrand();
	val[6] = analysisTag;
	val[7] = isl->estimateUniscore();

	Exon * newExon = new Exon(val);
	this->addExonToGene(newExon);
	
	//protocol added annotation
	
	//project unique reads
	for(imap::const_iterator it = isl->posmap.begin(); it != isl->posmap.end(); ++it) {
		if (it->second->getRead()->isUnambiguous()) this->projection(it->second, rl);
		else it->second->getRead()->resetForRemap(); 
	}
	//calc prexpress
	newExon->getGene()->calcPrexpressPerExon(rl);
	//return new annotation
	return newExon; // should be added to protocol
}

//OK (estim from mappings on sister gene)
int Island::estimateStrand() {
	int sameStrand = 0;
	int oppositeStrand = 0;
	for(imap::const_iterator it = posmap.begin(); it != posmap.end(); ++it){
		for (unsigned int mp = 0; mp < it->second->getRead()->mappings.size(); mp++) {
			for(int a=0; a < it->second->getRead()->readMapping(mp)->assignments(); a++) {
				if (it->second->getRead()->readMapping(mp)->getAssignment(a)->getGene()->gid() == bGene->gid()) {
					if (strandcmp(it->second->getRead()->readMapping(mp)->strand, it->second->strand)) ++sameStrand;
					else ++oppositeStrand;
				}
			}
		}
	}
	
	return (sameStrand > oppositeStrand) ? bGene->strand() : ( (sameStrand == oppositeStrand) ? 0 : (bGene->strand() * -1) );
}

//OK (calc from sister exons)
int Island::estimateUniscore() {
	exonset es;
	double uniscoreSum = 0.0, lengthSum = 0.0;
	for(imap::const_iterator it = posmap.begin(); it != posmap.end(); ++it){
		for (unsigned int mp = 0; mp < it->second->getRead()->mappings.size(); mp++) {
			for(int a=0; a < it->second->getRead()->readMapping(mp)->assignments(); a++) {
				if (it->second->getRead()->readMapping(mp)->getAssignment(a)->getGene()->gid() == bGene->gid()) {
					es.insert(it->second->getRead()->readMapping(mp)->getAssignment(a)->getExon());
				}
			}
		}
	}
	
	for(exonset::const_iterator it = es.begin(); it != es.end(); ++it){
		uniscoreSum += (double) (*it)->uniscore();
		lengthSum   += (double) (*it)->length();
	}	
	return rounddouble(uniscoreSum * ((double) this->length() / lengthSum ));
}

//OK (distribute ambiguous reads)
void Read::fractionate() {
	// only works if there are assignments (should be checked by calling routine)
	// get gene assignments and ther prexpress
	double expressum = 0;
	for (unsigned int m = 0; m < mappings.size(); m++) {
		if (this->readMapping(m)->assignments() != 0) {
			expressum += this->readMapping(m)->assignedGenesSum(); 
		} // skip if mapping has no gene assignment
	}
	if (expressum == 0) { // no genes have xpression evidence
		//prioritize spliced gene assignments (keep only highest splice class)
		
		// select highest splice class
		int maxSpliceClass = 0;
		for (unsigned int m = 0; m < mappings.size(); m++) {
			if (this->readMapping(m)->assignments() == 0) continue; // skip if mapping has no gene assignment
			if (this->readMapping(m)->gaps() > maxSpliceClass) maxSpliceClass = this->readMapping(m)->gaps();
		}
		
		// distribute evenly among them
		int scc = 0; // spliceclasscount
		for (unsigned int m = 0; m < mappings.size(); m++) {
			if (this->readMapping(m)->assignments() == 0) continue; // skip if mapping has no gene assignment
			if (this->readMapping(m)->gaps() == maxSpliceClass) scc++;
		}
		//fractionate in highest splice class
		for (unsigned int m = 0; m < mappings.size(); m++) {
			if (this->readMapping(m)->assignments() == 0) continue; // skip if mapping has no gene assignment
			if (this->readMapping(m)->gaps() == maxSpliceClass) {
				for (int a=0; a < this->readMapping(m)->assignments(); a++) {
					this->readMapping(m)->getAssignment(a)->setExpressfrac(1.0/double(scc));
				}
			}
			else { // will not be used -> set expressFrac to 0
				for (int a=0; a < this->readMapping(m)->assignments(); a++) {
					this->readMapping(m)->getAssignment(a)->setExpressfrac(0.0);
				}
			}
		}	
	}
	else { // raparm distribution  //OK
		for (unsigned int m = 0; m < mappings.size(); m++) {
			if (this->readMapping(m)->assignments() == 0) continue; // skip if mapping has no gene assignment -> is downweighted to zero
			double thisMappingWeight = this->readMapping(m)->assignedGenesSum() / expressum;
			for (int a=0; a < this->readMapping(m)->assignments(); a++) {
				this->readMapping(m)->getAssignment(a)->setExpressfrac(thisMappingWeight);
			}
		}	
	}
	return;
}

//OK
double Mapping::assignedGenesSum() {
	if (!(assign.size() > 0)) return 0;
	double ag = 0;
	bool redu;
	if (assign.size() > 1) { // more than 1 assignment
		for (unsigned int o = 0; o < (assign.size() - 1); o++) {
			redu = false;
			for (unsigned int b = o; b < assign.size(); ++b) {
				if (o == b) continue;
				if (assign[o]->getGene()->compatibleWith(assign[b]->getGene())) redu = true;
			}
			if (!(redu)) ag += assign[o]->getWeight();
		}
	}
	// add last gene anyway
	ag += assign[assign.size() - 1]->getWeight();
	return ag;
}

//OK (writes also to unicov via radd)
bool Read::calculateCoverage(Region ** regions, int rs) {
	double f  = 0;
	for(unsigned int p=0; p < mappings.size(); p++) f += mappings[p]->highestExpressFrac();
	
	if (this->isGenic() && (f > 0)) {
		// distribute among genes
		double rd = closed_interval_rand(0,f);
		f = 0;
		for(unsigned int p=0; p < mappings.size(); p++) {
			f += mappings[p]->highestExpressFrac();
			if (f >= rd) {
				return mappings[p]->radd(regions, rs);
			}
		}
	}
	else {
		// distribute among all possible mappings (even if genes are at 0 expressfrac)
		double rd = closed_interval_rand(0,1);
		f = 0;
		for(unsigned int p=0; p < mappings.size(); p++) {
			f += 1.0 / double(mappings.size());
			if (f >= rd) {
				return mappings[p]->radd(regions, rs);
			}
		}
	}
	return false;
}

//OK
double Mapping::highestExpressFrac() {
	if (!(assign.size() > 0)) return 0;
	double ag = 0;
	for (unsigned int o = 0; o < assign.size(); ++o) {
		ag = (assign[o]->getExpressFrac() > ag) ? assign[o]->getExpressFrac() : ag;
	}
	return ag;
}

//OK
bool Mapping::radd(Region ** regions, int rs) {
	for (int r = 0; r < rs; r++) {
		if (regionid == regions[r]->region()) {
			if (strand < 0) {
				for (intmap::const_iterator sl = slices->begin(); sl != slices->end(); ++sl) {
					regions[r]->revcovadd(sl->first, sl->second - sl->first + 1);
					if (this->getRead()->mappings.size() == 1) regions[r]->revunicovadd(sl->first, sl->second - sl->first + 1);
				}
				return true;
			}
			else {
				for (intmap::const_iterator sl = slices->begin(); sl != slices->end(); ++sl) {
					regions[r]->covadd(sl->first, sl->second - sl->first + 1);
					if (this->getRead()->mappings.size() == 1) regions[r]->unicovadd(sl->first, sl->second - sl->first + 1);
				}
				return true;
			}
		}
	}
	return false;
}






//OK
double Gene::rpkm(int total, int tg) {
	double ex = 0.0;
	if (this->tagged_length(tg) == 0) {
		ex = -1;
	}
	else if (mapscore > 0) {
		double el = (double) this->tagged_length(tg) / 1000.0;
		double rd = (double) total / 1000000.0;
		ex = mapscore / (el * rd);
	}
	return ex;
}



