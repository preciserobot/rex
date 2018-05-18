/*
 *  universe.hpp
 *  ames-ms/raparm-ms Object classes (with base methods)
 *
 *  Created by David Brawand on 30.03.10.
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
#include <algorithm>
#include <utility>
#include <map>
#include <math.h>
#include <assert.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/incremental_components.hpp> 
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>



#define MAXREGIONEXONS 1000000;
#define MAXGENEEXONS 1200;
#define MAXREADLENGTH 150;

using namespace std ;


/* ------------------------------------------------------------- */
/* Class declarations                                            */
/* ------------------------------------------------------------- */
class Read;
class Exon;
class Gene;
class Region;
class Mapping;
class Splice;

/* ------------------------------------------------------------- */
/* Type declarations                                             */
/* ------------------------------------------------------------- */
typedef map<int, int> intmap;
typedef map<int, Splice*> splicemap;
typedef vector<Mapping*> mapvec;
typedef set<int> intset;

class Exon {
public:
	Exon(int values[]) {
		geneid    = values[0];
		exonid    = values[1];
		regionid  = values[2];
		first     = values[3];
		last      = values[4];
		strand    = values[5];
		tag       = values[6];
		uniq      = 0.0;
		elength   = abs(last - first + 1);
		partmap   = new int [maxreadl];
		for (int i = 0; i < maxreadl; ++i) partmap[i] = 0;
		leftanchor  = 0;
		rightanchor = 0;
	}
	~Exon() {};

	bool assignOverlap(Mapping* mp, int readl); //AMES
	int uniqEstimate(int readl); //AMES
	int uniqEstimateBrute(int readl); //old method

	
	int gid() {
		return geneid;
	}
	int eid() {
		return exonid;
	}
	int length() {
		return elength;
	}
	int region() {
		return regionid;
	}
	int str() {
		return strand;
	}
	int start() {
		return first;
	}
	int end() {
		return last;
	}
	int tags() {
		return tag;
	}
	void addAnchor(int ancr) {
		anchors.insert(ancr);
	}
	//int uniscore() { return static_cast<int>(uniq + 0.5); }
	//float uniscore() { return uniq; }
	void setUniscore(float u) {
		uniq = u;
		return;
	}
	int mappings(int m) {
		return partmap[m];
	}
	void show() {
		std::cerr << "E " << regionid << "\t" << strand << "\t" << geneid << "\t" << exonid << std::endl;
		std::cerr << "\t" << first << "-" << last << std::endl;
		return;
	}
	friend class Read;
	friend class Region;
	unsigned int rightanchor;
	unsigned int leftanchor;
protected:
	int regionid;
	int geneid;      //the gene identifier (integer)
	int exonid;
	int first;
	int last;
	int strand;
	int tag;
	float uniq; // will be undefined -> constructor should give 0 value
	int elength;
	int *partmap;
	intset anchors;
	const static int maxreadl = MAXREADLENGTH;
};

class Gene {
public:
	int exonCount;
	int start;
	int end;
	Gene(int gid) {
		geneid    = gid;
		exons = new Exon *[maxexon];
		start = 9;
		end = 0;
		exonCount = 0;
		generegion = 0;
	}
	~Gene() {};
	int gid() {
		return geneid;
	}
	int length() {
		return (end - start + 1);
	}
	int region() {
		return generegion;
	}
	void addExon (Exon * ex) {
		if (start > end) {
			start = ex->start();
			end   = ex->end();
		}
		if (ex->start() < start) start = ex->start();
		if (ex->end()   > end)   end   = ex->end();
		exons[exonCount] = ex;
		exonCount++;
		if (maxexon <= exonCount) {
			std::cerr << "ERROR: Exceeded maximum exon number ("<< maxexon <<") for gene " << geneid << std::endl;
			exit(1);
		}		
		// verify that all exons have same region when added to gene
		for (int e=0; e < exonCount; e++) {
			if (generegion == 0) {
				generegion = ex->region();
			}
			if (generegion != ex->region()) {
				std::cerr << "ERROR: Inconsistent exon " << exons[e]->eid() << std::endl;
				exit(1);
			}
		}
		return;
	}
	void show() {
		std::cerr << "G " << geneid << "\t" << generegion << "\t" << exonCount << std::endl;
		std::cerr << "\t" << start << "-" << end << std::endl;
		return;
	}
	Exon * getExon(int number) {
		return exons[number];
	}
	int exonNumber() {
		return exonCount;
	}
	friend class Exon;
protected:
	int geneid;      //the gene identifier (integer)
	const static int maxexon = MAXGENEEXONS;
	Exon **exons;
	int generegion;
};

class Region {
public:
	Region(){};
	Region(int val){
		rid       = val;
	};
	~Region(){};
	
	bool projection(Mapping * mp, int readl); // AMES/RAPARM
	bool spliceProjection(int a, int z); // AMES (to weigh splcied mappings)
	int dicoindex(Mapping * mp); //COMMON
	int dicoindex(int i) {
		int begin = 0;
		int end = gsize;
		int mid = (gsize+begin)/2;
		
		while (begin != end) {
			if (garray[mid]->end < i) {
				begin = mid+1;
			}
			else if (garray[mid]->end > i) {
				end = mid; 
			}
			else {
				while (mid >= 0 && garray[mid]->end == i) {
					mid--;
				}
				mid++;
				end = mid;
				begin = mid;
			}
			mid = (end+begin)/2;
		}
		return mid;
	}
	
	
	int region() {
		return rid;
	}
	int genenumber() {
		return regiongenes.size();
	}
	int exonnumber() {
		int ec = 0;
		for(gmap::const_iterator it = regiongenes.begin(); it != regiongenes.end(); ++it){
			ec += it->second->exonCount;
		}
		return ec;
	}
	void show() {
		std::cerr << "R " << rid << std::endl;
		return;
	}
	
	void addExonToGene(Exon * ex) {
		gmap::iterator iter = regiongenes.find(ex->geneid);
		if( iter == regiongenes.end() ) regiongenes.insert(pair<int, Gene*>(ex->geneid, new Gene(ex->geneid)));
		
		if (regiongenes[ex->geneid]->gid() != ex->geneid || ex->region() != rid) {
			if (regiongenes[ex->geneid]->gid() != ex->geneid) std::cerr << "HUH?";
			if (ex->region() != rid) std::cerr << "HAH?";
			
			this->show();
			regiongenes[ex->geneid]->show();
			ex->show();

			abort();
		}
		
		regiongenes[ex->geneid]->addExon(ex);
	}
	void orderGenesByStart() {
		garray = new Gene * [regiongenes.size()];
		gsize = 0;
		maxGeneLength = 0;
		if (regiongenes.size() == 0) {
			std::cerr << "WARNING: No gene in region " << rid << std::endl;
			return;
		}
		
		mgmap orderedregiongenes;
		for(gmap::const_iterator it = regiongenes.begin(); it != regiongenes.end(); ++it){
			orderedregiongenes.insert(pair<int,Gene*>(it->second->end, it->second));
		}
		
		for(mgmap::const_iterator it = orderedregiongenes.begin(); it != orderedregiongenes.end(); ++it){
			garray[gsize++] = it->second;
			if (it->second->length() > maxGeneLength) maxGeneLength = it->second->length();
		}
		return;
	}
	
	friend class Exon;
	intmap introns;
protected:
	int	rid;
	typedef map<int, Gene*> gmap;
	typedef multimap<int, Gene*> mgmap;
	gmap regiongenes;
	Gene **garray;
	int gsize;
	int maxGeneLength;
};

class Splice {
public:
	Splice(){};
	Splice(int msid, int sr, int st){
		spliceID  = msid;
		seqRegion = sr;
		strand    = st;
	};
	~Splice(){};
	
	void addSlice(int from, int to) {
		slices.insert(pair <int,int> (from, to));
		return;
	}
	int region() {
		return seqRegion;
	}
	intmap * resolvePosition(int p, int l) {
		intmap * mslc = new intmap;
		int offset = 0;
		for (intmap::const_iterator s = slices.begin(); s != slices.end(); ++s) {
			if (offset < p && p <= offset + s->second - s->first + 1) {
				do {
					int start = s->first + (p - offset - 1);
					int end   = (l > s->second - start + 1) ? (s->second) : (start + l - 1);
					p        += end - start + 1;
					l        -= end - start + 1;
					offset   += s->second - s->first + 1;
					mslc->insert(pair <int,int> (start,end));
					++s;
				} while (l > 0);
				if (l == 0) return mslc;
				else abort();
			}
			offset += s->second - s->first + 1;
		}
		std::cerr << "ERROR: Mapping and Multisplice do not match!" << std::endl;
		abort();
	}
	
	bool centralintron(int * intr) {
		intr[0] = seqRegion;
		int splicesum = 0;
		for (intmap::const_iterator sl = slices.begin(); sl != slices.end(); ++sl) {
			splicesum += (sl->second - sl->first + 1);
		}
		int overh = splicesum / 2;
		splicesum = 0;
		for (intmap::const_iterator sl = slices.begin(); sl != slices.end(); ++sl) {
			if (splicesum == overh) {
				intr[2] = sl->first - 1;
				return true;
			}
			splicesum += (sl->second - sl->first + 1);
			intr[1] = sl->second + 1;
		}
		return false;
	}
	
protected:
	int spliceID;
	int seqRegion;
	int strand;
	intmap slices;
};

class Mapping {
public:
	Mapping(){};
	Mapping(int regid, int str, intmap * slc){
		regionid = regid;
		strand   = str;
		slices   = slc;
		accepted = false;
		//std::cerr << " *";
	};
	Mapping(int regid, int str, int pos, int readl){
		regionid = regid;
		strand   = str;
		slices   = new intmap;
		slices->insert(pair<int, int> (pos,pos+readl-1));
		accepted = false;
		//std::cerr << " x";
	};
	~Mapping(){};
	bool operator==(const Mapping& other) const {
		if (regionid != other.regionid) return false;
		for(intmap::const_iterator aa = slices->begin(); aa != slices->end(); ++aa) {
			if (other.slices->find(aa->first) == other.slices->end() || other.slices->find(aa->first)->second != aa->second) return false;
		}
		return true;
	}
	void show() {
		std::cerr << "M " << regionid << "\t" << strand << std::endl;
		for (intmap::const_iterator xx = slices->begin(); xx != slices->end(); ++xx) {
			std::cerr << "\t" << xx->first << "-" << xx->second << std::endl;
		}
		return;
	}
	bool overlaps(Gene * g) {
		// works also for spliced mappings
		if (((g->start <= this->start()) && (this->start() <= g->end)) || ((g->start <= this->end()) && (this->end() <= g->end))) {
			return true;
		}
		if ((this->start() < g->start) && (g->end < this->end())) {
			return true;
		}
		return false;
	}
	bool regiomap(Region* r) {
		if (regionid == r->region()) return true;
		return false;
	}
	int maplength() {
		int len = 0;
		for(intmap::const_iterator it = slices->begin(); it != slices->end(); ++it) {
			len += it->second;
		}
		return len;
	}
	int region() {
		return regionid;
	}
	int start() {
		return slices->begin()->first;
	}
	int end() {
		intmap::iterator it = slices->end();
		--it;
		return it->second;
	}
	bool spliced() {
		if (slices->size() > 1) {
			return true;
		}
		return false;
	}
	bool overlaps(Mapping* m) {
		for (intmap::const_iterator xx = this->slices->begin(); xx != this->slices->end(); ++xx) {
			for (intmap::const_iterator yy = m->slices->begin(); yy != m->slices->end(); ++yy) {
				if ((yy->first <= xx->first && xx->first <= yy->second) ||
					(yy->first <= xx->second && xx->second <= yy->second) ||
					(xx->first <= yy->first && yy->second <= xx->second)) {
					return true;
				}
			}
		}
		return false;
	}
	void accept() {
		accepted = true;
		return;
	}
	int gaps() {
		return slices->size() - 1;
	}
	friend class Read;
	friend class Exon;
	bool accepted;	
protected:
	int	regionid;
	int strand;
	intmap * slices;
};

class Read {
public:
	Read(){};
	Read(int identifier) {
		readid = identifier;
		uniq = 0.0;
	};
	~Read(){
		// destruct all mappings
	}
	
	bool unify(); // AMES-MS
	
	bool isUnique() {
		if (mappings.size() > 1) return false;
		else if (mappings.size() == 1) {
			mappings[0]->accept();
			return true;
		}
		else abort();
	}
	int countMappings() {
		return mappings.size();
	}
	void addFullMapping(int r, int s, int p, int l) {
		Mapping * map_n = new Mapping(r,s,p,l);
		if (!(this->contains(map_n))) mappings.push_back(map_n);
		else abort();
		return;
	}	
	void addSpliceMapping(int s, int p, int l, Splice * sp) { // region, strand, position, readlength, splice
		int r = sp->region();
		intmap * sm = sp->resolvePosition(p,l);
		Mapping * map_n = new Mapping(r,s,sm);
		if (!(this->contains(map_n))) mappings.push_back(map_n);
		else abort();
		return;
	}
	int rid() {
		return readid;
	}
	Mapping * trustedMapping() {
		Mapping * ret = NULL;
		for (unsigned int m = 0; m != mappings.size(); ++m) {
			if (mappings[m]->accepted) {
				if (ret != NULL) {
					cerr << "FATAL: Read has more than one accepted mapping" << endl;
					abort();
				}
				ret = mappings[m]; 
			}
		}
		if (ret == NULL) cerr << "WARNING: Called trustedMapping on unresolved Read" << endl;
		return ret;
	}
	bool contains(Mapping * mp) {
		for(unsigned int v = 0; v != mappings.size(); ++v) {
			if (mp == mappings[v]) return true;
		}
		return false;
	}
	friend class Mapping;
protected:
	int readid;
	float uniq;
	mapvec mappings;
};


  //////////////////////////////
 /////// COMMON METHODS ///////
//////////////////////////////



int Region::dicoindex(Mapping * mp) {
	int begin = 0;
	int end = gsize;
	int mid = (gsize+begin)/2;
	
	while (begin != end) {
		if (garray[mid]->end < mp->start()) {
			begin = mid+1;
		}
		else if (garray[mid]->end > mp->start()) {
			end = mid; 
		}
		else {
			while (mid >= 0 && garray[mid]->end == mp->start()) {
				mid--;
			}
			mid++;
			end = mid;
			begin = mid;
		}
		mid = (end+begin)/2;
	}
	return mid;
}


