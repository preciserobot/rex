/*
 *  universe.hpp
 *  ames-ms/raparm-ms/rm-ms Object classes (with base methods)
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
#include <sstream>
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

using namespace std;

/* ------------------------------------------------------------- */
/* DEFINITIONS                                                   */
/* ------------------------------------------------------------- */
#define MAXTAG 3;
#define MAXREADLENGTH 150;
#define MAXREGIONEXONS 1000000;
#define MAXGENEEXONS 1000;
#define	ISLANDIST 0;            // max distance between reads to be merged into an island
/* ------------------------------------------------------------- */
/* Class declarations                                            */
/* ------------------------------------------------------------- */
class Exon;
class Gene;
class Assignment;
class Region;
class Read;
class Mapping;
class Splice;
class Island;
/* ------------------------------------------------------------- */
/* Type declarations                                             */
/* ------------------------------------------------------------- */
typedef map<int, int> intmap;
typedef map<int, double> intdoublemap;
typedef map<int, Splice*> splicemap;
typedef vector<Mapping*> mapvec;
typedef vector<Assignment*> assignvec;
typedef vector<int> intvec;
typedef map<int, Gene*> gmap;
typedef multimap<int, Gene*> mgmap;
typedef multimap<int, Mapping*> imap;
typedef multimap<int, Island*> jmap;
typedef multimap<int, int> intmultimap;
typedef map<int, int> covmap;
typedef set<Mapping*> mapset;
typedef set<Exon*> exonset;
typedef vector<Exon*> exonvec;
typedef vector<Island*> istack;
/* ------------------------------------------------------------- */
/* inline functions                                              */
/* ------------------------------------------------------------- */
inline bool sameStrand(int a, int b) {
	if (abs(a-b) == 2) return false;
	else return true;
}
/* ------------------------------------------------------------- */
/* Class declarations                                            */
/* ------------------------------------------------------------- */
//OK
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
		uniq      = values[7];
		elength   = abs(last - first + 1);
		partmap   = new int [maxreadl];
		prexpress = 0;
		for (int i = 0; i < maxreadl; ++i) partmap[i] = 0;
	}
	~Exon() {};
	
	void updateRead(Read* rd); // RAPARM
	bool maps(Mapping* mp); // RAPARM
	void addMap(Mapping* mp, int readl, int selectedtag); // RAPARM
	void addWeightedMap(Mapping* mp, int readl); // RAPARM
	double overlappingReads(int readl); // RAPARM
	bool assignOverlap(Mapping* m, int readl); //RAPARM (replaces 2 older functions [maps&&addMap])

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
	void end(int e) {
		last = e;
		return;
	}
	int tags() {
		return tag;
	}
	int uniscore() {
		return uniq;
	}
	int mappings(int m) {
		return partmap[m];
	}
	void setGene (Gene * g) {
		genepointer = g;
	}
	Gene * getGene () {
		return genepointer;
	}
	void setPrexpress(double pex) {
		prexpress = pex;
		return;
	}
	double getPrexpress() {
		return prexpress;
	}
	bool joint (Exon * ex) {
		if (ex->start() < this->end() && ex->end() > this->start() ) {
			return true;
		}
		else if (this->start() - 1 == ex->end() || ex->start() - 1 == this->end()) {
			return true;
		}
		return false;
	}
	friend class Read;
	friend class Region;
protected:
	int regionid;
	int geneid;      //the gene identifier (integer)
	int exonid;
	int first;
	int last;
	int strand;
	int tag;
	int uniq;
	int elength;
	int *partmap;
	double prexpress; // uniqmap expression
	Gene * genepointer;
	const static int maxreadl = MAXREADLENGTH;
};

//OK
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
		mapscore   = 0.0;
	}
	~Gene() {};
	
	void addExon (Exon * ex); // RAPARM
	void calcPrexpressPerExon(int readl); // RAPARM // WiP OK?
	void calcPrexpressPerExonUsingTag(int readl, int activetag); // REPLACEMENT FOR PREXPRESS??? (if exons has no information use tagged exons as reference (average out on lowly expressed genes)

	int exonAddOverlaps(Exon * x); // RAPARM
	double rpkm(int total, int tg); // RAPARM
	
	bool compatibleWith (Gene * g); //CHECK

	int strand() {
		int str = 0;
		for (int e = 0; e != exonCount; ++e) {
			if ( strandcmp(exons[e]->str(), str) ) {
				if (exons[e]->str() != 0) str = exons[e]->str();
			}
			else {
				cerr << "ABORT TRAP: Gene has Exons on opposite strands!" << endl;
				abort();
			}
		}
		return str;
	}
	int gid() {
		return geneid;
	}
	int length() {
		return (end - start + 1);
	}
	int tagged_length (int t) {
		if (geneid == 0) return 0;
		int tl = 0;
		for (int e=0; e < exonCount; e++) if (exons[e]->tags() >= t) tl += exons[e]->length();
		// IF tl is zero gene will be ignored!
		//if (tl == 0) std::cerr << std::endl << "WARNING: zero tag length for gene " << geneid << " with " << exonCount << " exons" << std::endl;
		return tl;
	}
	int region() {
		return generegion;
	}
	Exon * getExon(int number) {
		return exons[number];
	}
	int exonNumber() {
		return exonCount;
	}
	int mappedReads() {
		return (int) (mapscore + 0.5);
	}
	void addMapscore(double o, double e) {
		mapscore += (double)o*(double)e;
		return;
	}
	double getPrexpressSum() {
		double sum = 0;
		for (int e=0; e < exonCount; e++) {
			sum += exons[e]->getPrexpress();
		}
		return sum;
	}
	string location () {
		stringstream loc;
		loc << start << "\t" << end << "\t" << this->strand();
		return loc.str();
	}
	friend class Exon;
protected:
	int geneid;      //the gene identifier (integer)
	const static int maxexon = MAXGENEEXONS;
	Exon **exons;
	int generegion;
	double mapscore; // the mapping score (sum of all overlapFrac*expressFrac assignments for this gene) >>(divided by tagged length and mapped reads ->RPKM)
};

//OK
class Assignment {
public:
	Assignment(Gene * g, Exon * exo, double o, double e) {
		gene = g;
		exon = exo;
		overlapFrac = o;
		expressFrac = e; //is -1 for ambig mapping if not fractionated yet, 1 for unambiguous
	};
	~Assignment(){};
	int getGeneID () {
		return gene->gid();
	}
	Gene * getGene () {
		return gene;
	}
	int getExonID () {
		return exon->eid();
	}
	Exon * getExon () {
		return exon;
	}
	double getExpressFrac () {
		return expressFrac;
	}
	double getWeight() {
		return exon->getPrexpress() * overlapFrac;
	}
	double getOverlapFrac (int tt) {
		return (exon->tags() >= tt) ? overlapFrac : 0.0;
	}
	void setExpressfrac(double f) {
		expressFrac = f;
		return;
	}
	bool commitExpression(int t) {
		double ovf = this->getOverlapFrac(t);
		if (expressFrac < 0.0) {
			std::cerr << "WARNING: Read may not have been fractinated yet (" << expressFrac << ") => setting to 0 !";
			gene->addMapscore(ovf, 0);
		}
		else gene->addMapscore(ovf, expressFrac);
		if (ovf > 0.0) return true;
		else return false;
	}
	bool singleExonAssign() {
		if (gene->exonNumber() == 1) return true;
		return false;
	}
	bool multiExonAssign() {
		if (gene->exonNumber() != 1) return true;
		return false;
	}
//?	friend class Mapping;
protected:
	Gene * gene;
	Exon * exon;
	double expressFrac;
	double overlapFrac;
};

//OK
class Region {
public:
	Region(){};
	Region(int val){
		rid   = val;
		gsize = 0;
		csize = -1;
		maxGeneLength = 0;
		nullgene = new Gene(0);
	};
	~Region(){};
	
	void projection(Mapping * mp, int readl); // OKOKOK
	
	int buildIslands(); // RAPARM
	int connectIslands();  // DEVELOPMENT
	int verifyIslands(bool extend); //sets the bestgene if requested
	int typeIslands(int tag, int readlength); // RAPARM
	Exon * addIslandAnnotationAndResetReads(Island * isl, int newExonId, int newGeneId, int analysisTag, int rl); //RETRO_ADD
	int dicoindex(Mapping * mp); //COMMON
	
	void covadd(int position, int length) {
		cov[position]++;
		cov[position+length]--;
		return;
	}
	void revcovadd(int position, int length) {
		revcov[position]++;
		revcov[position+length]--;
		return;
	}
	void unicovadd(int position, int length) {
		unicov[position]++;
		unicov[position+length]--;
		return;
	}
	void revunicovadd(int position, int length) {
		revunicov[position]++;
		revunicov[position+length]--;
		return;
	}
	
	int region() {
		return rid;
	}
	int genenumber() {
		return regiongenes.size();
	}
	Gene * getGene(int g) {
		return garray[g];
	}
	int exonnumber() {
		int ec = 0;
		for(gmap::const_iterator it = regiongenes.begin(); it != regiongenes.end(); ++it){
			ec += it->second->exonCount;
		}
		return ec;
	}
	void addExonToGene(Exon * ex) {
		gmap::iterator iter = regiongenes.find(ex->geneid);
		if( iter == regiongenes.end() ) regiongenes.insert(pair<int, Gene*>(ex->geneid, new Gene(ex->geneid)));
		regiongenes[ex->geneid]->addExon(ex);
		return;
	}
	int orderGenesByEnd(int tg) {
		int untagged = 0;
		if (regiongenes.size() == 0) {
			return untagged;
		}
		mgmap orderedregiongenes;
		for(gmap::const_iterator it = regiongenes.begin(); it != regiongenes.end(); ++it){
			if (it->second->tagged_length(tg) == 0) untagged++;
			orderedregiongenes.insert(pair<int,Gene*>(it->second->end, it->second));
		}
		garray = new Gene * [regiongenes.size()];
		gsize = 0;
		maxGeneLength = 0;
		for(mgmap::const_iterator it = orderedregiongenes.begin(); it != orderedregiongenes.end(); ++it){
			garray[gsize++] = it->second;
			if (it->second->length() > maxGeneLength) maxGeneLength = it->second->length();
		}
		return untagged;
	}
	
	int gCount() {
		return gsize;
	}

	int iCount () {
		return islands.size();
	}
	Island * getIsland(int i) {
		return islands[i];
	}
	void addIsland(Island * il) {
		islands.push_back(il);
		return;
	}
	friend class Exon;
	
	covmap cov;       // Coverage (pos,delta?)
	covmap revcov;    // Coverage (pos,delta?)
	covmap unicov;    // Coverage for unique reads
	covmap revunicov; // Coverage for unique reads
protected:
	Gene * nullgene; // the NULL gene
	int	rid;
	gmap regiongenes;
	imap allMaps;
	Gene **garray;
	int gsize;
	int maxGeneLength;
	istack islands; // island containing vector
	int csize; //connected islands
};

//OK
class Splice {
public:
	Splice(){};
	Splice(int msid, int sr, int st){
		spliceID  = msid;
		seqRegion = sr;
		strand    = st;
	};
	~Splice(){};
	
	void addSlice(int from, int to, int intronid) {
		if (from > to) {
			std::cerr << "Reversed slice";
			abort();
		}
		slices.insert(pair <int,int> (from, to));
		if (intronid != 0) introns[to] = intronid;
		return;
	}
	int region() {
		return seqRegion;
	}
	int str() {
		return strand;
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
		std::cerr << "ERROR: Mapping and Multisplice do not match! (" << spliceID << "|" << p << "|" << l << ")" << std::endl;
		for (intmap::const_iterator s = slices.begin(); s != slices.end(); ++s) cerr << "\t" << s->first << "\t" << s->second << endl;
		abort();
	}
	intvec * spannedIntrons(intmap * sli) {
		intvec * spi = new intvec;
		// get the last slice
		intmap::iterator lastslice = sli->end();
		--lastslice;
		for (intmap::const_iterator s = sli->begin(); s != sli->end(); ++s) {
			if (introns.find(s->second) != introns.end() && lastslice->second != s->second) {
				spi->push_back(introns.find(s->second)->second);
			}
		}
		return spi;
	}
	
protected:
	int spliceID;
	int seqRegion;
	int strand;
	intmap slices;
	intmap introns;
};

//OK
class Mapping {
public:
	Mapping(){};
	Mapping(Read * rd, int regid, int str, intmap * slc, int mm){
		regionid = regid;
		strand   = str;
		slices   = slc;
		read     = rd;
		mism     = mm;
		accepted = false;
		splice = NULL;
		flag = -1; // -1 sets not evaluated (or intergenic)
	};
	Mapping(Read * rd, int regid, int str, int pos, int readl, int mm){
		regionid = regid;
		strand   = str;
		slices   = new intmap;
		slices->insert(pair<int, int> (pos,pos+readl-1));
		read     = rd;
		mism     = mm;
		accepted = false;
		splice = NULL;
		flag = -1; // -1 sets not evaluated (or intergenic)
	};
	~Mapping(){
		delete slices;
		assign.clear();
	};
	
	Gene * assignedGene();
	double highestExpressFrac();
	double assignedGenesSum();
	void addAssignment(Gene * g, Exon * e, double ovp, double ef);
	bool radd(Region ** regions, int rs);

	bool operator==(const Mapping * other) const {
		if (regionid != other->regionid) return false;
		if (!(sameStrand(strand, other->strand))) return false;
		for(intmap::const_iterator aa = slices->begin(); aa != slices->end(); ++aa) {
			if (other->slices->find(aa->first) == other->slices->end() || other->slices->find(aa->first)->second != aa->second) return false;
		}
		return true;
	}
	void addSplice(Splice * spl) {
		splice = spl;
		return;
	}
	bool isUnspliced() {
		if (splice == NULL) return true;
		else return false;
	}
	Splice * getSplice() {
		return splice;
	}
	int str() {
		return strand;
	}
	void show() {
		std::cerr << regionid << "\t" << strand << std::endl;
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
	bool regiomap(Region * r) {
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
	bool notAccepted() {
		return !(accepted);
	}
	int gaps() {
		return slices->size() - 1;
	}
	void setFlag(int f) {
		// flags will never be downgraded
		if (flag < f)flag = f;
		return;
	}
	int getFlag() {
		return flag;
	}
	Read * getRead () {
		return read;
	}
	int assignments() {
		return assign.size();
	}
	Assignment * getAssignment(int a) {
		return assign[a];
	}
	int singleExon() {
		if (!(assign.size() > 0)) return 0;
		int se = 0;
		for (unsigned int o = 0; o < assign.size(); ++o) {
			if (assign[o]->singleExonAssign()) se++;
		}
		return se;
	}
	bool multiExonAssign() {
		if (!(assign.size() > 0)) return false;
		for (unsigned int o = 0; o < assign.size(); ++o) if (assign[o]->multiExonAssign()) return true;
		return false;
	}
	int sliceEndWithStart(int startsWith) {
		for (intmap::const_iterator it = slices->begin(); it != slices->end(); ++it) {
			if (it->first == startsWith) return it->second;
		}
		int rr = startsWith + this->maplength() - 1;
		std::cerr << "WARNING: Asked for a nonexistent mapping slice, returned " << rr << std::endl;
		return rr;
	}
	double sliceFracWithStart(int startsWith) {
		double sliceLength = 0.0;
		double totalLength = 0.0;
		for (intmap::const_iterator it = slices->begin(); it != slices->end(); ++it) {
			if (it->first == startsWith) sliceLength = it->second - it->first + 1;
			totalLength += it->second - it->first + 1;
		}
		return sliceLength / totalLength;
	}
	string cigarline() {
		stringstream ss;
		ss << "(";
		for (intmap::const_iterator it = slices->begin(); it != slices->end(); ++it) {
			ss << "." << it->first << "-" << it->second << ".";
		}
		ss << ")";
		
		return ss.str();
	}
	void removeAssignments() {
		for (unsigned int a = 0; a != assign.size(); ++a) delete assign[a];
		assign.clear();
	}
	intmap * getSlices() {
		return slices;
	}
	friend class Exon;   //OK
	friend class Region; //OK
	friend class Island; //OK
	friend class Read;   //OK
	bool accepted;	
protected:
	intmap * slices; // move to protected again
	int	regionid;
	int strand;
	int flag; // -1:intergenic 0:genic 1:exonic
	int mism; // number of mismatches
	Read * read;
	Splice * splice;
	assignvec assign;
};

//OK
class Read {
public:
	Read(){};
	Read(int identifier) {
		readid    = identifier;
		remapFlag = false;
	}
	~Read(){
		// destruct all mappings
		//delete mappings;
	}

	void fractionate(); // OK
	void equalize();    // OBSOLETE AS COVERAGE ROUTINE WILL TAKE CARE
	int  unify(bool leastMismatch); // WiP
	bool calculateCoverage(Region ** regions, int rs);
	
	//FUNCTIONS TO CHECK FOR NECESSITY
	//FUNCTIONS TO CHECK FOR NECESSITY
	//FUNCTIONS TO CHECK FOR NECESSITY
	int totalAssignments() {
		int as = 0;
		for(unsigned int m=0; m < mappings.size(); m++) as += mappings[m]->assignments();
		return as;
	}
	//FUNCTIONS TO CHECK FOR NECESSITY
	//FUNCTIONS TO CHECK FOR NECESSITY
	//FUNCTIONS TO CHECK FOR NECESSITY
	
	bool isUnambiguous() {
		if (mappings.size() == 1) {
			return true;
		}
		return false;
	}
	int countMappings() {
		return mappings.size();
	}
	void addFullMapping(int r, int s, int p, int l, int mm) { // region, strand, start, readlength
		Mapping * map_n = new Mapping(this,r,s,p,l,mm);
		if (!(this->contains(map_n))) mappings.push_back(map_n);
		else abort();
		return;
	}	
	void addSpliceMapping(Splice * sp, int s, int p, int l, int mm) { // region, strand, position, readlength, splice
		//cerr << " A" << sp ;
		int r = sp->region();
		//cerr << " B";
		intmap * sm = sp->resolvePosition(p,l);
		//cerr << " C";
		Mapping * map_n = new Mapping(this,r,s,sm,mm);
		//cerr << " D";
		map_n->addSplice(sp);
		//cerr << " E";
		if (!(this->contains(map_n))) mappings.push_back(map_n);
		else abort();
		//cerr << " F";
		return;
	}
	int rid() {
		return readid;
	}
	Mapping * readMapping(int mid) {
		return mappings[mid];
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

	bool isGenic() {
		int exom = 0;
		for (unsigned int m = 0; m < mappings.size(); m++) if (mappings[m]->assignments() > 0) exom++;
		if (exom > 0) return true;
		return false;
	}
	void resetForRemap() { // marks reads for remapping (after retro detection)
		remapFlag = true;
		for(unsigned int i = 0; i != mappings.size(); ++i) mappings[i]->removeAssignments();
		return;
	}
	bool flaggedForRemap() {
		return remapFlag;
	}
	friend class Mapping;
	friend class Island;
protected:
	bool remapFlag;
	int readid;
	double uniq;
	mapvec mappings;
};

//OK
class Island {
public:
	Island(Region * r, int firstpos, Mapping * m){
		region     = r;
		totalReads = m->sliceFracWithStart(firstpos);
		allReads   = 1;
		bGene      = NULL;
		bGeneReads = 0;
		from       = firstpos;
		to         = m->sliceEndWithStart(firstpos);
		leftmost   = m->start();
		rightmost  = m->end();
		connex     = -1;
		strand     = m->strand;
		flag       = m->getFlag();
		islandID   = r->iCount();
		rFlag      = '-';
		posmap.insert(imap::value_type(firstpos, m));
		// defined by other methods
		uniscore      = -1;
		dominanceTest = -1;
		lengthTest    = -1;
		multiexonTest = -1;
		//linkability
		linkage = false;
		updateLinkage(m);
	}
	~Island() {}
		
	bool canBeAddedthenAdd(int first, Mapping * m);
	bool evaluateGeneConverage(Gene * ng); //return true if bestgene found
	
	// merged 3 test methods into one
	void runTests();
	// annotation methods
	int estimateStrand();
	int estimateUniscore();
	
	string GFFline(int readlength) {
		stringstream ss;
		ss << setprecision(3);
		ss << region->region() << "\t"; //CHROMOSOME
		ss << "RNASEQ\t"; //SOURCE
		if      (rFlag == 'E') ss << "EXON\t";
		else if (rFlag == 'R') ss << "RETROGENE\t";
		else if (rFlag == 'F') ss << "FRAGMENT\t";
		else if (rFlag == 'P') ss << "PSEUDOGENE\t";
		else if (rFlag == '-') ss << "UNKNOWN\t";
		ss << from << "\t"; //START
		ss << to << "\t"; //END
		ss << this->coverage(readlength) << "\t"; //SCORE (COVERAGE)
		if (strand != 0) ss << strand << "\t"; //STRAND
		else             ss << "."    << "\t"; //STRAND
		ss << "." << "\t"; //FRAME
		ss << "Group " << connex; //GROUPINGID
		ss << " ; Island "        << islandID ;  //ISLAND_ID
		ss << " ; Type "          << flag ;      //TYPE (-1|0|+1)
		ss << " ; Length "        << to - from + 1 ; //LENGTH
		ss << " ; Reads "         << this->readCount() ; //READS
		ss << " ; Distinct "      << this->distinctReads() ; //DISTINCT_READS
		ss << " ; Unique "        << this->uniqueReads() ; //UNIQUE_READS
		ss << " ; Linkable "      << (this->linkable() ? "YES" : "NO"); //UNIQUE_READS
		if (this->bestGeneID() != 0) {
			if (flag != 1) {
				ss << " ; Gene "          << this->bestGeneID() ; //BESTGENE
				ss << " ; GeneReads "     << this->bestGeneReads() ; //BESTGENE_READS
				ss << " ; GeneSpliced "   << splicedGeneMaps ; //BESTGENE_SPLICED
				ss << " ; dominanceTest " << dominanceTest ; //FRAC OF AMBIGUOUS MAPPINGS IN BESTGENE
				ss << " ; lengthTest "    << lengthTest ; //LENGTH FRAC BESTGENE/ISLAND
				ss << " ; multiexonTest " << multiexonTest ; //MULTIEXONS
			}
		}
		return ss.str();
	}
	void updateLinkage(Mapping * mm) {
		if (!(mm->isUnspliced())) linkage = true;
		return;
	}
	bool linkable() {
		return linkage;
	}
	double coverage(int rl) {
		return double(rl) * totalReads / double(this->length());
	}
	double readCount() {
		return totalReads;
	}
	int allReadCount() {
		return allReads;
	}
	int start() {
		return from;
	}
	int end() {
		return to;
	}
	int length() {
		return to - from + 1;
	}
	int getIslandID() {
		return islandID; 
	}
	int getFlag() {
		return flag;
	}
	void raiseFlag(int fl) {
		if (fl > flag) flag = fl;
		return;
	}
	int distinctReads() {
		int dRead = 0;
		int last = 0;
		for(imap::const_iterator it = posmap.begin(); it != posmap.end(); ++it){
			if (it->first != last) dRead++;
			last = it->first;
		}
		return dRead;
	}
	Region * getRegion() {
		return region;
	}
	
	int bestGeneReads() {
		return bGeneReads;
	}
	void bestGene(Gene * g) {
		bGene = g;
		return;
	}
	Gene * getBestGene() {
		return bGene;
	}
	int bestGeneID() {
		int returner = 0;
		if (bGene != NULL) returner =  bGene->gid();
		return returner;
	}
	int uniqueReads () {
		int u = 0;
		for(imap::const_iterator it = posmap.begin(); it != posmap.end(); ++it) if (it->second->getRead()->isUnambiguous()) u++;
		return u;
	}
	int ambiguousReads () {
		int u = 0;
		for(imap::const_iterator it = posmap.begin(); it != posmap.end(); ++it) if (!(it->second->getRead()->isUnambiguous())) u++;
		return u;
	}
	
	char getRflagchar () {
		return rFlag;
	}
	void setRflag (char r) {
		rFlag = r;
		return;
	}
	int connexID() {
		return connex;
	}
	void connexID(int cid) {
		connex = cid;
	}
	int leftend() {
		return leftmost;
	}
	int rightend() {
		return rightmost;
	}
	
	//type tests
	double dominanceTest;
	double lengthTest;
	int multiexonTest;
	int splicedGeneMaps;
	
	friend class Region; // for island annotation build
protected:
	const static int idist    = ISLANDIST;
	typedef map<Gene*, int> gmap;
	
	gmap genemap;
	gmap gm; // temporary
	
	imap posmap;
	
	char rFlag;
	bool linkage; // true if island contains a spliced read
	double totalReads;
	int allReads;
	Region * region;
	int from;
	int to;
	int leftmost; // leftmost slice (even if not in island)
	int rightmost; // rightmost slice (even if not in island)
	int connex; // connex id (same ids are connected
	int strand; // 0:both 1:fwd -1:rev
	int flag;
	int islandID; // allows 1000000 islands
	int uniscore;
	Gene * bGene;
	int bGeneReads;
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


