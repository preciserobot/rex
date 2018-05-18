/*
 *  geco.hpp
 *  geco - the gene coverage evaluator
 *
 *  Created by David Brawand on 03.02.10.
 *  Copyright 2010 UNIL. All rights reserved.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <assert.h>
#include <math.h>

using namespace std ;

#define READLENGTH 76 // size of buckets in which we sort genes (by start coordinates)

/* ------------------------------------------------------------- */
/* GLOBAL TYPEDEFS                                               */
/* ------------------------------------------------------------- */
typedef map<int, int> covmap;
/* ------------------------------------------------------------- */
/* Class declarations                                            */
/* ------------------------------------------------------------- */
class Exon;       // Exons (no flattening needed)
class Transcript; // trancript refers to exons
class Gene;       // contains transcripts and all exons
class Bucket;     // a Gene Bucket making Regions "smaller" by splitting into junks
class Region;     // contains gene  (coveragemaps in exons)
class Linear;     // Linear model
typedef map<int, Exon*> exonmap;
typedef map<int, Gene*> genemap;
typedef multimap<int, Gene*> multigenemap;
typedef map<int, Bucket*> bucketmap;
typedef map<int, Transcript*> transcriptmap;
typedef set<int> intset;
typedef vector<int> intvec;
typedef vector<double> doublevec;
const static int rl = READLENGTH;

class Exon {
public:
	Exon(int i, int j, int k, int l) {
		exonid = i;
		first  = j;
		last   = k;
		brin   = l;
		leng   = k-j+1;
		ecov   = new int [leng];
		covinfo = 0; // the number of bases that have coverage info
	}
	~Exon() { 
		delete[] ecov;
	};
	int eid() {
		return exonid;
	}
	int length() {
		return leng;
	}
	int strand() {
		return brin;
	}
	int start() {
		return first;
	}
	int end() {
		return last;
	}
	int noInfoBases() {
		return leng - covinfo;
	}
	int coverageAtPosition (int p) {
		if (p >= leng) { cerr << "EXONOOPSIE (coverageAtPosition)"; exit(1); }
		return ecov[p];
	} 

	void exonCoverageAdd(int f, int t, int o);

protected:
	int exonid;
	int first;
	int last;
	int brin;
	int leng;
	int covinfo;
	intset ci;
	int * ecov; // coverage array
};

class Transcript {
public:
	Transcript(int val, int str, int trg, int tge) {
		transcriptid   = val;
		direction      = str;
		tregion        = trg;
		tgene          = tge;
	}
	~Transcript(){ /*Exons are cleared in gene destructor*/ };
	int tid() {
		return transcriptid;
	}
	int exonCount () {
		return rankedexons.size();
	}
	bool isReverse() {
		if (direction < 0) return true;
		else return false;
	}
	bool addExon(Exon * e);
	void evalPrintCoverage(int fc);

	friend class Gene;
	exonmap rankedexons;
protected:
	int transcriptid;
	int direction;
	int tregion;
	int tgene;
	Linear * linmod;
};

class Gene {
public:
	Gene(int val, int s, Region * rp) {
		geneid   = val;
		from     = -1;
		to       = -1;
		genestr  = s;
		regionpointer = rp;
	}
	~Gene(){
		transcripts.clear();
		exons.clear();
	};
	int gid() {
		return geneid;
	}
	int transcriptCount () {
		return transcripts.size();
	}
	int exonCount () {
		return exons.size();
	}
	int getStart() {
		return from;
	}
	int geneLength() {
		return to - from + 1;
	}
	bool addExon(int val[]);
	bool overlaps(int a, int z);
	bool evalAddCoverage(int a, int b, int c);
	void evalPrintCoverage(int fc, ofstream * fout);
	int getFirstBucketNumber();
	int getSecondBucketNumber();
	
	transcriptmap transcripts;
protected:
	int geneid;
	int from;
	int to;
	int genestr;
	exonmap exons;
	Region * regionpointer;
};

class Bucket {
public:
	Bucket(int rid, int bucketnumber){
		region = rid;
		bucket = bucketnumber;
	}
	~Bucket(){} // destructor must be empty to preserve genes
	void addGene (Gene * g) {
		int initialSize = bucketgenes.size();
		bucketgenes.insert(pair<int,Gene*>(g->gid(),g));
		if (bucketgenes.size() == initialSize) {
			cerr << "FATAL: Gene was not added to bucket" << endl;
			exit(1);
		}
		return;
	}
	int geneCount() {
		return bucketgenes.size();
	}
	friend class Region;
protected:
	int region, bucket;
	genemap bucketgenes;
};

class Region {
public:
	Region(int val){
		regionid      = val;
		annotationend = 0;
		coverageend   = 0;
		bsize         = 0;
	}
	~Region() {
		genes.clear();
		delete finishedGene;
	};
	int rid() {
		return regionid;
	}
	int geneCount () {
		return genes.size();
	}
	int annotend () {
		return annotationend;
	}
	void setGeneCount() {
		finishedGeneNumber = 0;
		finishedGene = new int[genes.size()];
		return;
	}
	void addFinishedGene(int fg) {
		finishedGene[finishedGeneNumber++] = fg;
		if (finishedGeneNumber > genes.size()) {
			for (int i = 0; i < finishedGeneNumber; i++) std::cerr << endl << regionid << "\t F " << finishedGene[i];
			for (genemap::const_iterator it = genes.begin(); it != genes.end(); ++it) {
				std::cerr << endl << regionid << "\t G " << it->first;
			}
			std::cerr << std::endl << "FATAL: More genes finished than available in region (" << regionid << ")?!? ";
			std::cerr << finishedGeneNumber << " > " << genes.size() << std::endl;
			exit(1);
		}
		return;
	}
	int finishedGeneCount() {
		return finishedGeneNumber;
	}
	int getFinishedGene (int g) {
		return finishedGene[g];
	}
	int getBucketSize() {
		return bsize;
	}
	bool isFollowingSlice(int covstart, int covend) {
		if (covend == coverageend) {
			return false;
		}
		else if (covend >= coverageend && covstart - 1 == coverageend) return true; 
		else {
			cerr << "ERROR: Coverage File is not ordered" << endl;
			exit(1);
		}
		return false;
	}
	void bucketizeGenes() {
		// select bucketsize by gene length
		bsize = 0;
		for (genemap::const_iterator it = genes.begin(); it != genes.end(); ++it) if (it->second->geneLength() > bsize) bsize = it->second->geneLength(); 
		// print bucketsize that was selected;
		//std::cerr << "REGION " << regionid << " has BucketSize of " << bsize << std::endl;
		//sort genes by start
		int bnum = 0; // bucket number
		for (genemap::const_iterator it = genes.begin(); it != genes.end(); ++it) {
			bnum = it->second->getFirstBucketNumber();
			if (genebuckets.find(bnum) == genebuckets.end()) {
				genebuckets.insert(pair<int, Bucket*> (bnum, new Bucket(regionid, bnum)  ));// create bucket of not exists
			}
			genebuckets.find(bnum)->second->addGene(it->second);
			/*
			if (bnum != it->second->getSecondBucketNumber()) {
				bnum = it->second->getSecondBucketNumber();
				if (genebuckets.find(bnum) == genebuckets.end()) {
					genebuckets.insert(pair<int, Bucket*> (bnum, new Bucket(regionid, bnum)  ));// create bucket of not exists
				}
				genebuckets.find(bnum)->second->addGene(it->second);
			}
			*/
		}
		return;
	}
	
	bool addExon(int val[]);
	intset addCoverageSlice(int val[]);
	intset addCoverageSliceUsingBuckets(int val[]);
	void evaluateAndPrint(int fc, ofstream * output); // takes fullCoverage for RPKM calc

	genemap genes;
protected:
	int bsize;
	int regionid;
	int annotationend;
	int coverageend;
	int finishedGeneNumber;
	int * finishedGene;
	bucketmap genebuckets;
};

class Linear {
public:		
	Linear(doublevec dv, int ttt) {
		
		// calculate the averages of arrays x and y
		covec = dv;
		double xa = 0, ya = 0;
		for (int i = 0; i < dv.size(); ++i) {
			xa += i;     //x[i]
			ya += dv[i]; //y[i]
		}
		xa /= dv.size();
		ya /= dv.size();
		// calculate auxiliary sums
		double xx = 0, yy = 0, xy = 0;
		for (int i = 0; i < dv.size(); ++i) {
			double tmpx = i - xa;
			double tmpy = dv[i] - ya;
			xx += tmpx * tmpx;
			yy += tmpy * tmpy;
			xy += tmpx * tmpy;
		}
		// make sure slope is not infinite
		if (fabs(xx) == 0) {
			//cerr << endl << "ERROR: vector has size " << dv.size() << " in transcript " << ttt << endl;
			//exit(1);
			//cerr << "WARNING: transcript has insufficient length -> setting slooe to 0" << endl;
			m_b = 0; //slope
			m_a = dv[0]; //intercept
			m_coeff = 1; //regression_coeff
		}
		else {
			assert(fabs(xx) != 0);
			// calculate regression line parameters
			m_b = xy / xx;
			m_a = ya - m_b * xa;
			m_coeff = (fabs(yy) == 0) ? 1 : xy / sqrt(xx * yy);
		}
	}
	~Linear(){};
	
	
	////*********////
	// find best fit transcripts for each gene (iterate with rapram???
	// make global statistic in rnadb
	////*********////
	
	double avgCov() {
		int covsum = 0;
		for (int p = 0; p < covec.size(); p++) covsum += covec[p];
		return (double)covsum / (double)covec.size();
	}
	int modCov() {
		int peak = 0;
		for (int p = covec.size() - 2; p >= 0; p--) { // skip last base
			if (peak < covec[p]) {
				peak = covec[p];
			}
		}
		return peak;
	}
	double getValue(double x) {
		return m_a + m_b * x;
	}
	//! Returns the slope of the regression line
	double getSlope() {
		return m_b;
	}
	//! Returns the intercept on the Y axis of the regression line
	double getIntercept() {
		return m_a;
	}
	//! Returns the linear regression coefficient
	double getCoefficient() {
		return m_coeff;
	}
protected:
	double m_a, m_b, m_coeff;
	doublevec covec;
};

/* ------------------------------------------------------------- */
/* Function declarations (to avoid forward declaration)          */
/* ------------------------------------------------------------- */

//OK
// 0       1     2           3      4      5       6               7
// region, gene, transcript, exon,  start, end,    strand,         rank
// gene,   exon, region,     start, end,   strand, tag/transcript, score
bool Region::addExon(int val[]){
	if (genes.find(val[0]) == genes.end()) {
		Gene * newGene = new Gene(val[0], val[5], this);
		genes.insert(pair<int,Gene*>(val[0], newGene));
	}
	if (genes.find(val[0]) != genes.end()) {
		if (annotationend < val[4]) annotationend = val[4];
		return genes.find(val[0])->second->addExon(val);
	}
	return false;
}
bool Gene::addExon(int val[]){
	// create exon
	Exon * ex = new Exon(val[1], val[3], val[4], val[5]);

	// add to global exon map if not already there
	if (exons.find(val[1]) == exons.end()) exons.insert(pair<int,Exon*>(val[1], ex));

	// update gene coordintes
	if (from < 0 || from > val[3]) from = val[3];
	if (to   < 0 || to   < val[4]) to   = val[4];
	
	// add exon to transcript using rank
	bool added = false;
	
	// add to all transcripts with inferior tag
	for (unsigned int i = 0; i <=  val[6]; ++i) {
		// create transcript if necessary
		if (transcripts.find(i) == transcripts.end()) {
			// create new transcript
			Transcript * tr = new Transcript(i, val[5], val[2], val[0]);
			transcripts.insert(pair<int,Transcript*>(i, tr));
		}
		// get transcript
		if (transcripts.find(i)->second->addExon(ex)) added = true;
	}
	return added;
}
bool Transcript::addExon(Exon * e) {
	// add to exon map by start
	rankedexons.insert(pair<int,Exon*> (e->start(),e));
	return true;
}



//coverage add
intset Region::addCoverageSliceUsingBuckets(int val[]) {
	if (coverageend < val[2]) coverageend = val[2];
	intset affectedGenes;
	
	int startbucket = (val[1]/bsize)-1;
	int endbucket   = val[2]/bsize;
	
	for(int bu = startbucket; bu <= endbucket; ++bu) {
		if (genebuckets.find(bu) != genebuckets.end()) {
			Bucket * gb = genebuckets.find(bu)->second;
			for(genemap::const_iterator it = gb->bucketgenes.begin(); it != gb->bucketgenes.end(); ++it) {
				if (it->second->evalAddCoverage(val[1], val[2], val[4])) {
					affectedGenes.insert(it->second->gid());
					//std::cout << regionid << "\t" << 2 << "\t" << val[1] << "\t" << val[2] << "\t" << val[4] << "\t" << it->second->gid() << endl;
				}
			}
		}
		else {
			//cerr << "DEBUG: GeneBucket (" << *bu << ") does not exist in region " << regionid << endl;
			//exit(1);
		}
	}
	return affectedGenes;
}

bool Gene::evalAddCoverage(int a, int b, int c) {
	if ( (from <= a && a <= to)  || (from <= b && b <= to) || (a < from  && to < b ) ) {
		for(exonmap::const_iterator it = exons.begin(); it != exons.end(); ++it) it->second->exonCoverageAdd(a,b,c);
		return true;
	}
	return false;
}

void Exon::exonCoverageAdd(int f, int t, int o) {
	if (t < first || last < f) return; // quickskip
	int a=0, b=0;
	if (first <= f && f <= last) {
		a = f - first;
		b = (t > last) ? (last - first) : (t - first);
		for (int i = a; i <= b; i++) {
			ecov[i] = o;
			ci.insert(i);
		}
		covinfo += b - a + 1;
		
	}	
	else if (f < first && t >= first) {
		a = 0;
		b = (t > last) ? (last - first) : (t - first); 
		for (int i = a; i <= b; i++){
			ecov[i] = o;
			ci.insert(i);
		}
		covinfo += b - a + 1;
	}
	if (covinfo > leng) {
		cerr << "WARNING: (ex) " << covinfo << ">" << leng << "? " << exonid << endl;
		//exit(1);
	}
	return;
}

// degradation check 	// direct print to stdout (cout)
void Region::evaluateAndPrint(int fc, ofstream * output) {
	if (coverageend < annotationend) {
		//std::cerr << "\rNOTICE: Annotation (" << annotationend << ") exceeds coverage (" << coverageend << ") in region " << regionid << std::endl;
		int endcoverage[5] = {regionid, coverageend + 1, annotationend, (annotationend - coverageend + 2), 0};
		//std::cerr << "  -> WILL ADD " << endcoverage[0] << " " << endcoverage[1] << " " << endcoverage[2] << " " << endcoverage[3] << " " << endcoverage[4] << endl;
		this->addCoverageSliceUsingBuckets(endcoverage);
		//std::cerr << "  ->    ADDED " << endcoverage[0] << " " << endcoverage[1] << " " << endcoverage[2] << " " << endcoverage[3] << " " << endcoverage[4] << endl;
	}
	// foreach gene calc
	for(genemap::const_iterator it = genes.begin(); it != genes.end(); ++it) it->second->evalPrintCoverage(fc, output);
	return;
}

void Gene::evalPrintCoverage(int fc, ofstream * fout) {
	doublevec tcov;
	for(transcriptmap::const_iterator tr = transcripts.begin(); tr != transcripts.end(); ++tr) {
		// build coverage string
		for (exonmap::const_iterator ex = tr->second->rankedexons.begin(); ex != tr->second->rankedexons.end() ; ++ex) {
			if (ex->second->noInfoBases()) {
				// checks if exon has info from coveragefile
				cerr << "EXON_noinfo_" << ex->second->noInfoBases() << "_" << ex->second->length() << "_" << ex->second->eid() << endl;
				exit(1);
			}
			// reverse if necessary
			if (tr->second->isReverse()) {
				for (int z = ex->second->length() - 1; z >= 0 ; z--) {
					tcov.push_back(log2(ex->second->coverageAtPosition(z)+1));
				}
			} else {
				for (int z = 0; z < ex->second->length(); z++) {
					tcov.push_back(log2(ex->second->coverageAtPosition(z)+1));
				}
			}
		}
		int tiddy = tr->second->tgene;
		tr->second->linmod = new Linear(tcov, tiddy);
		*fout << tr->second->tregion << "\t";
		*fout << tr->second->tgene << "\t";
		*fout << tr->second->transcriptid << "\t";
		*fout << tcov.size() << "\t";
		*fout << tr->second->linmod->avgCov() << "\t";
		*fout << tr->second->linmod->modCov() << "\t";
		*fout << log2(1000 * tr->second->linmod->avgCov() / (0.000001 * fc / rl)) << "\t"; //// RPKM based on coverage-----------------------<<<<<<<<<<<<<<<<<<<<<<<<< CONTINUE HERE
		*fout << tr->second->linmod->getIntercept() << "\t";
		*fout << tr->second->linmod->getSlope() << "\t";
		*fout << tr->second->linmod->getCoefficient() << std::endl;
		tcov.clear();
	}
	return;
}


// AVOIDING USAGE OF INCOMPLETE TYPES

int Gene::getFirstBucketNumber() {
	return from / regionpointer->getBucketSize();
}
int Gene::getSecondBucketNumber() {
	return to / regionpointer->getBucketSize();
}
