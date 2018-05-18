/*
 *  smack-fe.hpp
 *  smack-fe
 *
 *  Created by David Brawand on 16.03.10.
 *  Copyright 2010 Université de Lausanne. All rights reserved.
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
#include <sstream>
//#include "universe.hpp"

using namespace std;

/* ------------------------------------------------------------- */
/* GLOBAL TYPEDEFS                                               */
/* ------------------------------------------------------------- */
typedef map<int, int> intmap;
typedef multimap<int, int> intmultimap;
typedef set<int> intset;
/* ------------------------------------------------------------- */
/* Class declarations                                            */
/* ------------------------------------------------------------- */
class Splice;
class Region;     // contains gene  (coveragemaps in exons)
class Read;
class Mapping;
/* ------------------------------------------------------------- */
/* GLOBAL TYPEDEFS                                               */
/* ------------------------------------------------------------- */
typedef map<int,Read*>          readmap;
typedef map<int,Region*>        regiomap;
typedef vector<Mapping*>        mapvec;
typedef multimap<int, Mapping*> imap;

typedef map<int, Splice*>       splicemap;
typedef vector<int>             intvec;
/* ------------------------------------------------------------- */
/* GLOBAL TYPEDEFS                                               */
/* ------------------------------------------------------------- */
/* ------------------------------------------------------------- */
/* Classes                                                       */
/* ------------------------------------------------------------- */

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
		slices.insert(pair <int,int> (from, to));
		if (intronid != 0) introns.push_back(intronid);
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
	
	intvec introns;
protected:
	int spliceID;
	int seqRegion;
	int strand;
	intmap slices;
};


class Mapping {
public:
	Mapping(){};
	Mapping(Read * rd, int regid, int str, intmap * slc, int mm, int mte, int * values){
		if (str > 1 || str < -1) abort();
		regionid = regid;
		strand   = str;
		slices   = slc;
		read     = rd;
		mism     = mm;
		accepted = (mte == 0) ? 1 : -1; // lonely reads are automatically accepted
		for (int ii = 0; ii != fieldcount; ++ii) vv[ii] = values[ii];
	};
	Mapping(Read * rd, int regid, int str, int pos, int readl, int mm, int mte, int * values){
		if (str > 1 || str < -1) abort();
		regionid = regid;
		strand   = str;
		slices   = new intmap;
		slices->insert(pair<int, int> (pos,pos+readl-1));
		read     = rd;
		mism     = mm;
		accepted = (mte == 0) ? 1 : -1; // lonely reads are automatically accepted
		for (int ii = 0; ii != fieldcount; ++ii) vv[ii] = values[ii];
	};
	~Mapping(){};
	
	bool operator==(const Mapping * other) const {
		if (regionid != other->regionid) return false;
		if (strand != other->strand) return false;
		for(intmap::const_iterator aa = slices->begin(); aa != slices->end(); ++aa) {
			if (other->slices->find(aa->first) == other->slices->end() || other->slices->find(aa->first)->second != aa->second) return false;
		}
		return true;
	}
	void show() {
		std::cerr << regionid << "\t" << strand << std::endl;
		for (intmap::const_iterator xx = slices->begin(); xx != slices->end(); ++xx) {
			std::cerr << "\t" << xx->first << "-" << xx->second << std::endl;
		}
		return;
	}
	string getString() {
		stringstream ret;
		ret << vv[0] << "\t" << vv[1] << "\t" << vv[2] << "\t" << vv[3] << "\t" << vv[4];
		return ret.str();
	}
	int len() {
		if (slices->size() == 0) {
			cerr << endl << "OOPS: Empty mapping ";
			abort();
		}
		int len = 0;
		for(intmap::const_iterator it = slices->begin(); it != slices->end(); ++it) {
			len += it->second - it->first + 1;
		}
		return len;
	}
	int region() {
		return regionid;
	}
	int from() {
		return slices->begin()->first;
	}
	int to() {
		intmap::iterator it = slices->end();
		--it;
		return it->second;
	}
	int str() {
		return strand;
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
	int posFromStart(int offset) {
		int offsetsum = 0;
		for (intmap::const_iterator sl = slices->begin(); sl != slices->end(); ++sl) {
			if      (offset <= offsetsum + sl->second - sl->first) return sl->first + (offset - offsetsum); 
			else if (offset >  offsetsum + sl->second - sl->first) offsetsum += (sl->second - sl->first + 1);
		}
		cerr << endl << "ERROR: Asked for position outside of mapping" << endl;
		cerr << offset << "\tin " << endl;
		for (intmap::const_iterator er = slices->begin(); er != slices->end(); ++er) cerr << "\t" << er->first << "-" << er->first << endl;
		abort();
	}
	int posFromStartContinued(int offset) {
		if (offset < this->len()) return this->posFromStart(offset);
		else                      return this->to() - this->len() + offset + 1;
	}
	
	bool isWithin(int pos) {
		for (intmap::const_iterator sl = slices->begin(); sl != slices->end(); ++sl) {
			if (sl->first <= pos && pos <= sl->second) return true;
		}
		return false;
	}
	void accept(int acc) {
		if (accepted != -1) cerr << "WARNING: Redecision on mapping! (" << accepted << "->" << acc << ")" << endl;
		accepted = acc;
	}
	bool isAccepted() {
		return (accepted > 0) ? true : false;
	}
	bool notAccepted() {
		return (accepted < 1) ? true : false;
	}
	bool decided() {
		return (accepted > -1) ? true : false;
	}
	int gaps() {
		return slices->size() - 1;
	}
	Read * getRead () {
		return read;
	}
	friend class Region; //OK
	friend class Read;   //OK
	int accepted;
protected:
	int vv[5];
	const static int fieldcount = 5;
	intmap * slices; // move to protected again
	Read * read;
	int	regionid;
	int strand;
	int mism; // number of mismatches
};

//OK
class Read {
public:
	Read(){};
	Read(int identifier, int mate) {
		readid = identifier;
		mateid = mate;
	};
	~Read(){
		// destruct all mappings
	}
	int unify() {
		int removedMappings = 0;
		if (mappings.size() == 1) return removedMappings;
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
		int lowest = 1000; // so high it must be lowered
		for (unsigned int i=0; i < mappings.size(); i++) if (lowest > mappings[i]->mism) lowest = mappings[i]->mism;  
		// => FILTER ALL MAPPINGS WITH MORE MISMATCHES
		for (unsigned int i=0; i < mappings.size(); i++) {
			if (mappings[i]->mism > lowest) {
				mappings.erase(mappings.begin() + i);
				removedMappings++;
				--i; // step one back as vector gets shortenend by previous operation
			}
		}
		return removedMappings;
	}
	
	int mid() {
		return mateid;
	}
	int rid() {
		return readid;
	}
	bool isUnambiguous() {
		if (mappings.size() == 1) return true;
		else return false;
	}
	int countMappings() {
		return mappings.size();
	}
	void addFullMapping(int r, int s, int p, int l, int mm, int mate, int * mval) { // region, strand, start, readlength
		mappings.push_back(new Mapping(this,r,s,p,l,mm,mate,mval));
	}	
	void addSpliceMapping(Splice * sp, int s, int p, int l, int mm, int mate, int * mval) { // region, strand, position, readlength, splice
		mappings.push_back(new Mapping(this,sp->region(),s,sp->resolvePosition(p,l),mm,mate,mval));
	}
	int readMappingSize() {
		return mappings.size();
	}
	Mapping * readMapping(int mapid) {
		return mappings[mapid];
	}
	bool contains(Mapping * mp) {
		for(unsigned int v = 0; v != mappings.size(); ++v) {
			if (mp == mappings[v]) return true;
		}
		return false;
	}
protected:
	int readid;
	int mateid;
	mapvec mappings;
};
class Region {
public:
	Region(int val, int nd, int md){
		regionid = val;
		acc = 0;
		mind = nd;
		maxd = md;
		singletons = 0;
		for (int i=0; i <= abs(md-nd); ++i) statvec.push_back(0);
	}
	~Region() {};
	int rid() {
		return regionid;
	}
	void addIntron(int left, int right){
		introns.insert(pair<int,int> (left,right));
		return;
	}
	
	void addMapping(int * val, int rl, int mate) {
		if (rmap.find(val[0]) == rmap.end()) rmap.insert(pair<int,Read*>(val[0], new Read(val[0], mate)));
		rmap.find(val[0])->second->addFullMapping(val[2], val[1], val[3], rl, val[4], mate, val);
	}
	
	void addMapping(int * val, int rl, int mate, Splice * spl) {
		if (rmap.find(val[0]) == rmap.end()) rmap.insert(pair<int,Read*>(val[0], new Read(val[0], mate)));
		rmap.find(val[0])->second->addSpliceMapping(spl, val[1], val[3], rl, val[4], mate, val);
	}
	
	void indexMappings() {
		mmm.clear();
		endmap.clear();
		for(readmap::const_iterator rd = rmap.begin(); rd != rmap.end(); ++rd) {
			for (int mp = 0; mp != rd->second->readMappingSize(); ++mp) {
				mmm.insert(pair<int,Mapping*>(rd->second->rid(), rd->second->readMapping(mp)));
				endmap.insert(pair<int,Mapping*>(rd->second->readMapping(mp)->to(), rd->second->readMapping(mp)));
			}
		}
		return;
	}
	
	void markInRangeMappings(int readl, intmap* mapcount) { /***************************/
		boost::progress_display show_progress( endmap.size() );
		
		int a = mind + readl; // min from left read's start
		int z = maxd + readl; // max from left read's start
		
//		int counter = 0;
		for (imap::const_iterator it = endmap.begin(); it != endmap.end(); ++it) {
			++show_progress;
			
			// quickskip
//			cerr << endl;
//			cerr << "\r                               \r" << ++counter ;
			
			if (it->second->decided()) continue; // skips mappings which have already been decided (acc flag 0 or 1)
		
			// MATE CHECK (counts a priori accepted reads)
			if (it->second->getRead()->mid() == 0) {
				singletons++;
				continue;
			}
		/*	
			// NO MATE MAPPING CASE (ACCELERATES)
			if (mmm.find(it->second->getRead()->mid()) == mmm.end()) { // MATE IS NOT MAPPED IN CHROMOSOME
				if (mapcount->find(it->second->getRead()->rid()) != mapcount->end() && mapcount->find(it->second->getRead()->rid())->second == 1) {
					// high qual with dead brother
					it->second->accept(1);
					acc++;
				}
				else {
					it->second->accept(0);
				}
				continue;
			}
		*/	
			// get accepted starts for mapping (using intron table)
//			cerr << "\tX";
			intmap possible;
//			cerr << "1";
			possible.clear();
//			cerr << "2" << endl;
//			cerr << ".>" << a << endl;
//			cerr << "A>" << it->first << endl;
//			cerr << "B>" << it->second->getRead()->rid() << endl;
//			cerr << "C>" << it->second->getRead()->mid() << endl;
//			cerr << ":>" << it->second->len() << endl;
//			cerr << "|>" << z << endl;
			// add immediate neighbor positions
			for(int d = a;  d < it->second->len(); d++) { // if overlap suggested by parameter
//				cerr << "=>" << d << "\t" << z-d << endl;
				possible.insert(pair<int,int> (it->second->posFromStart(d), z-d));
			}
//			cerr << "3";
			for(int d = it->second->len();  d <= z; d++) { // after read without introns
				possible.insert(pair<int,int> (it->second->posFromStartContinued(d), z-d));
			}
//			cerr << "4";
			// find possible positions using introns
			for (intmap::const_iterator pp = possible.find(it->second->posFromStartContinued(it->second->len())); pp != possible.end(); ++pp) {
				//cerr << endl << pp->first << "~" << pp->second << ":";
				//check if intron to add
				if (introns.find(pp->first) != introns.end()) {
					pair<intmultimap::iterator,intmultimap::iterator> ret = introns.equal_range(pp->first);
					for (intmultimap::iterator rr = ret.first; rr != ret.second; ++rr) {
						if (possible.find(rr->second + 1) != possible.end() && possible.find(rr->second + 1)->second < pp->second) { // exists already
							possible.insert(pair<int,int>(rr->second + 1, pp->second));
						}
						else  possible.insert(pair<int,int>(rr->second + 1, pp->second));
					}
				}
				// add following position if doesn't exist already
				else if (pp->second > 0) { 
					if (possible.find(pp->first + 1) == possible.end() || possible.find(pp->first + 1)->second < pp->second - 1) {
						// add followup
						possible.insert(pair<int,int>(pp->first + 1, pp->second - 1));
					}
				}
			}

			// get all paired mappings (first success is kept)
//			cerr << "\tZ";
			pair<imap::iterator,imap::iterator> ret = mmm.equal_range(it->second->getRead()->mid());
//			cerr << "\r(" << counter << ":1)";
//			cerr << "lower bound points to: " << ret.first->first << " => " << ret.first->second << "   upper bound points to: " << ret.second->first << " => " << ret.second->second << endl;
			for (imap::iterator rr = ret.first; rr != ret.second; ++rr) {
//				cerr << endl << "{{" << rr->first << "||";
//				cerr << endl << rr->second << "}}";
//				cerr << endl << rr->second->getRead()->rid() << ":";
//				cerr << rr->second->getRead()->rid();
//				cerr << endl << "<<" << rr->second->accepted << ">>";
//				rr->second->show();
//				cerr << "\r(" << counter << ":2.1)";
				if (rr->second->decided()) continue;
//				cerr << "\r(" << counter << ":2.2)";
				if (rr->second->str() == it->second->str()) continue;
//				cerr << "\r(" << counter << ":3)";
				if (possible.find(rr->second->from()) != possible.end()) { // find mate!
						// randomize accept
					if (( rand() % 2 ) == 0) {
						it->second->accept(1);
						rr->second->accept(0);
					}
					else {
						it->second->accept(0);
						rr->second->accept(1);
					}
					// count acc
					acc++;
//					cerr << "\r(" << counter << ":4)";
					statvec[maxd - mind - possible.find(rr->second->from())->second]++;
//					cerr << "\r(" << counter << ":5)";
					break;
				}
			}
//			cerr << endl << "\t_";
			if (!(it->second->decided())) { // nothing was accepted? not decided
				// check single 
				if (mapcount->find(it->second->getRead()->rid()) != mapcount->end() && mapcount->find(it->second->getRead()->rid())->second == 1) {
					// check brother (when lo qual then accept)
					if (mapcount->find(it->second->getRead()->mid()) != mapcount->end() || mapcount->find(it->second->getRead()->mid())->second > 1) {
						// high qual with dead/loqaul brother
						it->second->accept(1);
						acc++;
					}
				}
			}
		}
		return;
	}
	
	int totalReads() {
		return rmap.size();
	}
	int totalMappings() {
		if (!(mmm.size() == endmap.size())) {
			cerr << "FATAL ERROR: Mapping index size not consistent! (mmm != startmap || mmm != endmap)";
			abort();
		}
		return mmm.size();
	}
	int acceptedMappings() {
		return acc;
	}
	int singles() {
		return singletons;
	}
	int getDistanceMode(int sts) {
		int mode = 0;
		int max  = 0;
		for (int i=0; i < sts; ++i) {
			if (statvec[i] > max) {
				max = statvec[i];
				mode = i;
			}
		}
		return mode;
	}
	intvec getDistanceDistributon() {
		return statvec;
	}
	
	readmap rmap;
	imap mmm; // mapping multi map
	imap endmap; // mappings sorted by end
protected:
	int singletons; // reads without partner
	int acc;  // accepted mappings
	int maxd; // maximum distance
	int mind; // minimum distance
	intvec statvec;
	int regionid;
	intmultimap introns;
};


