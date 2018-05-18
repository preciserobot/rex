/*
 *  splicigar.hpp
 *  splicigar
 *
 *  Created by David Brawand on 10.03.10.
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

#define MINEXONSIZE 9;

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
class Region;
class Splicigar;
class Intron;
/* ------------------------------------------------------------- */
/* GLOBAL TYPEDEFS                                               */
/* ------------------------------------------------------------- */
typedef map<int,Region*> regiomap;
typedef multimap<int,Splicigar*> cigarmultimap;
typedef multimap<int,Intron*> intronmultimap;
typedef map<int,Intron*> intronmap;
typedef set<Splicigar*> splicigarset;
/* ------------------------------------------------------------- */
/* Classes                                                       */
/* ------------------------------------------------------------- */
class Intron {
public:
	Intron(int intronid, int first, int last, int str) {
		iid   = intronid;
		left  = first;
		right = last; 
		strand = str;
	}
	~Intron() {}
	int from() {
		return left;
	}
	int to() {
		return right;
	}
	int intronid() {
		return iid;
	}
	int intronstrand() {
		return strand;
	}
protected:
	int iid;
	int left;
	int right;
	int strand;
};

class Splicigar {
public:
	Splicigar() {}
	Splicigar(int givenid, Intron * central) {
		sid = givenid;
		strand = central->intronstrand();
		extremes[0] = central->from();
		extremes[1] = central->to();
		offsets[0] = 0;
		offsets[1] = 0;
		slices.insert(pair<int,Intron*>(0,central));
		startsat = 0;
		endsat   = 0;
	}
	~Splicigar() {}
	int cigarid() {
		return sid;
	}
	void cigarid(int cid) {
		sid = cid;
		return;
	}
	int cigarstrand() {
		return strand;
	}
	void starts(int d) {
		startsat = d - offsets[0] + extremes[0];
		return;
	}

	void ends(int d) {
		endsat = extremes[1] + d - offsets[1];
		return;
	}
	int starts() { return startsat; }
	int ends() { return endsat; }
	void addIntron(Intron * add) {
		// add intron to splices by calculation correct index (remainder)
		if (add->from() < extremes[0]) {
			offsets[0]  = offsets[0] + add->to() - extremes[0] + 1;
			slices.insert(pair<int,Intron*>(offsets[0],add));
			extremes[0] = add->from();
		}
		else if (add->to()   > extremes[1]) {
			offsets[1]  = offsets[1] + add->from() - extremes[1] - 1;
			slices.insert(pair<int,Intron*>(offsets[1],add));
			extremes[1] = add->to();
		}
		else {
			bool checker = false;
			for (intronmap::const_iterator it = slices.begin(); it != slices.end(); ++it) if (add->intronid() == it->second->intronid()) checker = true;
			if (!(checker)) {
				cerr << "I'd like to add " << add->intronid() << "\t" << add->from() << "\t" << add->to() << endl;
				for (intronmap::const_iterator it = slices.begin(); it != slices.end(); ++it) {
					cerr << "\t" << it->second->intronid() << "\t" << it->second->from() << "\t" << it->second->to() << endl;
				}
				abort();
			}
		}
		return;
	}
	void mergeleft(Splicigar * sc, int ovrh) {
		for (intronmap::const_iterator sl = sc->slices.begin(); sl != sc->slices.end(); ++sl) {
			this->addIntron(sl->second);
		}
		this->ends(ovrh);
		delete sc;
		return;
	}
	int cigarLength() {
		int len = 0;
		int lt = startsat;
		for(intronmap::const_iterator sl = slices.begin(); sl != slices.end(); ++sl){
			len += sl->second->from() - lt;
			lt = sl->second->to() + 1;
		}
		len += endsat - lt + 1;
		return len;
	}
	bool checkLength(int l) {
		int cigarlength = this->cigarLength();
		if (cigarlength == l) return true;
		else return false;
	}
	int centerIntron(int dis) {
		int len = 0;
		int lt = startsat;
		for(intronmap::const_iterator sl = slices.begin(); sl != slices.end(); ++sl){
			len += sl->second->from() - lt;
			if (len == dis) return sl->second->intronid();
			lt = sl->second->to() + 1;
		}
		abort();
	}
	
	intronmap slices; // bases left (-/+), intron
	int offsets[2];   // last<->first base in intron [1234]444444[567]=>
	int extremes[2];  // last/first base in intron
private:
	int sid;
	int strand;
	int startsat;
	int endsat;
};

class Region {
public:
	Region(int val, int dist) {
		regionid = val;
		overhang = dist;
	}
	~Region() {
		left.clear();
		right.clear();
	}
	int rid() {
		return regionid;
	}
	void addIntron(int * i) {
		Intron * intron = new Intron(i[0],i[3],i[4],i[2]);
		left.insert(pair<int,Intron*>(intron->from(),intron));
		right.insert(pair<int,Intron*>(intron->to(),intron));
		return;
	}
	void addExon(int * i) {
		intset * addto;
		if (i[3] == 1) addto = &exonicplus;
		else if (i[3] == -1) addto = &exonicminus;
		else abort();
		for (unsigned int j = i[1]; j <= i[2]; j++) addto->insert(j);
		return;
	}
	int makeCigars (int cigarnumber, bool stringent, int readlength) {
		// set cigar storage
		// activate
		if (stringent) activecigars = &ssc;
		else activecigars = &sc;
		
		// should crosscheck each if non-exonic is at same time an intron 
		int initial;
		splicigarset newcigars;
		pair<cigarmultimap::iterator,cigarmultimap::iterator> ret;
		for (intronmultimap::const_iterator it = left.begin(); it != left.end(); ++it) {
			
			//find closest splice sites and evaluate distance
			int upNext = readlength, downNext = readlength;
			for (int i = 0; i < readlength; ++i) {
				if (right.find(it->second->from() - (1 + i)) != right.end()) {
					upNext = i;
					break;
				}
			}
			for (int i = 0; i < readlength; ++i) {
				if (left.find(it->second->to() + (1 + i)) != left.end()) {
					downNext = i;
					break;
				}
			}
			//bool avoidRedundancy = (upNext+downNext >= readlength) ? true : false; // will be set to TRUE if closest neigboring introns are to far to be different
			//cout << upNext << "<->" << downNext << "\t" << (noRedundancy ? "1" : "0") << endl;

			// select exon strand (by starting intron)
			intset * exoncover;
			if (it->second->intronstrand() == 1) exoncover = &exonicplus;
			else if (it->second->intronstrand() == -1) exoncover = &exonicminus;
			else {
				std::cerr << "ERROR: Strand for Intron " << it->second->intronid() << " (" << it->second->from();
				std::cerr << "-" << it->second->to() << ") is " << it->second->intronstrand() << std::endl;
				abort();
			}
			Splicigar * tmpcigar = new Splicigar(++cigarnumber, it->second);
			int tmpintronid = it->second->intronid();
			
			activecigars->insert(pair<int,Splicigar*> (tmpintronid, tmpcigar));
			//UPSTREAM
			do {
				initial = activecigars->count(it->second->intronid());
				ret = activecigars->equal_range(it->second->intronid());
				for (cigarmultimap::const_iterator cc = ret.first; cc != ret.second; ++cc) {
					if (cc->second->starts() > 0) continue; // has already come to an end (eg has finished)
					// go left on same strand
					for (int t = cc->second->offsets[0]; t >= -1*overhang; t--) {
						int check = cc->second->extremes[0] + t - cc->second->offsets[0] - 1;
						//std::cerr << "{" << t << "|" << check << "}";
						if (t == -1*overhang) cc->second->starts(-1*overhang);// set done by setting the start coordinate
						else if ((t <= cc->second->offsets[0] - minex) && (right.find(check) != right.end()) && (right.find(check)->second->intronstrand() == it->second->intronstrand())) {
							// found intron
							if (stringent && exoncover->find(check) == exoncover->end()) { // no exon at this place
								// just add intron, do not clone as intron is mandatory (no more exon after intron position)
								cc->second->addIntron(right.find(check)->second); //add intron to current
								--initial; //decrement initial count to be sure that a further iteration will be done nevertheless the actual cigar count has not changed yet
							}
							else {
								// has intron and is exonic -> clone
								Splicigar * copy = new Splicigar(* cc->second); //clone
								copy->addIntron(right.find(check)->second); //add intron
								copy->cigarid(++cigarnumber); // set new cigarid
								newcigars.insert(copy);//store
							}
						}
					}
				}
				// merge with new ssc
				for (splicigarset::iterator nc = newcigars.begin(); nc != newcigars.end(); ++nc) {
					activecigars->insert(pair<int,Splicigar*>(it->second->intronid(), *nc));
				}
				newcigars.clear();
			} while (initial != activecigars->count(it->second->intronid())); // if initial decrement has not desired effect replace by check if all are finished (also below, would be much less time-efficient

			
			//DOWNSTREAM WITH MERGE
			//if (avoidRedundancy) {
			//} 
			if (false) {
				int dummycigarnumber = 0;
				cigarmultimap partcigars;
				partcigars.insert(pair<int,Splicigar*> (it->second->intronid(), new Splicigar(++dummycigarnumber, it->second)));
				do {
					initial = partcigars.count(it->second->intronid());
					ret = partcigars.equal_range(it->second->intronid());
					for (cigarmultimap::const_iterator cc = ret.first; cc != ret.second; ++cc) {
						if (cc->second->ends() > 0) continue; // has already come to an end
						// go right on same strand
						for (int t = cc->second->offsets[1]; t <= overhang; t++) {
							int check = cc->second->extremes[1]+t-cc->second->offsets[1]+1;

							// IS DONE?
							if (t == overhang) cc->second->ends(overhang);// set done by setting the end coordinate
							// IS THERE AN INTRON?
							else if ((t >= cc->second->offsets[1] + minex) && (left.find(check) != left.end()) && (left.find(check)->second->intronstrand() == it->second->intronstrand())) {
								// 
								if (stringent && exoncover->find(check) == exoncover->end()) {
									// just add intron, do not clone as intron is mandatory (no more exon after intron position)
									cc->second->addIntron(left.find(check)->second); //add intron to current
									--initial; //decrement initial count to be sure that a further iteration will be done nevertheless the actual cigar count has not changed yet
								}
								else {
									Splicigar * copy = new Splicigar(* cc->second); //clone
									copy->addIntron(left.find(check)->second); //add intron
									copy->cigarid(++cigarnumber); // set new cigarid
									newcigars.insert(copy);//store
								}
							}
						}
					}
					// merge with new ssc
					for (splicigarset::iterator nc = newcigars.begin(); nc != newcigars.end(); ++nc) {
						partcigars.insert(pair<int,Splicigar*>(it->second->intronid(), *nc));
					}
					newcigars.clear();
				} while (initial != partcigars.count(it->second->intronid())); // if initial decrement has not desired effect replace by check if all are finished (also below, would be much less time-efficient
				
				// do an intelligent merge
				if (activecigars->count(it->second->intronid()) > partcigars.count(it->second->intronid())) {
					// multiply partcigars
					do {
						Splicigar * copy = new Splicigar(* partcigars.begin()->second); //clone
						copy->cigarid(++dummycigarnumber); // set new cigarid
						partcigars.insert(pair<int,Splicigar*> (it->second->intronid(), new Splicigar(* partcigars.begin()->second))); // multiply first element to get equal numbers
					} while (activecigars->count(it->second->intronid()) > partcigars.count(it->second->intronid()));
				}
				else if (activecigars->count(it->second->intronid()) < partcigars.count(it->second->intronid())) {
					// multiply activecigars
					do {
						Splicigar * copy = new Splicigar(* activecigars->begin()->second); //clone
						copy->cigarid(++cigarnumber); // set new cigarid
						activecigars->insert(pair<int,Splicigar*> (it->second->intronid(), new Splicigar(* activecigars->begin()->second))); // multiply first element to get equal numbers
					} while (activecigars->count(it->second->intronid()) < partcigars.count(it->second->intronid()));
				}
				if (activecigars->count(it->second->intronid()) != partcigars.count(it->second->intronid())) {
					cerr << "ERROR: There's a big problem (unequal sizes)!" << endl;
				}
				
				// do the merge
				pair<cigarmultimap::iterator,cigarmultimap::iterator> partc   = partcigars.equal_range(it->second->intronid());
				pair<cigarmultimap::iterator,cigarmultimap::iterator> activec = activecigars->equal_range(it->second->intronid());
				cigarmultimap::const_iterator dd = activec.first;
				for (cigarmultimap::const_iterator cc = partc.first; cc != partc.second; ++cc) {
					dd->second->mergeleft(cc->second, overhang); // also frees memory from argument (calls delete)
					// check if end concurrently and also increment dd at every iteration
					if (++dd == activec.second) {
						if (++cc != partc.second) {
							cerr << "ERROR: There's an even bigger problem (unbalanced merge)!" << endl;
							abort();
						}
						break; // okay it ends (needed as otherwise cc goes past partc.second
					}
					//if (dd->second->slices.size() > 3) cerr << "(" << dd->second->centerIntron(overhang) << ")";
				}
				partcigars.clear();
			}
			//DOWNSTREAM WITH REDUNDANCY (A BIT SLOWER BUT NOT RELEVANT)
			else {
				do {
					initial = activecigars->count(it->second->intronid());
					ret = activecigars->equal_range(it->second->intronid());
					for (cigarmultimap::const_iterator cc = ret.first; cc != ret.second; ++cc) {
						if (cc->second->ends() > 0) continue; // has already come to an end
						// go right on same strand
						for (int t = cc->second->offsets[1]; t <= overhang; t++) {
							int check = cc->second->extremes[1]+t-cc->second->offsets[1]+1;
							// IS DONE?
							if (t == overhang) cc->second->ends(overhang);// set done by setting the end coordinate
							// IS THERE AN INTRON?
							else if ((t >= cc->second->offsets[1] + minex) && (left.find(check) != left.end()) && (left.find(check)->second->intronstrand() == it->second->intronstrand())) {
								// 
								if (stringent && exoncover->find(check) == exoncover->end()) {
									// just add intron, do not clone as intron is mandatory (no more exon after intron position)
									cc->second->addIntron(left.find(check)->second); //add intron to current
									--initial; //decrement initial count to be sure that a further iteration will be done nevertheless the actual cigar count has not changed yet
								}
								else {
									Splicigar * copy = new Splicigar(* cc->second); //clone
									copy->addIntron(left.find(check)->second); //add intron
									copy->cigarid(++cigarnumber); // set new cigarid
									newcigars.insert(copy);//store
								}
							}
						}
					}
					// merge with new ssc
					for (splicigarset::iterator nc = newcigars.begin(); nc != newcigars.end(); ++nc) {
						activecigars->insert(pair<int,Splicigar*>(it->second->intronid(), *nc));
					}
					newcigars.clear();
				} while (initial != activecigars->count(it->second->intronid()));
			}
		}
		
		
		// stats && multipslices
		for (cigarmultimap::const_iterator aa = activecigars->begin(); aa != activecigars->end(); ++aa) {
			if (aa->second->checkLength(2*overhang)) {
				stats[aa->second->slices.size()]++;
			}
			else {
				cerr << "OOPSIE -> " << aa->second->cigarLength()  << " != " << (2*overhang) << endl;
				cerr << "\t\t" << aa->second->starts() << endl;
				for (intronmap::const_iterator it = aa->second->slices.begin(); it != aa->second->slices.end(); ++it) {
					cerr << "\t" << it->second->intronid() << "\t" << it->second->from() << "\t" << it->second->to() << endl;
				}
				cerr << "\t\t" << aa->second->ends() << endl;
				abort();
			}
		}
		return cigarnumber;
	}
	cigarmultimap * activecigars;
	cigarmultimap sc; // startingintron,cigar
	cigarmultimap ssc; // startingintron,stringentcigar
	intmap stats; // count how many junctions are in cigar
protected:
	int regionid;
	int overhang;
	intronmultimap left;
	intronmultimap right;
	intset exonicplus;
	intset exonicminus;
	const static int minex = MINEXONSIZE;
};
