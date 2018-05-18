/*
 *  universe.hpp
 *  smack-ms
 *
 *  Created by David Brawand on 06.05.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

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

/* ------------------------------------------------------------- */
/* Class declarations                                            */
/* ------------------------------------------------------------- */
class Region;
class Mapping;
class Splice;
/* ------------------------------------------------------------- */
/* Type declarations                                             */
/* ------------------------------------------------------------- */
typedef map<int, int> intmap;
typedef map<int, float> intfloatmap;
typedef vector<Mapping*> mapvec;
typedef vector<int> intvec;
typedef multimap<int, Mapping*> imap;
typedef multimap<int, int> intmultimap;
typedef set<Mapping*> mapset;


//LEGACY (OLD SMACK)
class Mapping {
public:
	Mapping(int m, int val[]) {
		readid  = val[0];
		region  = val[1];
		start1  = val[2];
		length1 = val[3];
		start2  = val[4];
		length2 = val[5];
		strand  = val[6];
		acceptance = -1;
		mateid = m;
	}
	int from() {
		return start1;
	}
	int to() {
		if (length2 > 0) return start2+length2-1;
		return start1+length1-1;
	}
	int posFromStart(int pos) {
		// pos is an offset!
		
		if (pos < length1 + length2) {
			if (pos < length1) return start1 + pos;
			else               return pos - length1 + start2;
		}
		else {
			cerr << "ERROR: Asked for position outside of mapping" << endl;
			cerr << pos << "\tin " << start1 << " " << length1 << " " << start2 << " " << length2 << endl;
			exit(1);
		}
	}
	int posFromStartContinued(int pos) {
		if (pos < length1+length2) return this->posFromStart(pos);
		else                       return this->to() - length1 - length2 + pos + 1;
	}
	
	
	
	int offsetPos (int off) { // returns the offset from an absolute position
		if (off < start1) {
			cerr << "ERROR: Offset out of range (before the start)!" << endl;
			exit(1);
		}
		else if (off <= start1 + length1 - 1) {
			return off - start1;
		}
		else if (off <= start2 + length2 - 1) {
			return off - start1;
		}
		cerr << "ERROR: Offset out of range (after the end)!" << endl;
		exit(1);
	}
	bool isWithin(int pos) {
		if (start1 <= pos && pos <= start1+length1-1) return true;
		if (start2 <= pos && pos <= start2+length2-1) return true;
		return false;
	}
	bool seamlessCover(Mapping * mp) {
		int offset = this->offsetPos(mp->from());
		for (int i = offset; i < this->len(); i++) {
			if (mp->posFromStart(i-offset) != this->posFromStart(i)) return false;
		}
		return true;
	}
	bool equals(Mapping * m) {
		if (this->rid() == m->rid() && this->s1() == m->s1() && this->s2() == m->s2() && this->l1() == m->l1() && this->l2() == m->l2() && this->str() == m->str() ) {
			return true;
		}
		else {
			return false;
		}
	}
	int rid()         { return readid; }
	int mid()         { return mateid; }
	int regionid()    { return region; }
	int len()         { return length1+length2; }
	int accept(int i) { return acceptance = i; }
	int accepted()    { return acceptance; }
	int mate()        { return mateid; }
protected:
	int readid;
	int mateid;
	int region;
	int start1;
	int start2;
	int length1;
	int length2;
	int strand;
	int acceptance;
};

