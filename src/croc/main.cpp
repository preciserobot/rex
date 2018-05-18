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
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/incremental_components.hpp> 
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>

using namespace std;

/* CROCO - Complete Redundant Orthology Clusters */
namespace boost {
	enum vertex_component_t { vertex_component = 111 };
	BOOST_INSTALL_PROPERTY(vertex, component);
}
using namespace boost;

typedef map<string, int> GeneID;
typedef map<string, int> stringcounts;
typedef map<int, string> invGeneID;
typedef set<string> setstring;
typedef set<int> setint;
template <typename ComponentMap>
struct vertexComponent {
	
	vertexComponent() {}
	
	vertexComponent(ComponentMap component, int f_component) : m_component(component), m_f_component(f_component) {}
	
	template <typename Vertex>
	bool operator()(const Vertex& v) const {
		return (get(m_component, v) == m_f_component);
	}
	
	ComponentMap m_component;
	int m_f_component;
};

/* MAIN */
int main(int argc, char* argv[]) {
	/* USAGE */
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " genelinks (results to STDOUT)" << std::endl;
		exit(1);
	}
	
	/* DECLARATIONS */
	char      line[200];
	char      scanL[20];
	char      scanR[20];
	
	int lico = 0;
	GeneID ids;
	invGeneID name;
	int runningid = 0;
	
	typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_component_t, int> > Graph;
	
	typedef property_map<Graph, vertex_component_t>::type ComponentMap;
	typedef filtered_graph<Graph, keep_all,	vertexComponent<ComponentMap> > FilteredGraph;
	
	graph_traits < Graph >::vertex_iterator vi, vi_end;
	graph_traits < Graph >::out_edge_iterator oei, oei_end;
	graph_traits < FilteredGraph >::vertex_iterator fvi, fvi_end;
	graph_traits < FilteredGraph >::out_edge_iterator foei, foei_end;
	
	Graph * g;
	FilteredGraph * fg;
	
	g = new Graph;
	
	// READ LINK FILE and add edges
	std::ifstream linkFile(argv[1]);
	if (linkFile.fail()) {
		std::cerr << "Error: could not read from file " << argv[1] << std::endl;
		exit(1);
	}
	std::cerr << "Reading gene links " << argv[1] << std::endl;
	while (!linkFile.eof()) {
		linkFile.getline(line, 200);
		if (linkFile.eof()) break;
		
		if (++lico % 1000 == 0) std::cerr << "\r" << lico;
		sscanf(line, "%s%s", scanL, scanR);
		string scanLstring = string(scanL);
		string scanRstring = string(scanR);
		
		if (ids.find(scanLstring) == ids.end()) {
			ids.insert(pair<string,int>(scanLstring,++runningid));
			name.insert(pair<int,string>(runningid,scanLstring));
		}
		if (ids.find(scanRstring) == ids.end()) {
			ids.insert(pair<string,int>(scanRstring,++runningid));
			name.insert(pair<int,string>(runningid,scanRstring));
		}
		// add edge
		add_edge(ids.find(scanLstring)->second, ids.find(scanRstring)->second, *g);
	}
	linkFile.close();
	cerr << "\r" << lico << " gene pairs read" << std::endl;
	
	
	//// SHOW GRAPH
	/*	std::cerr<<std::endl<<"Graph out-edges:"<<std::endl;
	 for (tie(vi, vi_end) = vertices(*g); vi != vi_end; ++vi) {
	 std::cerr<<name[*vi] << " " << *vi <<" outedges - ";
	 for(tie(oei, oei_end)=out_edges(*vi, *g); oei != oei_end; ++oei) {
	 std::cerr<<name[source(*oei, *g)]<<"<->"<<name[target(*oei, *g)]<<"  ";
	 }
	 std::cerr<<std::endl;
	 }
	 std::cerr<<std::endl;
	 */	//// SHOW GRAPH
	
	int numComponents = connected_components(*g, get(vertex_component, *g));
	
	std::cerr<<"Graph has "<<numComponents<<" components."<<std::endl;


	string familyprefix = "GF";

	int connex = 0;
	int compo = numComponents;
	int clusterNumber = 0;
	for(int i=0; i<numComponents; i++) {
		if (--compo % 100 == 0) std::cerr << "\r" << compo << " ";
		keep_all efilter;
		vertexComponent<ComponentMap> vfilter(get(vertex_component, *g), i);
		fg = new FilteredGraph(*g, efilter, vfilter);
		
		
		//std::cerr<<"Filtered graph (component "<<i<<") ";

		setstring Genes;
		int elements = 0;
		for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
			Genes.insert(name[*fvi]);
			++elements;
		}
		// printing
		++connex;
		if (elements > 1 && elements == Genes.size()) {
			++clusterNumber;
			for (setstring::const_iterator it = Genes.begin(); it != Genes.end(); ++it) {
				std::cout << familyprefix << clusterNumber << "\t" << *it << std::endl;
			}
		}
		else std::cerr << "D'oh!" << endl;

		delete fg;
		fg = NULL;
	}
	
	std::cerr << "\rFound " << clusterNumber << " (" << connex << ") connected components (gene families)."<<std::endl;
	
	delete g;
	g = NULL;
	
	return EXIT_SUCCESS;
	
	
	
} // END OF MAIN

