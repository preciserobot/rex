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


namespace boost {
	enum vertex_component_t { vertex_component = 111 };
	BOOST_INSTALL_PROPERTY(vertex, component);
}
using namespace boost;

typedef map<string, int> GeneID;
typedef map<string, int> stringcounts;
typedef map<int, string> invGeneID;
typedef map<string,string> stringstring;
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
	setstring prefixes; // for the ordered output
	
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
		prefixes.insert(scanLstring.substr(0,8));
		prefixes.insert(scanRstring.substr(0,8));
		
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
	

	
	setint doneVertex;

	
	
	int numComponents = connected_components(*g, get(vertex_component, *g));
	
	std::cerr<<"Graph has "<<numComponents<<" components."<<std::endl;
	
	
	int yippies  = 0;
	int complete = 0;
	int puzzling = 0;
	int crappy   = 0;
	int runningclusterid = 0;
	int compo = numComponents;
	for(int i=0; i<numComponents; i++) {
		if (--compo % 100 == 0) std::cerr << "\r" << compo << " ";
		keep_all efilter;
		vertexComponent<ComponentMap> vfilter(get(vertex_component, *g), i);
		fg = new FilteredGraph(*g, efilter, vfilter);
		
		//std::cerr<<"Filtered graph (component "<<i<<") ";
		//count members (vertices) and edges
		int subgraphMembers = 0;
		int subgraphEdges = 0;
		for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
			++subgraphMembers;
			for(tie(foei, foei_end)=out_edges(*fvi, *fg); foei != foei_end; ++foei) {
				++subgraphEdges;
			}
		}
		int ll = subgraphMembers;
		if (subgraphMembers * --ll == subgraphEdges){
			// check if no ensembl crappy one2one
			stringstring orderedGenes;
			int elements = 0;
			for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
				orderedGenes.insert(pair<string,string>(name[*fvi].substr(0,8),name[*fvi]));
				++elements;
				doneVertex.insert(*fvi);
			}
			// printing
			++complete;
			if (elements == orderedGenes.size()) {
				++yippies;
				std::cout << ++runningclusterid;
				for (setstring::const_iterator it = prefixes.begin(); it != prefixes.end(); ++it) {
					std::cout << "\t" << orderedGenes[*it];
				}
				std::cout << "\n";
			}
			else std::cerr << "D'oh!" << endl;
		}		
		else {
			// crappy ensembl or incomplete
			++puzzling;
/*			
			
			
			// make vertice name prefix list

			stringcounts prefixes;
			int totalPrefixes = 0;
			for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
				++prefixes[name[*fvi].substr(0,6)];
				//std::cerr << name[*fvi].substr(0,6) << "\t";
				++totalPrefixes;
			}
			//std::cerr << endl;
			
			if (prefixes.size() != totalPrefixes) {
				++crappy;
				// delete anything which has more than one prefix count
				
				for (stringcounts::const_iterator it = prefixes.begin(); it != prefixes.end(); ++it) {
					if (it->second > 1) {
						// multigene
						string killme = it->first;
						// remove vertices with this name
						//std::cerr << it->first << "\t" << it->second << "\n";
						for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
							if (killme.compare(name[*fvi].substr(0,6)) == 0) {
								// equal strings
								doneVertex.insert(*fvi);
							}
							//std::cerr << name[*fvi].substr(0,6) << "\t";
							++totalPrefixes;
						}
					}
				}
				
				
				// check vertices if there are duplicate substrings
				//std::cerr<<std::endl<<"Crappy SubGraph out-edges:"<<std::endl;
				//for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
					//std::cerr<<name[*fvi]<<" outedges - ";
					//for(tie(foei, foei_end)=out_edges(*fvi, *fg); foei != foei_end; ++foei) std::cerr<<name[source(*foei, *fg)]<<"<->"<<name[target(*foei, *fg)]<<"  ";
					//std::cerr<<std::endl;
				//}
				//std::cerr<<std::endl;
			}
 */
		}
		
		delete fg;
		fg = NULL;
	}
	
	//ITERATION STEP
/*
	while (doneVertex.size() > 0) {
		//iterate
		int vnum = doneVertex.size();
		std::cerr << "Killing " << vnum << " vertices" << endl;
		// clear used vertices
		for (setint::const_iterator it = doneVertex.begin(); it != doneVertex.end(); ++it) {
			if (--vnum % 100 == 0) std::cerr << "\r" << vnum;
			clear_vertex(*it, *g);
		}
		vnum = doneVertex.size();
		for (setint::const_iterator it = doneVertex.begin(); it != doneVertex.end(); ++it) {
			if (--vnum % 100 == 0) std::cerr << "\r" << vnum;
			remove_vertex(*it, *g);
		}
		doneVertex.clear();

		// of the remaining connected components list lowest connection number and add to kill list
		std::cerr << "lala";
		numComponents = connected_components(*g, get(vertex_component, *g));
		std::cerr << "lolo";
		std::cerr<<"Graph has "<<numComponents<<" components."<<std::endl;
		compo = numComponents;
		/*
		for(int i=0; i<numComponents; i++) {
			if (--compo % 100 == 0) std::cerr << "\r" << compo << " ";
			keep_all efilter;
			vertexComponent<ComponentMap> vfilter(get(vertex_component, *g), i);
			fg = new FilteredGraph(*g, efilter, vfilter);
			
			//std::cerr<<"Filtered graph (component "<<i<<") ";
			//count members (vertices) and edges
			int subgraphMembers = 0;
			int subgraphEdges = 0;
			for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
				++subgraphMembers;
				for(tie(foei, foei_end)=out_edges(*fvi, *fg); foei != foei_end; ++foei) {
					++subgraphEdges;
				}
			}
			int ll = subgraphMembers;
			if (subgraphMembers * --ll == subgraphEdges){
				// check if no ensembl crappy one2one
				stringstring orderedGenes;
				int elements = 0;
				for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
					orderedGenes.insert(pair<string,string>(name[*fvi].substr(0,8),name[*fvi]));
					++elements;
					doneVertex.insert(*fvi);
				}
				// printing
				++complete;
				if (elements == orderedGenes.size()) {
					++yippies;
					std::cout << ++runningclusterid;
					for (setstring::const_iterator it = prefixes.begin(); it != prefixes.end(); ++it) {
						std::cout << "\t" << orderedGenes[*it];
					}
					std::cout << "\n";
				}
				else std::cerr << "D'oh!" << endl;
			}		
			else {
				// crappy ensembl or incomplete
				++puzzling;
				
				
				
				// make vertice name prefix list
				
				stringcounts prefixes;
				int totalPrefixes = 0;
				for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
					++prefixes[name[*fvi].substr(0,6)];
					//std::cerr << name[*fvi].substr(0,6) << "\t";
					++totalPrefixes;
				}
				//std::cerr << endl;
				
				if (prefixes.size() != totalPrefixes) {
					++crappy;
					// delete anything which has more than one prefix count
					
					for (stringcounts::const_iterator it = prefixes.begin(); it != prefixes.end(); ++it) {
						if (it->second > 1) {
							// multigene
							string killme = it->first;
							// remove vertices with this name
							//std::cerr << it->first << "\t" << it->second << "\n";
							for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
								if (killme.compare(name[*fvi].substr(0,6)) == 0) {
									// equal strings
									doneVertex.insert(*fvi);
								}
								//std::cerr << name[*fvi].substr(0,6) << "\t";
								++totalPrefixes;
							}
						}
					}
					
					
					// check vertices if there are duplicate substrings
					//std::cerr<<std::endl<<"Crappy SubGraph out-edges:"<<std::endl;
					//for (tie(fvi, fvi_end) = vertices(*fg); fvi != fvi_end; ++fvi) {
					//std::cerr<<name[*fvi]<<" outedges - ";
					//for(tie(foei, foei_end)=out_edges(*fvi, *fg); foei != foei_end; ++foei) std::cerr<<name[source(*foei, *fg)]<<"<->"<<name[target(*foei, *fg)]<<"  ";
					//std::cerr<<std::endl;
					//}
					//std::cerr<<std::endl;
				}
			}
			
			delete fg;
			fg = NULL;
		}
	}
	
	*/
	
	
	
	
	std::cerr << "\rFound " << complete << " complete subgraphs."<<std::endl;
	std::cerr << yippies << " YIPPIES shouted (means no ensembl error)" << endl;
	std::cerr << puzzling << " INCOMPLETE or PUZZLING subgraphs " << endl;
	std::cerr << crappy << " CRAPPY subgraphs (unnamed ensembl paralogies)" << endl;
	
	delete g;
	g = NULL;
	
	return EXIT_SUCCESS;
	
	
		
} // END OF MAIN

