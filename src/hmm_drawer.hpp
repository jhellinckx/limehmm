#ifndef __GRAPH_DRAWER_HPP
#define __GRAPH_DRAWER_HPP

#include "hmm.hpp"
#include "graph.hpp"

namespace graph_drawer {
	extern void draw_to_file(const Graph& graph, std::string filename){

	}
	extern void draw_to_file(const HiddenMarkovModel& hmm, std::string filename){
		draw_to_file(hmm.get_graph(), filename);
	}
}

#endif