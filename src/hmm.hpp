#ifndef __HIDDENMARKOVMODEL_HPP
#define __HIDDENMARKOVMODEL_HPP

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include "constants.hpp"
#include "state.hpp"
#include "graph.hpp"

template<typename T>
using Matrix = std::vector<std::vector<T>>;

class HiddenMarkovModel{
private:
	/* Name of this hmm */
	std::string _name;

	/* Begin and end states */
	State* _begin;
	State* _end;

	/* Holds the states and the transitions when building the hmm */
	Graph<State> _graph;

	/* Fields initialized after brew(). 
	These fields are mainly used by hmm algorithms. */
	/* Transition matrix */
	Matrix<double> _A;
	/* Distributions set */
	std::vector<Distribution> _B;

public:
	/* Default constructor. Inits an empty hmm with default values. */
	HiddenMarkovModel() : 
		HiddenMarkovModel(std::to_string((ptrdiff_t)this)) {}

	HiddenMarkovModel(const std::string& name) : 
		HiddenMarkovModel(name, State(hmm_config::kDefaultStartStateLabel + name), 
			State(hmm_config::kDefaultEndStateLabel + name)) {}

	HiddenMarkovModel(const State& begin, const State& end) :
		HiddenMarkovModel(std::to_string((ptrdiff_t)this), begin, end) {}

	/* Complete constructor. */
	HiddenMarkovModel(const std::string& name, const State& begin, const State& end) : 
		_name(name), _begin(nullptr), _end(nullptr), _graph() {
			_graph.add_vertex(begin);
			_begin = _graph.get_vertex(begin);
			_graph.add_vertex(end);
			_end = _graph.get_vertex(end);
	}

	std::string name() const { return _name; }
	void set_name(const std::string& name) { _name = name; } 
	std::size_t num_states() const { return _graph.num_vertices(); }
	std::size_t num_transitions() const { return _graph.num_edges(); }	

	bool contains(const State& state) const {
		return _graph.contains(state);
	}

	State& begin() { 
		if(_begin != nullptr){
			return *_begin;
		}
		else{
			throw std::logic_error("No beginning");
		}
	}

	State& end() {
		if(_end != nullptr){
			return *_end;
		}
		else{
			throw std::logic_error("No end");
		}
	}

	/* See behavior of Graph::add_vertex() */	
	void add_state(const State& state){
		_graph.add_vertex(state);
	}

	/* See behavior of Graph::remove_vertex() */
	void remove_state(const State& state){
		if(state == *_begin) _begin = nullptr;
		else if(state == *_end) _end = nullptr;
		_graph.remove_vertex(state);
	}

	/* See behavior of Graph::add_edge() */
	void add_transition(const State& from, const State& to, double probability){
		_graph.add_edge(from, to, probability);
	}

	/* See behavior of Graph::remove_edge() */
	void remove_transition(const State& from, const State& to){
		_graph.remove_edge(from, to);
	}

	/* Prepares the hmm before calling algorithms on it */
	void brew() {
		// /* Init raw emission and transition matrices */
		// std::vector<State*>& states = _graph.get_vertices();
		// std::for_each(states.begin(), states.end(),
		// 				[&_A, &_B]
		// /* Normalize transitions */
		// std::vector<State*>& states = _graph.get_vertices();
		// std::for_each(states.begin(), states.end(),
		// 				[_A, _B]

		// /* Normalize emissions */
	}

	void sample() {

	}

	void decode() {

	}

	void train() {

	}

	void save(){

	}

	void load(){

	}
};


#endif