#ifndef __HIDDENMARKOVMODEL_HPP
#define __HIDDENMARKOVMODEL_HPP

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include "constants.hpp"
#include "state.hpp"
#include "graph.hpp"
#include "utils.hpp"

template<typename Elem>
using Matrix = std::vector<std::vector<Elem>>;

/* <-------- Exceptions --------> */

class HMMException : public std::logic_error {
protected:
	HMMException(const std::string& message) :
		std::logic_error(message) {}
};

class StateNotFoundException : public HMMException {
public:
	template<typename T>
	StateNotFoundException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("StateNotFoundException: " + msg, t)) {}

	StateNotFoundException(const std::string& msg) : 
		HMMException("StateNotFoundException: " + msg) {}
};

class StateExistsException : public HMMException {
public:
	template<typename T>
	StateExistsException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("StateExistsException: " + msg, t)) {}
};

class TransitionNotFoundException : public HMMException {
public:
	template<typename T>
	TransitionNotFoundException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("TransitionNotFoundException: " + msg, t)) {}

	TransitionNotFoundException(const std::string& msg) : 
		HMMException("TransitionNotFoundException: " + msg) {}
};

class TransitionExistsException : public HMMException {
public:
	template<typename T>
	TransitionExistsException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("TransitionExistsException: " + msg, t)) {}
};

class TransitionLogicException : public HMMException {
public:
	template<typename T>
	TransitionLogicException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("TransitionLogicException: " + msg, t)) {}
};

/* <----------------------------> */

class HiddenMarkovModel{
private:
	/* Name of this hmm */
	std::string _name;

	/* Begin and end states */
	State* _begin;
	State* _end;

	/* Holds the states and the transitions when building the hmm */
	Graph<State> _graph;

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

	bool has_state(const State& state) const {
		return _graph.has_vertex(state);
	}

	bool has_transition(const State& from_state, const State& to_state) const {
		return _graph.has_edge(from_state, to_state);
	}

	State& begin() { 
		if(_begin != nullptr){
			return *_begin;
		}
		else{
			throw StateNotFoundException(error_message::kHMMHasNoBeginState);
		}
	}

	State& end() {
		if(_end != nullptr){
			return *_end;
		}
		else{
			throw StateNotFoundException(error_message::kHMMHasNoEndState);
		}
	}

	/* See behavior of Graph::add_vertex() */	
	void add_state(const State& state){
		try{
			_graph.add_vertex(state);	
		}
		catch(const VertexExistsException<State>& e){
			throw StateExistsException(e.trigger(), error_message::kHMMAddStateExists);
		}
	}

	/* See behavior of Graph::remove_vertex() */
	void remove_state(const State& state){
		if(state == *_begin) _begin = nullptr;
		else if(state == *_end) _end = nullptr;
		try{
			_graph.remove_vertex(state);	
		}
		catch(const VertexNotFoundException<State>& e){
			throw StateNotFoundException(e.trigger(), error_message::kHMMRemoveStateNotFound);
		}
	}

	std::string transition_string(const State& from, const State& to) const {
		return from.to_string() + " -> " + to.to_string();
	}

	/* See behavior of Graph::add_edge() */
	void add_transition(const State& from, const State& to, double probability){
		if(from == end()) throw TransitionLogicException(transition_string(from, to), error_message::kAddedTransitionFromEndState);
		if(to == begin()) throw TransitionLogicException(transition_string(from, to), error_message::kAddedTransitionToBeginState);
		if(probability < 0) throw TransitionLogicException(transition_string(from, to), error_message::kAddedTransitionNegativeProbability);
		try{
			_graph.add_edge(from, to, probability);	
		} 
		catch(const EdgeExistsException<Edge<State>>& e){
			throw TransitionExistsException(transition_string(*(e.trigger().from()), *(e.trigger().to())), error_message::kHMMAddTransitionExists);
		}
		catch(const IncidentVertexNotFoundException<State>& e){
			throw StateNotFoundException(e.trigger(), error_message::kAddTransitionStateNotFound);
		}
	}

	void begin_transition(const State& state, double probability) {
		add_transition(begin(), state, probability);
	}

	void end_transition(const State& state, double probability) {
		add_transition(state, end(), probability);
	}

	/* See behavior of Graph::remove_edge() */
	void remove_transition(const State& from, const State& to){
		try{
			_graph.remove_edge(from, to);	
		}
		catch(const EdgeNotFoundException<Edge<State>>& e){
			throw TransitionNotFoundException(transition_string(*(e.trigger().from()), *(e.trigger().to())), error_message::kHMMRemoveTransitionNotFound);
		}
		
	}

	/* Prepares the hmm before calling algorithms on it. */
	void brew() {
		/* Get the states from graph. */
		std::vector<State*> states = _graph.get_vertices();
		/* Keep track of the matrix index of each state. */
		std::map<const State*, std::size_t> states_indices;

		/* Raw transition matrix. */
		Matrix<double> A(states.size());
		/* Init one row per state and and map state to matrix index. */
		for(std::size_t i = 0; i < states.size(); ++i) {
			A[i] = std::vector<double>(states.size());
			states_indices[states[i]] = i;
		}
		/* Fill transition matrix with log probabilities. */
		for(std::size_t i = 0; i < states.size(); ++i) {
			std::vector<Edge<State>*> out_edges = _graph.get_out_edges(*(states[i]));
			for(Edge<State>* out_edge : out_edges){
				A[i][states_indices[out_edge->from()]] = (out_edge->weight() == nullptr) ? utils::kNegInf : log(*(out_edge->weight()));
			}
			/* Sum the probabilities to check wether a normalization is needed. */
			double prob_sum = std::accumulate(out_edges.begin(), out_edges.end(), double(0), 
				[](const double previous, const Edge<State>* edge) { 
					return (edge->weight() == nullptr) ? previous : previous + *(edge->weight());
				});
			/* Normalize with log probabilities. */
			if(prob_sum != 1.0){
				utils::log_normalize(A[i].begin(), A[i].end(), log(prob_sum));
			}
		}
		/* Fill emission matrix with PDFs. */
		std::vector<Distribution*> B;


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