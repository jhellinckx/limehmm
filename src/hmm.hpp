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
#include <utility>
#include "constants.hpp"
#include "state.hpp"
#include "graph.hpp"
#include "utils.hpp"

template<typename Elem>
using Matrix = std::vector<std::vector<Elem>>;

std::ostream& operator<<(std::ostream& out, const Matrix<double>& matrix){
	std::size_t longest_string = 0;
	for(std::size_t i = 0; i < matrix.size(); ++i){
		for(std::size_t j = 0; j < matrix[i].size(); ++j){
			std::string double_string = std::to_string(matrix[i][j]);
			if(double_string.length() > longest_string) longest_string = double_string.length();
		}
	}
	for(std::size_t i = 0; i < matrix.size(); ++i){
		for(std::size_t j = 0; j < matrix[i].size(); ++j){
			std::string double_string = std::to_string(matrix[i][j]);
			out << std::string(longest_string - double_string.length(), ' ');
			out << double_string;
			out << ' ';
		}
		out << std::endl;
	}
	return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<Distribution*>& vec){
	for(const Distribution* dist : vec){
		if(dist == nullptr) out << "Silent" << std::endl;
		else out << *dist << std::endl;
	}
	return out;
}

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
	/* Name of this hmm. */
	std::string _name;

	/* Begin and end states. */
	State* _begin;
	State* _end;

	/* Holds the states and the transitions when building the hmm. */
	Graph<State> _graph;

	/* Generated by brew. */
	std::map<std::string, std::size_t> _states_indices;
	std::size_t _begin_state_index;
	std::size_t _end_state_index;
	Matrix<double> _A;
	std::vector<Distribution*> _B;
	std::vector<double> _pi_begin;
	std::vector<double> _pi_end;

	void _clear_raw_data() {
		for(Distribution* dist : _B){
			if(dist != nullptr) delete dist;
		}
		_states_indices.clear();
		_pi_begin.clear();
		_pi_end.clear();
		_A.clear();
	}

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
		/* Get rid of previous data. */
		_clear_raw_data();
		/* Get the states from graph. */
		std::vector<State*> states = _graph.get_vertices();		
		/* Keep track of the matrix index of each state. */
		std::map<std::string, std::size_t> states_indices;
		/* Raw transition matrix. Init its size. Decrement by 2 since
		begin and end states transitions are stored in separated arrays. */
		Matrix<double> A(states.size() - 2);
		/* Begin/end states transitions. */
		std::vector<double> pi_begin(states.size() - 2);
		std::vector<double> pi_end(states.size() - 2);
		/* Init one row per state and map state to matrix index. 
		Default values set to negative infinity since we use log probabilites. */
		std::size_t i = 0;
		for(State* p_state : states){
			State& state = *p_state;
			if(state != begin() && state != end()){
				if(state.is_silent()) throw std::runtime_error("silent non-begin/end states not yet supported");
				A[i] = std::vector<double>(states.size() - 2, utils::kNegInf);
				states_indices[state.name()] = i;
				++i;
			}
		}
		/* Fill transitions with log probabilities and check whether a normalization is needed. */
		auto fill_normalize = [&states_indices] (const std::vector<Edge<State>*>& edges, std::vector<double>& prob_vec_to_fill, bool update_from) {
			std::function<std::size_t(const Edge<State>*)> state_index;
			if(update_from) state_index = [&states_indices](const Edge<State>* edge){ return states_indices[edge->from()->name()]; };
			else state_index = [&states_indices](const Edge<State>* edge){ return states_indices[edge->to()->name()]; };
			double prob_sum = 0;
			double prob;
			for(Edge<State>* edge : edges){
				prob = (edge->weight() == nullptr) ? 0 : *(edge->weight());
				prob_sum += prob;
				prob_vec_to_fill[state_index(edge)] = log(prob);
			}
			if(prob_sum != 1.0){
				utils::for_each_log_normalize(prob_vec_to_fill.begin(), prob_vec_to_fill.end(), log(prob_sum));
			}
		};
		for(std::size_t i = 0; i < states.size(); ++i){
			State& state = *states[i];
			/* Fill begin transitions. Throw exception if begin state has predecessors. */
			if(state == begin()){
				if(_graph.get_in_edges(state).size() > 0) throw std::logic_error("begin state cannot have predecessors");
				std::vector<Edge<State>*> out_edges = _graph.get_out_edges(state);
				fill_normalize(out_edges, pi_begin, false);
			}
			
			/* Fill end transitions. Throw exception if end state has successors. */
			else if(state == end()){
				if(_graph.get_out_edges(state).size() > 0) throw std::logic_error("end state cannot have successors");
				std::vector<Edge<State>*> in_edges = _graph.get_in_edges(state);
				fill_normalize(in_edges, pi_end, true);
			}
			/* Fill normal transitions aka matrix A. */
			else{
				std::vector<Edge<State>*> out_edges = _graph.get_out_edges(state);
				fill_normalize(out_edges, A[states_indices[state.name()]], false);
			}
		}
		/* Fill emission matrix with the states PDFs. */
		std::vector<Distribution*> B(states.size() - 2);
		for(const State* p_state : states){
			const State& state = *p_state;
			if(state != begin() && state != end()){
				if(state.is_silent()) throw std::runtime_error("silent non-begin/end states not yet supported");
				Distribution* distribution = state.distribution().clone();
				distribution->log_normalize();
				B[states_indices[state.name()]] = distribution;
			}
		}

		_A = std::move(A);
		_B = std::move(B);
		_states_indices = std::move(states_indices);
		_pi_begin = pi_begin;
		_pi_end = pi_end;

	}

	Matrix<double>& raw_transitions() { return _A; }
	std::vector<Distribution*>& raw_pdfs() { return _B; }
	std::vector<double>& raw_pi_begin() { return _pi_begin; }
	std::vector<double>& raw_pi_end() { return _pi_end; }
	std::map<std::string, std::size_t>& states_indices() { return _states_indices; }

	// template<typename Symbol> //TODO
	// std::vector<double> forward(const std::vector<Symbol>& symbols) {
	// 	std::pair<std::vector<double>, std::vector<double>> fwd = std::make_pair(std::vector<double>(_A.size()), std::vector<double>(_A.size()));
	// 	fwd.first;
	// }

	// template<typename Symbol>
	// std::vector<double> backward(const std::vector<Symbol>& symbols) {

	// }

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

	virtual ~HiddenMarkovModel(){
		_clear_raw_data();
	}
};


#endif